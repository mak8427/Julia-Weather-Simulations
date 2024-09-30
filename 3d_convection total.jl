# Import necessary packages
using Pkg
pkg"add Oceananigans, CairoMakie"  # Ensure dependencies are installed

# Import required modules
using Oceananigans
using Oceananigans.Units: minutes, hour, hours, day
using CairoMakie  # For visualization

# -----------------------------------
# GRID SETUP
# -----------------------------------
# Create a 2D grid with 64 points in x and z directions, flat in y direction.
grid = RectilinearGrid(size=(64, 64), extent=(64, 64), halo=(3, 3), topology=(Periodic, Flat, Bounded))

# The grid is 64x1x64 with 3 halo points for high-order advection accuracy.
# x is periodic, z is bounded, and y is flat.
println(grid)

# -----------------------------------
# BOUNDARY CONDITIONS
# -----------------------------------
# Define time-dependent surface buoyancy flux function that decays over time.
buoyancy_flux(x, t, params) = params.initial_buoyancy_flux * exp(-t^4 / (24 * params.shut_off_time^4))
buoyancy_flux_parameters = (initial_buoyancy_flux = 1e-8, shut_off_time = 2hours)

# Apply the flux boundary condition at the top of the domain.
buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, parameters = buoyancy_flux_parameters)

# Set the buoyancy gradient at the bottom (stabilizing condition).
N² = 1e-4  # Buoyancy frequency squared (stabilizing).
buoyancy_gradient_bc = GradientBoundaryCondition(N²)

# Combine the top and bottom boundary conditions into one set.
buoyancy_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc, bottom = buoyancy_gradient_bc)

# -----------------------------------
# PHOTOPLANKTON DYNAMICS
# -----------------------------------
# Define phytoplankton growth as a function of light and grazing by zooplankton.
growing_and_grazing(x, z, t, P, params) = (params.μ₀ * exp(z / params.λ) - params.m) * P
plankton_dynamics_parameters = (μ₀ = 1/day, λ = 5, m = 0.1/day)  # Growth and decay rates

# Create a Forcing object to represent phytoplankton dynamics.
plankton_dynamics = Forcing(growing_and_grazing, field_dependencies = :P, parameters = plankton_dynamics_parameters)

# -----------------------------------
# MODEL SETUP
# -----------------------------------
# Build the nonhydrostatic model for planktonic convection with the given grid and boundary conditions.
model = NonhydrostaticModel(; 
    grid, 
    advection = UpwindBiasedFifthOrder(),
    timestepper = :RungeKutta3,
    closure = ScalarDiffusivity(ν=1e-4, κ=1e-4),  # Viscosity and diffusivity
    coriolis = FPlane(f=1e-4),  # Coriolis parameter
    tracers = (:b, :P),  # Tracers: buoyancy (b) and phytoplankton (P)
    buoyancy = BuoyancyTracer(),
    forcing = (; P=plankton_dynamics),  # Apply phytoplankton forcing
    boundary_conditions = (; b=buoyancy_bcs)  # Apply buoyancy boundary conditions
)

# -----------------------------------
# INITIAL CONDITIONS
# -----------------------------------
# Define initial conditions for buoyancy (b) and phytoplankton (P).
mixed_layer_depth = 32  # Depth of the mixed layer in meters
stratification(z) = z < -mixed_layer_depth ? N² * z : - N² * mixed_layer_depth
noise(z) = 1e-4 * N² * grid.Lz * randn() * exp(z / 4)  # Random noise for initial buoyancy
initial_buoyancy(x, z) = stratification(z) + noise(z)

# Initialize buoyancy and phytoplankton concentration in the model.
set!(model, b=initial_buoyancy, P=1)

# -----------------------------------
# SIMULATION SETUP
# -----------------------------------
# Configure the simulation with a time-step of 2 minutes and a stop time of 24 hours.
simulation = Simulation(model, Δt=2minutes, stop_time=24hours)

# Adapt the time-step using a CFL number of 1.0, with a max time-step of 2 minutes.
conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=2minutes)

# Add a callback to print progress during the simulation.
using Printf
progress(sim) = @printf("Iteration: %d, time: %s, Δt: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt))
add_callback!(simulation, progress, IterationInterval(100))

# -----------------------------------
# OUTPUT WRITING
# -----------------------------------
# Set up an output writer to log velocities and plankton concentration every 20 minutes.
outputs = (w = model.velocities.w, P = model.tracers.P, avg_P = Average(model.tracers.P, dims=(1, 2)))
simulation.output_writers[:simple_output] = JLD2OutputWriter(
    model, outputs, 
    schedule = TimeInterval(20minutes), 
    filename = "convecting_plankton.jld2", 
    overwrite_existing = true
)

# -----------------------------------
# RUN THE SIMULATION
# -----------------------------------
run!(simulation)

# -----------------------------------
# VISUALIZATION (OPTIONAL)
# -----------------------------------
# To visualize, load the output file and create a time-series of the results.
# Extract the w and P fields as matrices to be used in heatmap!
w_field = interior(w_timeseries[1])  # Extract the interior part of the field (2D data)
P_field = interior(P_timeseries[1])  # Same for phytoplankton concentration

# Get the grid points for the x and z coordinates using the appropriate accessor functions
x_coords = model.grid.x  # x-coordinates of the grid points
z_coords = model.grid.z  # z-coordinates of the grid points

# Plotting with CairoMakie
fig = Figure(size = (1200, 1000))

# Vertical velocity heatmap
ax_w = Axis(fig[2, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, x_coords, z_coords, w_field; colormap = :balance, colorrange = (-maximum(abs, w_field), maximum(abs, w_field)))
Colorbar(fig[2, 1], hm_w; label = "Vertical velocity (m/s)", flipaxis = false)

# Phytoplankton concentration heatmap
ax_P = Axis(fig[3, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_P = heatmap!(ax_P, x_coords, z_coords, P_field; colormap = :matter, colorrange = (0.95, 1.1))
Colorbar(fig[3, 1], hm_P; label = "Phytoplankton concentration", flipaxis = false)

# Show the figure
fig