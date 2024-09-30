
# Import required modules
using Oceananigans
using Oceananigans.Units: minutes, hour, hours, day
using CairoMakie
using Printf

# Set up the grid
grid = RectilinearGrid(GPU(), size=(64, 64), extent=(64, 64), halo=(3, 3), topology=(Periodic, Flat, Bounded))

# Boundary conditions
buoyancy_flux(x, t, params) = params.initial_buoyancy_flux * exp(-t^4 / (24 * params.shut_off_time^4))
buoyancy_flux_parameters = (initial_buoyancy_flux = 1e-8, shut_off_time = 2hours)
buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, parameters = buoyancy_flux_parameters)

# Plot buoyancy flux over time
times = range(0, 12hours, length=100)
fig = Figure(size = (800, 300))
ax = Axis(fig[1, 1]; xlabel = "Time (hours)", ylabel = "Surface buoyancy flux (m² s⁻³)")
flux_time_series = [buoyancy_flux(0, t, buoyancy_flux_parameters) for t in times]
lines!(ax, times ./ hour, flux_time_series)
fig

# Define boundary conditions
N² = 1e-4 # s⁻²
buoyancy_gradient_bc = GradientBoundaryCondition(N²)
buoyancy_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc, bottom = buoyancy_gradient_bc)

# Phytoplankton dynamics
growing_and_grazing(x, z, t, P, params) = (params.μ₀ * exp(z / params.λ) - params.m) * P
plankton_dynamics_parameters = (μ₀ = 1/day, λ = 5, m = 0.1/day)
plankton_dynamics = Forcing(growing_and_grazing, field_dependencies = :P, parameters = plankton_dynamics_parameters)

# Set up the model
model = NonhydrostaticModel(
    grid = grid,
    advection = UpwindBiasedFifthOrder(),
    timestepper = :RungeKutta3,
    closure = ScalarDiffusivity(ν=1e-4, κ=1e-4),
    coriolis = FPlane(f=1e-4),
    tracers = (:b, :P),
    buoyancy = BuoyancyTracer(),
    forcing = (; P = plankton_dynamics),
    boundary_conditions = (; b = buoyancy_bcs)
)

# Initial conditions
mixed_layer_depth = 32 # m
stratification(z) = z < -mixed_layer_depth ? N² * z : -N² * mixed_layer_depth
noise(z) = 1e-4 * N² * grid.Lz * randn() * exp(z / 4)
initial_buoyancy(x, z) = stratification(z) + noise(z)
set!(model, b = initial_buoyancy, P = 1)

# Set up the simulation
simulation = Simulation(model, Δt = 2minutes, stop_time = 24hours)

# Adaptive time stepping
conjure_time_step_wizard!(simulation, cfl = 1.0, max_Δt = 2minutes)

# Progress callback
function progress(sim)
    @printf("Iteration: %d, time: %s, Δt: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt))
end
add_callback!(simulation, progress, IterationInterval(100))

# Output writer
outputs = (w = model.velocities.w, P = model.tracers.P, avg_P = Average(model.tracers.P, dims=(1, 2)))
simulation.output_writers[:simple_output] = JLD2OutputWriter(model, outputs, schedule = TimeInterval(20minutes), filename = "convecting_plankton.jld2", overwrite_existing = true)

# Run the simulation
run!(simulation)

# Visualization setup
filepath = simulation.output_writers[:simple_output].filepath
w_timeseries = FieldTimeSeries(filepath, "w")
P_timeseries = FieldTimeSeries(filepath, "P")
avg_P_timeseries = FieldTimeSeries(filepath, "avg_P")
buoyancy_flux_time_series = [buoyancy_flux(0, t, buoyancy_flux_parameters) for t in w_timeseries.times]

# Animation setup
n = Observable(1)
title = @lift @sprintf("t = %s", prettytime(w_timeseries.times[$n]))
wn = @lift w_timeseries[$n]
Pn = @lift P_timeseries[$n]
avg_Pn = @lift avg_P_timeseries[$n]
w_lim = maximum(abs, interior(w_timeseries))
P_lims = (0.95, 1.1)

fig = Figure(size = (1200, 1000))
ax_w = Axis(fig[2, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_P = Axis(fig[3, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_b = Axis(fig[2, 3]; xlabel = "Time (hours)", ylabel = "Buoyancy flux (m² s⁻³)", yaxisposition = :right)
ax_avg_P = Axis(fig[3, 3]; xlabel = "Plankton concentration (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_P, 0.85, 1.3)
fig[1, 1:3] = Label(fig, title, tellwidth=false)

hm_w = heatmap!(ax_w, wn; colormap = :balance, colorrange = (-w_lim, w_lim))
Colorbar(fig[2, 1], hm_w; label = "Vertical velocity (m s⁻¹)", flipaxis = false)
hm_P = heatmap!(ax_P, Pn; colormap = :matter, colorrange = P_lims)
Colorbar(fig[3, 1], hm_P; label = "Plankton 'concentration'", flipaxis = false)
lines!(ax_b, w_timeseries.times ./ hour, buoyancy_flux_time_series; linewidth = 1, color = :black, alpha = 0.4)
scatter!(ax_b, @lift Point2(w_timeseries.times[$n] / hour, buoyancy_flux_time_series[$n]); marker = :circle, markersize = 16, color = :black)
lines!(ax_avg_P, avg_Pn)

# Record the animation
frames = 1:length(w_timeseries.times)
record(fig, "convecting_plankton.mp4", frames, framerate=8) do i
    n[] = i
end
