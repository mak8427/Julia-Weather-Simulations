
using CairoMakie
using Oceananigans

grid = RectilinearGrid(size=128, z=(-0.5, 0.5), topology=(Flat, Flat, Bounded))

ScalarDiffusivity{}(ν=0.0, κ=1.0)

closure = ScalarDiffusivity(κ=1)
model = NonhydrostaticModel(; grid, closure, tracers=:T)

width = 0.1
initial_temperature(z) = exp(-z^2 / (2width^2))
set!(model, T=initial_temperature)


using CairoMakie
set_theme!(Theme(fontsize = 24, linewidth=3))

fig = Figure()
axis = (xlabel = "Temperature (ᵒC)", ylabel = "z")
label = "t = 0"

z = znodes(model.tracers.T)
T = interior(model.tracers.T, 1, 1, :)

lines(T, z; label, axis)


min_Δz = minimum_zspacing(model.grid)
diffusion_time_scale = min_Δz^2 / model.closure.κ.T

simulation = Simulation(model, Δt = 0.1 * diffusion_time_scale, stop_iteration = 1000)
run!(simulation)


using Printf

label = @sprintf("t = %.3f", model.clock.time)
lines!(interior(model.tracers.T, 1, 1, :), z; label)
axislegend()
fig