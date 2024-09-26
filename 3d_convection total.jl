# ================================================
# 3D Turbulence Visualization with Oceananigans and GLMakie
# ================================================

# ---------------------
# Step 1: Import Necessary Packages
# ---------------------
using Oceananigans      # For handling simulation data
using GLMakie           # For interactive 3D visualization
using JLD2              # For loading JLD2 files
using Statistics        # For data manipulation
using ColorSchemes      # For colormaps
using Colors            # For Colorant types and macros

# ---------------------
# Step 2: Load Simulation Data
# ---------------------

# Path to your JLD2 simulation output file
jld2_filepath = "3d_turbulence.jld2"

# Load the FieldTimeSeries for ω_mag and E
ω_mag_timeseries = FieldTimeSeries(jld2_filepath, "ω_mag")
E_timeseries = FieldTimeSeries(jld2_filepath, "E")

# Access the times at which data was saved
times = ω_mag_timeseries.times

println("Loaded data at times:")
println(times)

# ---------------------
# Step 3: Select Time Index and Extract Data
# ---------------------

# Choose the time index you want to visualize
# For example, the last time step
time_index = length(times)

# Extract the 3D vorticity magnitude field at the selected time
ω_mag_field = ω_mag_timeseries[time_index]

# Obtain grid coordinates from the field
xC, yC, zC = nodes(ω_mag_field)

# Convert OffsetArrays to standard arrays for compatibility
xC = collect(xC)
yC = collect(yC)
zC = collect(zC)

# Extract the 3D vorticity magnitude data and convert to a standard array
ω_mag_data = collect(ω_mag_field[:, :, :])

println("Grid sizes:")
println("xC: ", length(xC))
println("yC: ", length(yC))
println("zC: ", length(zC))
println("ω_mag_data size: ", size(ω_mag_data))

# ---------------------
# Step 4: Downsample the Data (Optional for Performance)
# ---------------------

# Define downsampling factor
# Adjust this based on your system's performance and desired resolution
downsample_factor = 2  # Use 1 for full resolution, 2 for half, etc.

# Downsample grid coordinates
xC_ds = xC[1:downsample_factor:end]
yC_ds = yC[1:downsample_factor:end]
zC_ds = zC[1:downsample_factor:end]

# Downsample vorticity data
ω_mag_ds = ω_mag_data[1:downsample_factor:end, 1:downsample_factor:end, 1:downsample_factor:end]

println("Downsampled grid sizes:")
println("xC_ds: ", length(xC_ds))
println("yC_ds: ", length(yC_ds))
println("zC_ds: ", length(zC_ds))
println("ω_mag_ds size: ", size(ω_mag_ds))

# ---------------------
# Step 5: Normalize the Data
# ---------------------

# Normalize the data for better visualization (scale between 0 and 1)
ω_mag_normalized = ω_mag_ds ./ maximum(ω_mag_ds)

# ---------------------
# Step 6: Visualization Functions
# ---------------------

# ---------------------
# Option A: Volume Rendering
# ---------------------
# Function for volume rendering
function volume_rendering(x, y, z, data, current_time)
    # Create a figure and 3D axis
    fig = Figure(resolution = (800, 600))
    ax = Axis3(fig[1, 1],
               xlabel = "x",
               ylabel = "y",
               zlabel = "z",
               title = "Vorticity Magnitude Volume at time = $(round(current_time, digits=2))")
    
    # Create the volume plot
    volume!(ax, x, y, z, data;
            algorithm = :absorption,      # Rendering algorithm
            colormap = :inferno,          # Color scheme
            transparency = true)          # Enable transparency
    
    # Adjust the camera for a better view
    ax.scene.camera.projection = :perspective
    cam3d!(ax.scene, eyeposition = Vec3f(2, 2, 2))
    
    # Display the figure
    display(fig)
    
    # Save the figure as an image
    save("vorticity_magnitude_volume_time_$(round(current_time, digits=2)).png", fig)
end


# ---------------------
# Option B: Isosurface Plotting
# ---------------------
function isosurface_plotting(x, y, z, data, current_time)
    # Define isovalue levels (adjust as needed)
    isovalue_levels = [0.3, 0.5, 0.7]
    
    # Create a figure and 3D axis
    fig = Figure(resolution = (800, 600))
    ax = Axis3(fig[1, 1],
               xlabel = "x",
               ylabel = "y",
               zlabel = "z",
               title = "Vorticity Magnitude Isosurfaces at time = $(round(current_time, digits=2))")
    
    # Plot isosurfaces
    for (i, level) in enumerate(isovalue_levels)
        # Ensure the index is within the range of the colormap
        colormap_length = length(ColorSchemes.viridis)
        color_index = clamp(Int(round(Float64(i) * (colormap_length / (length(isovalue_levels) + 1)))), 1, colormap_length)
        
        isosurface!(ax, x, y, z, data;
                    isovalue = level,
                    color = ColorSchemes.viridis[color_index],
                    opacity = 0.5)
    end
    
    # Add lighting for better depth perception
    light!(ax.scene, position = Vec3f(1, 1, 1), color = colorant"#FF0000", intensity = 1.0)
    
    # Adjust the camera for a better view
    ax.camera.projection = :perspective
    cam3d!(ax.scene, eyeposition = Vec3f(2, 2, 2))
    
    # Display the figure
    display(fig)
    
    # Save the figure as an image
    save("vorticity_magnitude_isosurfaces_time_$(round(current_time, digits=2)).png", fig)
end

# ---------------------
# Option C: Slice Plotting
# ---------------------
function slice_plotting(x, y, z, data, current_time)
    # Create a figure and 3D axis
    fig = Figure(resolution = (800, 600))
    ax = Axis3(fig[1, 1],
               xlabel = "x",
               ylabel = "y",
               zlabel = "z",
               title = "Vorticity Magnitude Slices at time = $(round(current_time, digits=2))")
    
    # Define slice positions (central slices)
    x_slice = Int(length(x) ÷ 2)
    y_slice = Int(length(y) ÷ 2)
    z_slice = Int(length(z) ÷ 2)
    
    # Extract slices
    slice_x = data[x_slice, :, :]
    slice_y = data[:, y_slice, :]
    slice_z = data[:, :, z_slice]
    
    # Normalize slices (already normalized globally, but ensure no division by zero)
    slice_x .= slice_x ./ maximum(slice_x)
    slice_y .= slice_y ./ maximum(slice_y)
    slice_z .= slice_z ./ maximum(slice_z)
    
    # Create the slices in the 3D plot
    heatmap!(ax, y, z, slice_x';
             transform = (:x, x[x_slice]),
             colormap = :inferno,
             transparency = true)
    heatmap!(ax, x, z, slice_y';
             transform = (:y, y[y_slice]),
             colormap = :inferno,
             transparency = true)
    heatmap!(ax, x, y, slice_z';
             transform = (:z, z[z_slice]),
             colormap = :inferno,
             transparency = true)
    
    # Add lighting for better depth perception
    light!(ax.scene, position = Vec3f(1, 1, 1), color = colorant"red", intensity = 1.0)
    
    # Adjust the camera for a better view
    ax.camera.projection = :perspective
    cam3d!(ax.scene, eyeposition = Vec3f(2, 2, 2))
    
    # Display the figure
    display(fig)
    
    # Save the figure as an image
    save("vorticity_magnitude_slices_time_$(round(current_time, digits=2)).png", fig)
end

# ---------------------
# Step 7: Run Visualization Options
# ---------------------

# Choose which visualization option to run by uncommenting the desired function call

# Option A: Volume Rendering
volume_rendering(xC_ds, yC_ds, zC_ds, ω_mag_normalized, times[time_index])

# Option B: Isosurface Plotting
# isosurface_plotting(xC_ds, yC_ds, zC_ds, ω_mag_normalized, times[time_index])

# Option C: Slice Plotting
# slice_plotting(xC_ds, yC_ds, zC_ds, ω_mag_normalized, times[time_index])

# ---------------------
# Step 8: Creating an Animation (Optional)
# ---------------------

# Uncomment the following block to create an animation of Volume Rendering over time

=begin
using Oceananigans, GLMakie, JLD2, Statistics

# Define downsampling factor
downsample_factor = 2  # Adjust based on performance

# Downsample grid coordinates
xC_ds = xC[1:downsample_factor:end]
yC_ds = yC[1:downsample_factor:end]
zC_ds = zC[1:downsample_factor:end]

# Create a figure and 3D axis
fig = Figure(resolution = (800, 600))
ax = Axis3(fig[1, 1],
           xlabel = "x",
           ylabel = "y",
           zlabel = "z",
           title = "Vorticity Magnitude Volume Rendering Over Time")

# Initialize the volume plot with empty data
initial_volume = NaN * ones(length(xC_ds), length(yC_ds), length(zC_ds))
volume_plot = volume!(ax, xC_ds, yC_ds, zC_ds, initial_volume;
                      algorithm = :absorption,
                      colormap = :inferno,
                      transparency = true)

# Add lighting for better depth perception
light!(ax.scene, position = Vec3f(1, 1, 1), color = Colorant"white", intensity = 1.0)

# Adjust the camera for a better view
ax.camera.projection = :perspective
cam3d!(ax.scene, eyeposition = Vec3f(2, 2, 2))

# Define the update function for each frame
function update_volume!(i)
    @info "Rendering frame $i of $(length(times))"
    
    # Extract the 3D data at time index i and downsample
    ω_mag_field = ω_mag_timeseries[i]
    ω_mag_data = collect(ω_mag_field[:, :, :])
    ω_mag_ds = ω_mag_data[1:downsample_factor:end, 1:downsample_factor:end, 1:downsample_factor:end]
    
    # Normalize the data
    ω_mag_normalized = ω_mag_ds ./ maximum(ω_mag_ds)
    
    # Update the volume plot data
    volume_plot[4][] = ω_mag_normalized
    
    # Update the title with the current time
    ax.title = "Vorticity Magnitude Volume at time = $(round(times[i], digits=2))"
end

# Create the animation by iterating over all time steps
record(fig, "vorticity_animation.mp4", 1:length(times); framerate = 10) do i
    update_volume!(i)
end

@info "Animation saved as vorticity_animation.mp4"
=end

# ---------------------
# End of Script
# ---------------------
