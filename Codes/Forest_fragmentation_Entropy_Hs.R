
# Load necessary libraries
library(terra)
library(ggplot2)

# Load the forest raster
raster_file <- "/gpfs1/data/satellite/sentinel/germany/forest-classification/v004/forest_mask_from_Blickensdoerfer_2022_Dominant_Species_Class_R20m.tif"
forest_raster <- rast(raster_file)
print("Loaded forest raster")

# Define tile size and overlap
tile_size <- 1000  # Size of the central area
overlap_size <- 16  # Overlap for the 31x31 pixel window

# Function to create tiles with overlap
create_tiles_with_overlap <- function(forest_raster, tile_size, overlap_size) {
  r_ext <- ext(forest_raster)
  r_res <- res(forest_raster)

  x_steps <- seq(r_ext[1], r_ext[2], by = tile_size * r_res[1])
  y_steps <- seq(r_ext[3], r_ext[4], by = tile_size * r_res[2])

  tiles <- list()
  for (i in seq_along(x_steps)[-length(x_steps)]) {
    for (j in seq_along(y_steps)[-length(y_steps)]) {
      tile_extent <- ext(
        x_steps[i] - overlap_size * r_res[1], x_steps[i + 1] + overlap_size * r_res[1],
        y_steps[j] - overlap_size * r_res[2], y_steps[j + 1] + overlap_size * r_res[2]
      )
      tiles[[length(tiles) + 1]] <- crop(forest_raster, tile_extent)
    }
  }
  print("Tiles created")
  return(tiles)
}

# Function to calculate spatial entropy based on provided formula
spatial_entropy_Hs <- function(values, pixel_size = 20, scaling_factor = 5, na.rm = TRUE) {
  if (length(values) != 31 * 31) {
    stop("Error: The window is not 31x31 pixels. Check the focal function configuration.")
  }
  
  values[is.na(values)] <- 0
  
  #Check if the central pixel is forest
  
    central_pixel <- values[ceiling(length(values) / 2)]
    if (central_pixel == 0) return(NA)  # If the central pixel is not forest, return NA
  
  forest_pixels <- sum(values == 1)
  non_forest_pixels <- sum(values == 0)
  if (forest_pixels == 0 || non_forest_pixels == 0) return(0)
  
  total_pixels <- forest_pixels + non_forest_pixels
  p_forest <- forest_pixels / total_pixels
  p_non_forest <- non_forest_pixels / total_pixels
  
  matrix_size <- sqrt(length(values))
  values_matrix <- matrix(values, nrow = matrix_size, byrow = TRUE)
  transitions <- 0
  for (i in 1:(matrix_size - 1)) {
    for (j in 1:(matrix_size - 1)) {
      if (values_matrix[i, j] != values_matrix[i + 1, j]) transitions <- transitions + 1
      if (values_matrix[i, j] != values_matrix[i, j + 1]) transitions <- transitions + 1
    }
  }
  
  max_transitions <- 2 * (matrix_size - 1) * matrix_size
  edge_length <- (transitions / max_transitions) * pixel_size * scaling_factor
  
  forest_indices <- which(values == 1, arr.ind = TRUE)
  non_forest_indices <- which(values == 0, arr.ind = TRUE)
  if (is.null(dim(forest_indices))) forest_indices <- matrix(forest_indices, nrow = 1)
  if (is.null(dim(non_forest_indices))) non_forest_indices <- matrix(non_forest_indices, nrow = 1)
  
  centroid_forest <- colMeans(forest_indices) * pixel_size
  centroid_non_forest <- colMeans(non_forest_indices) * pixel_size
  
  euclidean_dist <- sqrt(sum((centroid_forest - centroid_non_forest)^2))
  distance <- euclidean_dist
  
  Hs_forest <- (edge_length / distance) * p_forest * log2(p_forest)
  Hs_non_forest <- (edge_length / distance) * p_non_forest * log2(p_non_forest)
  fragmentation <- -(Hs_forest + Hs_non_forest)
  
  H_proportion <- -(p_forest * log2(p_forest) + p_non_forest * log2(p_non_forest))
  
  weight_base <- 0.5
  weight_fragmentation <- 1.0
  H_shannon_weighted <- weight_base * H_proportion + weight_fragmentation * fragmentation
  
  homogeneity_penalty <- ifelse(p_non_forest > 0.8, 0.2, 0)
  H_shannon_weighted <- H_shannon_weighted + homogeneity_penalty
  
  Hs <- max(0, min(H_shannon_weighted, 1))
  
  return(Hs)
}

# Process and save each tile
process_and_save_tile <- function(tile_with_overlap, tile_without_overlap, tile_index, window_size = 31) {
  total_forest_pixel_count <- global(tile_with_overlap, fun = function(values) sum(values == 1, na.rm = TRUE))[[1]]

  if (total_forest_pixel_count == 0) {
    print(paste("Skipping tile", tile_index, "- no forest pixels found"))
    return(NULL)
  }

  # Apply the focal function with a 31x31 window
  entropy <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = spatial_entropy_Hs)

  # Mask for areas without forest pixels
  pixel_count <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = function(x) sum(x == 1))
  entropy_masked <- mask(entropy, pixel_count, maskvalue = 0)

  entropy_sum <- global(entropy_masked, sum, na.rm = TRUE)[[1]]
  if (is.na(entropy_sum) || entropy_sum == 0) {
    print(paste("Skipping entropy for tile", tile_index, "- no valid entropy calculated"))
    gc()
    return(NULL)
  }

  entropy_masked_final <- crop(entropy_masked, tile_without_overlap)

  # Save the processed tile
  writeRaster(entropy_masked_final, paste0("/gpfs1/work/acostave/Entropy2/forest_entropy_masked_tile_", tile_index, ".tif"), overwrite = TRUE)

  print(paste("Processed and saved tile", tile_index))

  # Free memory
  rm(entropy, pixel_count, entropy_masked)
  gc()
}

# Main loop to create, process, and save tiles
tiles_with_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 16)
tiles_without_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 0)

tile_paths <- vector("list", length(tiles_with_overlap))
for (i in seq_along(tiles_with_overlap)) {
  tile_paths[[i]] <- process_and_save_tile(tiles_with_overlap[[i]], tiles_without_overlap[[i]], i, window_size = 31)
  gc()  # Garbage collection after each tile
}

print("Tile processing complete")

# Filter out NULL paths in case some tiles were skipped
tile_paths <- Filter(Negate(is.null), tile_paths)

# Load the processed tiles and create a mosaic
tile_rasters <- lapply(tile_paths, rast)
forest_entropy_mosaic <- do.call(mosaic, tile_rasters)

# Save the mosaic as a new raster file
output_mosaic_path <- "/gpfs1/work/acostave/Entropy2/forest_entropy_mosaic.tif"
writeRaster(forest_entropy_mosaic, output_mosaic_path, overwrite = TRUE)

print("Mosaic created and saved.")


#######################TILES WITH MOSAIC AT THE END##############################
# # Load necessary libraries
# library(terra)
# library(ggplot2)
# 
# # Load the forest raster
# raster_file <- "/gpfs1/data/satellite/sentinel/germany/forest-classification/v004/forest_mask_from_Blickensdoerfer_2022_Dominant_Species_Class_R20m.tif"
# forest_raster <- rast(raster_file)
# print("Loaded forest raster")
# 
# # Define tile size and overlap
# tile_size <- 1000  # Size of the central area
# overlap_size <- 16  # Overlap for the 31x31 pixel window
# 
# # Function to create tiles with overlap
# create_tiles_with_overlap <- function(forest_raster, tile_size, overlap_size) {
#   r_ext <- ext(forest_raster)
#   r_res <- res(forest_raster)
#   
#   x_steps <- seq(r_ext[1], r_ext[2], by = tile_size * r_res[1])
#   y_steps <- seq(r_ext[3], r_ext[4], by = tile_size * r_res[2])
#   
#   tiles <- list()
#   for (i in seq_along(x_steps)[-length(x_steps)]) {
#     for (j in seq_along(y_steps)[-length(y_steps)]) {
#       tile_extent <- ext(
#         x_steps[i] - overlap_size * r_res[1], x_steps[i + 1] + overlap_size * r_res[1],
#         y_steps[j] - overlap_size * r_res[2], y_steps[j + 1] + overlap_size * r_res[2]
#       )
#       tiles[[length(tiles) + 1]] <- crop(forest_raster, tile_extent)
#     }
#   }
#   print("Tiles created")
#   return(tiles)
# }
# 
# # Function to calculate spatial entropy based on provided formula
# spatial_entropy_Hs <- function(values, pixel_size = 20, edge_factor = 2, max_distance = 876, na.rm = TRUE) {
#   if (length(values) != 31 * 31) {
#     stop("Error: The window is not 31x31 pixels. Check the focal function configuration.")
#   }
#   
#   values[is.na(values)] <- 0  # Replace NAs with 0
#   
#   # Check if the central pixel is forest
#   central_pixel <- values[ceiling(length(values) / 2)]
#   if (central_pixel == 0) return(NA)  # If the central pixel is not forest, return NA
#   
#   # Count forested and non-forested pixels
#   forest_pixels <- sum(values == 1)
#   non_forest_pixels <- sum(values == 0)
#   
#   if (forest_pixels == 0 || non_forest_pixels == 0) return(0)  # Return 0 if there's no mix
#   
#   total_pixels <- forest_pixels + non_forest_pixels
#   p_forest <- forest_pixels / total_pixels
#   p_non_forest <- non_forest_pixels / total_pixels
#   
#   # Configure the matrix of values
#   matrix_size <- sqrt(length(values))
#   values_matrix <- matrix(values, nrow = matrix_size, byrow = TRUE)
#   
#   # Calculate the edge length (number of forest/non-forest transitions)
#   transitions <- 0
#   for (i in 1:(matrix_size - 1)) {
#     for (j in 1:(matrix_size - 1)) {
#       if (values_matrix[i, j] != values_matrix[i + 1, j]) transitions <- transitions + 1
#       if (values_matrix[i, j] != values_matrix[i, j + 1]) transitions <- transitions + 1
#     }
#   }
#   edge_length <- transitions * pixel_size * edge_factor
#   
#   # Calculate centroids of forested and non-forested pixels
#   forest_indices <- which(values == 1, arr.ind = TRUE)
#   non_forest_indices <- which(values == 0, arr.ind = TRUE)
#   
#   # Return 0 if we don't have both forested and non-forested pixels
#   if (length(forest_indices) == 0 || length(non_forest_indices) == 0) return(0)
#   
#   # Convert to matrices if single coordinates are returned (edge cases)
#   if (is.null(dim(forest_indices))) forest_indices <- matrix(forest_indices, nrow = 1)
#   if (is.null(dim(non_forest_indices))) non_forest_indices <- matrix(non_forest_indices, nrow = 1)
#   
#   # Coordinates of the centroids
#   centroid_forest <- colMeans(forest_indices) * pixel_size
#   centroid_non_forest <- colMeans(non_forest_indices) * pixel_size
#   
#   # Calculate the distance between centroids and apply max limit
#   euclidean_dist <- sqrt(sum((centroid_forest - centroid_non_forest)^2))
#   distance <- min(euclidean_dist, max_distance)
#   
#   # Calculate the weighted entropy according to the provided formula
#   Hs_forest <- (edge_length / distance) * p_forest * log2(p_forest)
#   Hs_non_forest <- (edge_length / distance) * p_non_forest * log2(p_non_forest)
#   
#   # Combine forest and non-forest entropy and apply the negative sign as per the formula
#   H_shannon_weighted <- -(Hs_forest + Hs_non_forest)
#   
#   # Ensure that entropy values are normalized between 0 and 1
#   Hs <- max(0, min(H_shannon_weighted, 1))
#   
#   return(Hs)
# }
# 
# # Process and save each tile
# process_and_save_tile <- function(tile_with_overlap, tile_without_overlap, tile_index, window_size = 31) {
#   total_forest_pixel_count <- global(tile_with_overlap, fun = function(values) sum(values == 1, na.rm = TRUE))[[1]]
#   
#   if (total_forest_pixel_count == 0) {
#     print(paste("Skipping tile", tile_index, "- no forest pixels found"))
#     return(NULL)
#   }
#   
#   # Apply the focal function with a 31x31 window
#   entropy <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = spatial_entropy_Hs)
#   
#   # Mask for areas without forest pixels
#   pixel_count <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = function(x) sum(x == 1))
#   entropy_masked <- mask(entropy, pixel_count, maskvalue = 0)
#   
#   entropy_sum <- global(entropy_masked, sum, na.rm = TRUE)[[1]]
#   if (is.na(entropy_sum) || entropy_sum == 0) {
#     print(paste("Skipping entropy for tile", tile_index, "- no valid entropy calculated"))
#     gc()
#     return(NULL)
#   }
#   
#   entropy_masked_final <- crop(entropy_masked, tile_without_overlap)
#   
#   # Save the processed tile
#   output_tile_path <- paste0("/gpfs1/work/acostave/Entropy/forest_entropy_masked_tile_", tile_index, ".tif")
#   writeRaster(entropy_masked_final, output_tile_path, overwrite = TRUE)
#   
#   print(paste("Processed and saved tile", tile_index))
#   
#   # Free memory
#   rm(entropy, pixel_count, entropy_masked)
#   gc()
#   
#   return(output_tile_path)  # Return the path to the saved tile
# }
# 
# # Main loop to create, process, and save tiles
# tiles_with_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 16)
# tiles_without_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 0)
# 
# # List to store the paths of all saved tiles
# tile_paths <- vector("list", length(tiles_with_overlap))
# 
# for (i in seq_along(tiles_with_overlap)) {
#   tile_path <- process_and_save_tile(tiles_with_overlap[[i]], tiles_without_overlap[[i]], i, window_size = 31)
#   if (!is.null(tile_path)) {
#     tile_paths[[i]] <- tile_path
#   }
#   gc()  # Garbage collection after each tile
# }
# 
# print("Tile processing complete")
# 
# # Create the mosaic from the saved tiles
# tile_rasters <- lapply(tile_paths, rast)  # Load each tile as a raster
# mosaic_raster <- do.call(mosaic, tile_rasters)  # Combine tiles into a single mosaic
# 
# # Save the mosaic
# output_mosaic_path <- "/gpfs1/work/acostave/Entropy/forest_entropy_mosaic.tif"
# writeRaster(mosaic_raster, output_mosaic_path, overwrite = TRUE)
# 
# print(paste("Mosaic saved to", output_mosaic_path))




# ####THIS IS THE FINALLY COODEEEEEEE#########
# 
# # Load necessary libraries
# library(terra)
# library(ggplot2)
# 
# # Load the forest raster
# raster_file <- "/gpfs1/data/satellite/sentinel/germany/forest-classification/v004/forest_mask_from_Blickensdoerfer_2022_Dominant_Species_Class_R20m.tif"
# forest_raster <- rast(raster_file)
# print("Loaded forest raster")
# 
# # Define tile size and overlap
# tile_size <- 1000  # Size of the central area
# overlap_size <- 16  # Overlap for the 31x31 pixel window
# 
# # Create tiles with overlap
# create_tiles_with_overlap <- function(forest_raster, tile_size, overlap_size) {
#   r_ext <- ext(forest_raster)
#   r_res <- res(forest_raster)
# 
#   x_steps <- seq(r_ext[1], r_ext[2], by = tile_size * r_res[1])
#   y_steps <- seq(r_ext[3], r_ext[4], by = tile_size * r_res[2])
# 
#   tiles <- list()
#   for (i in seq_along(x_steps)[-length(x_steps)]) {
#     for (j in seq_along(y_steps)[-length(y_steps)]) {
#       tile_extent <- ext(x_steps[i] - overlap_size * r_res[1], x_steps[i + 1] + overlap_size * r_res[1],
#                          y_steps[j] - overlap_size * r_res[2], y_steps[j + 1] + overlap_size * r_res[2])
#       tiles[[length(tiles) + 1]] <- crop(forest_raster, tile_extent)
#     }
#   }
#   print("Tiles created")
#   return(tiles)
# }
# 
# # Create tiles with and without overlap
# tiles_with_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 16)
# tiles_without_overlap <- create_tiles_with_overlap(forest_raster, tile_size = 1000, overlap_size = 0)
# 
# # Entropy function (# edge factor 2 -> amplifies the importance of edge lenght)
# spatial_entropy_Hs <- function(values, pixel_size = 20, edge_factor = 2, distance_factor = 1, max_distance = 876, na.rm = TRUE) {
#   if (length(values) != 31 * 31) {
#     stop("Error: The window is not 31x31 pixels. Check the focal function configuration.")
#   }
# 
#   values[is.na(values)] <- 0  # Replace NAs with 0
# 
#   # Check if the central pixel is forest
#   central_pixel <- values[ceiling(length(values) / 2)]
#   if (central_pixel == 0) return(NA)  # If the central pixel is not forest, return NA
# 
#   # Count forested and non-forested pixels
#   forest_pixels <- sum(values == 1)
#   non_forest_pixels <- sum(values == 0)
# 
#   if (forest_pixels == 0 || non_forest_pixels == 0) return(0)  # Return 0 if there's no mix
# 
#   total_pixels <- forest_pixels + non_forest_pixels
#   p_forest <- forest_pixels / total_pixels
#   p_non_forest <- non_forest_pixels / total_pixels
# 
#   # Configure the matrix of values
#   matrix_size <- sqrt(length(values))
#   values_matrix <- matrix(values, nrow = matrix_size, byrow = TRUE)
# 
#   # Calculate the edge length (number of forest/non-forest transitions)(fragementation or complexity)
#   transitions <- 0
#   for (i in 1:(matrix_size - 1)) {
#     for (j in 1:(matrix_size - 1)) {
#       if (values_matrix[i, j] != values_matrix[i + 1, j]) transitions <- transitions + 1
#       if (values_matrix[i, j] != values_matrix[i, j + 1]) transitions <- transitions + 1
#     }
#   }
#   edge_length <- transitions * pixel_size * edge_factor
# 
#   # Calculate centroids of forested and non-forested pixels
#   forest_indices <- which(values == 1, arr.ind = TRUE)
#   non_forest_indices <- which(values == 0, arr.ind = TRUE)
# 
#   if (length(forest_indices) == 0 || length(non_forest_indices) == 0) return(0)
# 
#   if (is.null(dim(forest_indices))) forest_indices <- matrix(forest_indices, nrow = 1)
#   if (is.null(dim(non_forest_indices))) non_forest_indices <- matrix(non_forest_indices, nrow = 1)
# 
#   # Coordinates of the centroids
#   centroid_forest <- colMeans(forest_indices) * pixel_size
#   centroid_non_forest <- colMeans(non_forest_indices) * pixel_size
# 
#   # Calculate the distance between centroids
#   euclidean_dist <- sqrt(sum((centroid_forest - centroid_non_forest)^2))
# 
#   # Limit the maximum distance
#   distance <- min(euclidean_dist, max_distance)
# 
#   # Shannon entropy for forest and non-forest
#   Hs_forest <- p_forest * log2(p_forest)
#   Hs_non_forest <- p_non_forest * log2(p_non_forest)
#   H_shannon <- -(Hs_forest + Hs_non_forest)
# 
#   # Calculate spatial entropy based on `edge_length`, `H_shannon`, and controlled distance
#   Hs_raw <- H_shannon * edge_length / (distance * distance_factor)
# 
#   # Normalize using a `max_Hs` based on the maximum edge and distance limit
#   max_edge_length <- 37200  # Maximum possible edge in meters for 31x31 window
#   max_Hs <- max_edge_length / max_distance
# 
#   # Final normalization to ensure values between 0 and 1
#   Hs_normalized <- Hs_raw / max_Hs
#   Hs <- max(0, min(Hs_normalized, 1))
# 
#   return(Hs)
# }
# 
# # Process each tile and calculate entropy
# process_and_save_tile <- function(tile_with_overlap, tile_without_overlap, tile_index, window_size = 31) {
#   total_forest_pixel_count <- global(tile_with_overlap, fun = function(values) sum(values == 1, na.rm = TRUE))[[1]]
# 
#   if (total_forest_pixel_count == 0) {
#     print(paste("Skipping tile", tile_index, "- no forest pixels found"))
#     return(NULL)
#   }
# 
#   # Apply the focal function with a 31x31 window
#   entropy <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = spatial_entropy_Hs)
# 
#   # Mask for areas without forest pixels
#   pixel_count <- focal(tile_with_overlap, w = matrix(1, nrow = window_size, ncol = window_size), fun = function(x) sum(x == 1))
#   entropy_masked <- mask(entropy, pixel_count, maskvalue = 0)
# 
#   entropy_sum <- global(entropy_masked, sum, na.rm = TRUE)[[1]]
#   if (is.na(entropy_sum) || entropy_sum == 0) {
#     print(paste("Skipping entropy for tile", tile_index, "- no valid entropy calculated"))
#     gc()
#     return(NULL)
#   }
# 
#   entropy_masked_final <- crop(entropy_masked, tile_without_overlap)
# 
#   # Save the processed tile
#   writeRaster(entropy_masked_final, paste0("/gpfs1/work/acostave/Entropy/forest_entropy_masked_tile_", tile_index, ".tif"), overwrite = TRUE)
# 
#   print(paste("Processed and saved tile", tile_index))
# 
#   # Free memory
#   rm(entropy, pixel_count, entropy_masked)
#   gc()
# }
# 
# # Loop to process and calculate entropy for each tile
# for (i in seq_along(tiles_with_overlap)) {
#   process_and_save_tile(tiles_with_overlap[[i]], tiles_without_overlap[[i]], i, window_size = 31)
#   gc()  # Garbage collection after each tile
# }
# 
# print("Tile processing complete")
