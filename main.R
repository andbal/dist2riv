# MPI option disabled, require more checks.
# ### MPI PARALLELIZATION ###
# is_mpi <- FALSE          # TRUE to execute in parallel with MPI
# if (is_mpi){is_parallel <- TRUE}
# ###########################

library("raster")

#setwd("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis")
source("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis/dist2riv.R")

### read dem
dem <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BrJ/01_Influence_Analysis/clip/DEM_clip1.tif")
#dem <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/LaD/03_Influence_Analysis/dtm_2pt5m.tif")
#plot(dem)

### read river tif
adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BrJ/01_Influence_Analysis/clip/Adige_clip1.tif")
#adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/LaD/03_Influence_Analysis/Adige.tif")
#plot(adige, add=T, legend=F, col="black")

### compute river influence
start <- proc.time()
rast <- dist2ind(dem = dem, river = adige, p_d = 500, p_h = 2, read_tab = FALSE,
                 is_parallel = FALSE, is_mpi = FALSE)
end <- proc.time()-start
print(end)