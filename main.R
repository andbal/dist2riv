# MPI option forced FALSE, require more checks.
# ### MPI PARALLELIZATION ###
# is_mpi <- FALSE          # TRUE to execute in parallel with MPI
# ###########################

library("raster")

#setwd("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis")
#source("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis/dist2riv.R")
source("GitHub/dist2riv/dist2riv.R")
# source("GitHub/dist2riv/dist2riv_ser.R")

### read dem
# dem <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis/data/dtm_2pt5m_new_clip.tif")
dem <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BrJ/01_Influence_Analysis/clip/DEM_clip1.tif")
# dem <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/LaD/03_Influence_Analysis/dtm_2pt5m.tif")
#plot(dem)

### read river tif
# adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis/data/adige_2pt5m_new_clip.tif")
# adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BaA/01_Influence_Analysis/data/adige_polygon_2pt5m.tif")
adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/BrJ/01_Influence_Analysis/clip/Adige_clip1.tif")
# adige <- raster("H:/Projekte/MONALISA/12_GIS/Workspace/LaD/03_Influence_Analysis/Adige.tif")
#plot(adige, add=T, legend=F, col="black")

### compute river influence
start <- proc.time()
rast <- dist2riv(dem = dem, river = adige, p_h = 4, read_tab = FALSE )
end <- proc.time()-start
print(end)

### Create tab time: (parallel on)
# DEM_clip1.tif             50 sec
# dtm_2pt5m_new_clip.tif    8 min
# dtm_5m.tif                2 h (4h su VSC?!)
# dtm_2pt5m.tif             unknown