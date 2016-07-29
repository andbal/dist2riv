### PROCEDURE
# require("raster")
# require("rgdal")
# source('dist2ind_BaA.R')
# dem <- raster("path/dem.tif")
# adige <- raster("path/river.tif")
# img <- dist2ind(dem=dem, river=adige)
# plot(img)

# TODO
#1) MPI implementation is disabled, require more checks

dist2ind <- function(dem=dem, river=adige,  p_d = 500, p_h = 1, read_tab = FALSE,
                     coordsys = "+proj=utm +zone=32 +ellps=WGS84",
                     is_parallel = FALSE, is_mpi = FALSE)
{
if (is_parallel){
    # if (is_mpi){
    # # Load library
    # library("Rmpi")
    # require("snow")
    # # Set master and slave cores
    # if (mpi.comm.rank(0) > 0) {
    # sink(file="/dev/null")
    # # runMPIslave()
    # slaveLoop(makeMPImaster())
    # mpi.quit()
    # }
    # cores <- mpi.universe.size()-1
    # }
library("parallel")
cores <- detectCores()-1
}
    
    # FUNC 1: raster processing
    proc_rst <- function(map,isdem=FALSE){
        # get raw value
        mapVal <- getValues(map)
        writeLines("Got values from raster\n")
        # eliminate NA values
        nonamap   <- !is.na(mapVal)
        mapVal_nona <- map[nonamap]
        writeLines("Eliminated NA from raster\n")
        # retrieve coordinates no NA cell numbers
        xy <- xyFromCell(object = map, cell = which(nonamap)) #object = map
        # bind coordinates and values
        maptab <- cbind(xy, mapVal_nona)
        # ONLY for DEM: create .RData
        if (isdem){
            nonadem <- nonamap
            # Nr. of rows and columns
            lr <- nrow(map)
            lc <- ncol(map)
            # extent of dem
            ex_xm = map@extent@xmin; ex_xM = map@extent@xmax
            ex_ym = map@extent@ymin; ex_yM = map@extent@ymax
            # save table
            save("lr", "lc", "ex_xm","ex_xM","ex_ym","ex_yM","nonadem",file="dem.RData")
        }
        return(maptab)
    }

    # FUNC 2: distance (horizontal+vertical) of each dem point to nearest river point
    calc_delta_m <- function(x, demtab, rivtab) {
    ind_dist <- which.min(sqrt((demtab[x,1] - rivtab[,1])^2 + (demtab[x,2] - rivtab[,2])^2 ))
    deltaxy <- sqrt((demtab[x,1] - rivtab[ind_dist,1])^2 + (demtab[x,2]- rivtab[ind_dist,2])^2 )
    deltaH <- demtab[x,3] - rivtab[ind_dist,3]
    return(c(deltaxy, deltaH))
    }

    # FUNC 3: distance as simple Euclidean Theorem
    distf <- function(x) {
    y <- unlist(x)
    if (length(y)>2) {stop("Incorrect length for distf() arguments. Required row length is 2.")}
    dist <- sqrt(y[1]^2 + y[2]^2)
    return(dist)
    }

    # FUNC 4: normalisation between 0 and 1
    # formula x_norm = (x_i - x_min) / (x_max - x_min)
    normf <- function(x,n_max,n_min){
    x_norm <- (x - n_min)/(n_max - n_min)
    return(x_norm)
    }
    
    # MAIN
    if (read_tab) {
        # read already computed table
        load("TAB_x_y_Altit_Dist_Risk.RData")
    } else {
        # postprocess data
        #- DEM ----------------
        writeLines("\nProcessing DEM file\n")
        demtab <- proc_rst(dem,isdem=TRUE)
        nvals <- length(demtab[,1])
        rm(dem)
        writeLines("DEM DONE & SAVED\n")
        
        #- RIVER --------------
        writeLines("Processing RIVER file\n")
        rivtab <- proc_rst(river)
        rm(river)
        writeLines("RIVER DONE & SAVED\n")
        #----------------------
    
        # run parallel or not
        if (is_parallel) {
            nvals_list <- as.list(1:nvals)
            cl <- makeCluster(cores)
            clusterExport(cl,varlist=c("calc_delta_m","demtab","rivtab","nvals_list"),
                          envir = environment())
            xy_h <- parLapply(cl, X = nvals_list, fun = calc_delta_m,
                              demtab = demtab, rivtab = rivtab)
            stopCluster(cl)
            xy_h <- matrix(unlist(x = xy_h), ncol = 2, byrow = T)
        } else {
            xy_h <- sapply(X = 1:nvals, FUN = calc_delta_m,
                           demtab = demtab, rivtab = rivtab)
            xy_h <- t(xy_h)
        }
        
        colnames(xy_h) <- c("delta_xy", "delta_h")
        disttab <- cbind(demtab, xy_h)
        writeLines("distance TAB DONE\n")
        
        # add 0 col
        disttab <- cbind(disttab, matrix(0, nvals, 1))
        colnames(disttab)[6] <- "RISK"
    
        # save tab in .Rdata file
        save("disttab", file="TAB_x_y_Altit_Dist_Risk.RData")
    }
    
    load("dem.RData")

    # calculate RISK (of ground water influence)
    # distance from river + apply normalisation
    writeLines("Calculate risk of ground water influence")
    risktab <- disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6]
    tab_xyh <- cbind(disttab[disttab[,4] < p_d & disttab[,5] < p_h,4],
                     disttab[disttab[,4] < p_d & disttab[,5] < p_h,5])
    if (is_parallel){
        start <- proc.time()
        cl <- makeCluster(cores)
        clusterExport(cl, varlist=c("tab_xyh","distf","normf"),
                      envir = environment())
        risktab_temp <- parRapply(cl, x = tab_xyh,FUN= distf)
        risktab <- parSapply(cl, X = risktab_temp, FUN = normf,
                             n_max=max(risktab_temp), n_min=min(risktab_temp))
        stopCluster(cl)
        end <- proc.time()-start
        print(end)
        remove(risktab_temp)
    } else {
        start <- proc.time()
        risktab_temp <- apply(X = tab_xyh,MARGIN = c(1),FUN = distf)
        risktab <- sapply(X = risktab_temp, FUN = normf,
                          n_max=max(risktab_temp), n_min=min(risktab_temp))
        end <- proc.time()-start
        print(end)
        remove(risktab_temp)
    }
    disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6] <- risktab
  
    # dummy raster
    rast <- raster(ncol = lc, nrow=lr, xmn=ex_xm, xmx=ex_xM, ymn=ex_ym, ymx=ex_yM)
    projection(rast) <- coordsys
    
    # dummy values = 0
    rast[] <- 0
    df <- data.frame(which(nonadem), disttab[,6])
    rast[which(nonadem)] <- df[,2]
    rast <- trim(x = rast, padding = 0, values = '0')
    rast[rast==0] <- NA
    
    # name of raster map
    namefile <- paste("rast_indicator_", p_d, "_", p_h, ".tif", sep="")
    
    # write raster map
    writeRaster(rast, namefile, format="GTiff", overwrite=TRUE)
    
    return(rast)
    
    # if (is_mpi) {
    # mpi.finalize()          # LAST COMMAND - necessary to properly execute MPI
    # }
}