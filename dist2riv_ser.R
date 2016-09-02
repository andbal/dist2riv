# Description:  Compute river influence on nearby soil
# Author(s):    Daniele La Cecilia
#                   original author
#               Johannes Brenner
#                   merged original scripts, added some functionality
#               Andrea Balotti - andrea.balotti@hotmail.it
#                   parallel implementation and code optimisation with functions
#
# Copyright (c) 2016 Andrea Balotti

# NB. 'parallel' give some boost in the tab computation (that account for ~ 87%
# of the whole computational time: tab + map).
# For the risk map computation, the serial version is still faster.

# TODO
#2) If possible parallelise raster reading (very slow for high resolution raster)
#1) MPI implementation is disabled, require more checks

dist2riv <- function(dem=dem, river=adige,  p_h = 3, read_tab = FALSE,
                     coordsys = "+proj=utm +zone=32 +ellps=WGS84", is_mpi = FALSE )
{
if (is_mpi){
    # Load library
    library("Rmpi")
    require("snow")
    # Set master and slave cores
    if (mpi.comm.rank(0) > 0) {
    sink(file="/dev/null")
    # runMPIslave()
    slaveLoop(makeMPImaster())
    mpi.quit()
    }
    cores <- mpi.universe.size()-1
} else {
    # Load library
    library("parallel")
    cores <- detectCores()-1
}
    
    # FUNC 1: raster processing
    proc_rst <- function(map,isdem=FALSE){
        # get raw value
        # s <- proc.time()
        mapVal <- getValues(map)
        # e <- proc.time()-start
        # print(e)
        writeLines("Got values from raster")
        # print(c("Required:",e))
        # eliminate NA values
        # s <- proc.time()
        nonamap   <- !is.na(mapVal)
        writeLines("map nona created.")
        mapVal_nona <- map[nonamap]
        # mapVal_nona <- getValues(map,)
        writeLines("raster nona created.")
        # e <- proc.time()-start
        # print(e)
        writeLines("Eliminated NA from raster.")
        # print(c("Required:",e))
        # retrieve coordinates no NA cell numbers
        # s <- proc.time()
        xy <- xyFromCell(object = map, cell = which(nonamap)) #object = map
        # e <- proc.time()-start
        # print(e)
        writeLines("xy coordinates written.")
        # print(c("Required:",e))
        # bind coordinates and values
        maptab <- cbind(xy, mapVal_nona)
        # ONLY for DEM: create .RData
        if (isdem) {
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
    # deltaxy <- numeric(length = length(demtab))
    deltaH <- demtab[x,3] - rivtab[ind_dist,3]
    return(c(deltaxy, deltaH))
    }
    
    # # FUNC 2a: distance (vertical) of each dem point to nearest river point
    # calc_delta_m_MOD <- function(x, demtab, rivtab) {
    #     ind_dist <- which.min(sqrt((demtab[x,1] - rivtab[,1])^2 + (demtab[x,2] - rivtab[,2])^2 ))
    #     # deltaxy <- sqrt((demtab[x,1] - rivtab[ind_dist,1])^2 + (demtab[x,2]- rivtab[ind_dist,2])^2 )
    #     deltaH <- demtab[x,3] - rivtab[ind_dist,3]
    #     # return(c(deltaxy, deltaH))
    #     return(deltaH)
    # }

    # # FUNC 3: distance as simple Euclidean Theorem
    # distf <- function(x) {
    # y <- unlist(x)
    # if (length(y)>2) {stop("Incorrect length for distf() arguments. Required row length is 2.")}
    # dist <- sqrt(y[1]^2 + y[2]^2)
    # return(dist)
    # }

    # # FUNC 4: normalisation between 0 and 1
    # # formula x_norm = (x_i - x_min) / (x_max - x_min)
    # normf <- function(x,n_max,n_min) {
    # x_norm <- (x - n_min)/(n_max - n_min)
    # return(x_norm)
    # }
    
    # MAIN
    if (read_tab) {
        # read already computed table
        load("TAB_x_y_Altit_Dist_Risk.RData")
    } else {
        # postprocess data
        #- DEM ----------------
        writeLines("\nProcessing DEM file\n")

        s0 <- proc.time()
        # cl <- makeCluster(cores)
        # clusterExport(cl,varlist=c("dem"), envir = environment())
        # clusterEvalQ(cl,expr = "proc_rst")
        demtab <- proc_rst(dem,isdem=TRUE)
        # demtab <- parApply(cl,X=dem,FUN=proc_rst,isdem=TRUE)
        # stopCluster(cl)
        e0 <- proc.time()-start
        print(e0)
        nvals <- length(demtab[,1])
        rm(dem)
        writeLines("DEM DONE & SAVED\n")
        
        #- RIVER --------------
        writeLines("Processing RIVER file\n")
        s0 <- proc.time()
        # cl <- makeCluster(cores)
        # clusterExport(cl,varlist=c("river"), envir = environment())
        # clusterEvalQ(cl,expr = "proc_rst")
        rivtab <- proc_rst(river)
        # rivtab <- parLapply(cl,X=river,fun=proc_rst,isdem=FALSE)
        # stopCluster(cl)
        e0 <- proc.time()-start
        print(e0)
        rm(river)
        writeLines("RIVER DONE & SAVED\n")
        #----------------------
    
        # run parallel or not
        # cl <- makeCluster(cores)
        # clusterExport(cl,varlist=c("demtab","rivtab"),
        #               envir = environment())
        # clusterEvalQ(cl,expr = "calc_delta_m")
        # xy_h <- parLapply(cl, X = 1:nvals, fun = calc_delta_m,
        #                   demtab = demtab, rivtab = rivtab)
        xy_h <- lapply(X = 1:nvals, FUN = calc_delta_m,
                       demtab = demtab, rivtab = rivtab)
        # stopCluster(cl)
        xy_h <- matrix(unlist(x = xy_h), ncol = 2, byrow = T)
        colnames(xy_h) <- c("delta_xy", "delta_h")
        disttab <- cbind(demtab, xy_h)
        writeLines("distance TAB DONE\n")
        
        # add 0 col
        disttab <- cbind(disttab, matrix(0, nvals, 1))
        colnames(disttab)[6] <- "RISK"
    
        # save tab in .Rdata file
        save("disttab", file="TAB_x_y_Altit_Dist_Risk.RData")
    }

    # calculate RISK (of ground water influence)
    # distance from river + apply normalisation
    writeLines("Calculate risk of ground water influence")
    
    ### Alternative Solution 2 - h as value and no normalization
    writeLines("Use heigh difference from river")
    # disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6] <- disttab[disttab[,4] < p_d & disttab[,5] < p_h ,5]
    disttab[disttab[,5] < p_h ,6] <- disttab[disttab[,5] < p_h ,5]
    
    ### Alternative Solution - h as value + normalization
    # writeLines("Use heigh difference from river")
    # risktab <- disttab[disttab[,4] < p_d & disttab[,5] < p_h ,5]
    # if (is_parallel) {
    #     start <- proc.time()
    #     cl <- makeCluster(cores)
    #     clusterExport(cl, varlist="risktab",envir = environment())
    #     clusterEvalQ(cl,expr ="normf")
    #     risktab_temp <- parSapply(cl, X = risktab, FUN = normf,
    #                               n_max=max(risktab), n_min=min(risktab))
    # } else {
    #     risktab_temp <- sapply(X = risktab, FUN = normf,
    #                            n_max=max(risktab), n_min=min(risktab))
    # }
    # disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6] <- risktab_temp
    # remove(risktab_temp)
    ###
    
    ### DON'T REMOVE - could be used in case of further improvements
    # risktab <- disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6]
    # tab_xyh <- cbind(disttab[disttab[,4] < p_d & disttab[,5] < p_h,4],
    #                  disttab[disttab[,4] < p_d & disttab[,5] < p_h,5])
    # if (is_parallel){
    #     start <- proc.time()
    #     cl <- makeCluster(cores)
    #     clusterExport(cl, varlist="tab_xyh",envir = environment())
    #     clusterEvalQ(cl,expr =c("distf","normf"))
    #     risktab_temp <- parRapply(cl, x = tab_xyh,FUN= distf)
    #     risktab <- parSapply(cl, X = risktab_temp, FUN = normf,
    #                          n_max=max(risktab_temp), n_min=min(risktab_temp))
    #     stopCluster(cl)
    #     end <- proc.time()-start
    #     print(end)
    #     remove(risktab_temp)
    # } else {
    #     start <- proc.time()
    #     risktab_temp <- apply(X = tab_xyh,MARGIN = c(1),FUN = distf)
    #     risktab <- sapply(X = risktab_temp, FUN = normf,
    #                       n_max=max(risktab_temp), n_min=min(risktab_temp))
    #     end <- proc.time()-start
    #     print(end)
    #     remove(risktab_temp)
    # }
    # disttab[disttab[,4] < p_d & disttab[,5] < p_h ,6] <- risktab
    ###
    
    # dummy raster
    load("dem.RData")
    rast <- raster(ncol = lc, nrow=lr, xmn=ex_xm, xmx=ex_xM, ymn=ex_ym, ymx=ex_yM)
    projection(rast) <- coordsys
    
    # dummy values = 0
    rast[] <- 0
    df <- data.frame(which(nonadem), disttab[,6])
    rast[which(nonadem)] <- df[,2]
    rast <- trim(x = rast, padding = 0, values = 0)
    rast[rast==0] <- NA
    
    # name of raster map
    namefile <- paste("rast_indicator_", p_h, "h.tif", sep="")
    
    # write raster map
    writeRaster(rast, namefile, format="GTiff", overwrite=TRUE)
    
    return(rast)
    
    # if (is_mpi) {
    # mpi.finalize()          # LAST COMMAND - necessary to properly execute MPI
    # }
}