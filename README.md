# dist2riv
Function for R environment that compute the river influence in the nearby soil, once maximum horizontal and vertical influence distances are given.

The area of influence is simply the square root of the sum of the square of horizontal distance and the square of vertical distance. Then, this value is normalised between 0 (max. influence) and 1 (min. influence).

Firstly is computed a table that contain the horizontal and vertical distance from the nearest river pixel for each _dtm_ cell. After that, only the cells inside the maximum distance limits are retained and normalised.

# How to use:
Load library and function:
```R
require("raster")
source('dist2riv.R')
```
Load the raster relative to terrain elevation and river:
```R
dem <- raster("path/dtm.tif")
adige <- raster("path/river.tif")
```
Compute the distance influence:
```R
img <- dist2ind(dem=dem, river=adige)
```
Visualise the influence raster:
```R
plot(img)
```
# Note
The option `is_parallel = TRUE` give some boost in the _table_ computation (that account for ~ 87% of the whole computational time: tab + map). For the map creation only, the serial version is faster so it is suggested to turn off `is_parallel` when `read_table = TRUE`.
