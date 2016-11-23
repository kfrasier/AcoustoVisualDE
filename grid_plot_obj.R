# given the argument fill (the covariate vector to use as the fill) and a name,
# return a geom_polygon object
# fill must be in the same order as the polygon data
library(plyr)
grid_plot_obj <- function(fill, name, sp){
  
  # what was the data supplied?
  names(fill) <- NULL
  row.names(fill) <- NULL
  data <- data.frame(fill)
  names(data) <- name
  
  spdf <- SpatialPolygonsDataFrame(sp, data)
  spdf@data$id <- rownames(spdf@data)
  spdf.points <- fortify(spdf, region="id")
  spdf.df <- join(spdf.points, spdf@data, by="id")
  
  # seems to store the x/y even when projected as labelled as "long" and "lat"
  spdf.df$x <- spdf.df$long
  spdf.df$y <- spdf.df$lat
  
  geom_polygon(aes_string(x="x",y="y",fill=name, group="group"), data=spdf.df)
}