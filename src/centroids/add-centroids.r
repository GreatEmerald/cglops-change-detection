library(pbapply)

# Add 300 m centroids to reference data
AddCentroids = function(df, grid=300)
{
    CentroidStrings = pbapply(df, 1, function(row){
        system(paste0("python centroids/identify_centroids_LC100_grid.py ", row[["x"]], " ", row[["y"]], " | grep -Po '(?<=", grid, "m PROBA-V UTM pixel: ).*'"), TRUE)
        }, cl=4)
    CentroidMatrix = apply(do.call(rbind, strsplit(c(CentroidString, CentroidString), ", ")), 2, as.numeric)
    return(cbind(df, centroid_x=CentroidMatrix[,1], centroid_y=CentroidMatrix[,2]))
}
