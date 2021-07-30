# Preprocessing for additional global data

source("centroids/add-centroids.r")

# Should link to the other preprocessing functions later
PreprocessAdditionalTraining = function(path="../data/training_data_100m_20200309_V4_no_time_gaps_additional.csv")
{
    AdditionalPointDF = read.csv(path)
    # Lots of duplicates for some reason, remove them
    AdditionalPointDF = AdditionalPointDF[!duplicated(AdditionalPointDF),]
    AdditionalPointDF = AddCentroids(AdditionalPointDF)
    write.csv(AdditionalPointDF, "../data/training_data_100m_20200309_V4_no_time_gaps_centroids_additional.csv")
    
    AdditionalPoints = st_as_sf(AdditionalPointDF, coords = c("centroid_x", "centroid_y"), dim="XY")
}
