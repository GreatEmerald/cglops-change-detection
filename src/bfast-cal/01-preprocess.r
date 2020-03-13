# Functions for preprocessing data prior to running break detection

# Import the validation CSVs; they contain one entry per year that gives the approximate time of the break in the year.
# Different CSVs are formatted differently and need to be handled differently, so one function per dataset.

# 1) Load the CSV. Ideally we only need an st_read() to get an sf object
LoadReferenceDataAfrica = function(path="../data/ValidationPoints-AfricaPriority.csv")
{
    # Normally if we have a good CSV, we just do this.
    #Data = st_read(path, options=c("X_POSSIBLE_NAMES=x", "Y_POSSIBLE_NAMES=y"))
    # Now our centroids are in another file, so we need to do more.
    DataBadCoords = read.csv(path)
    CentroidData = read.csv("../data/2019_11_06_training_data_100m.csv")
    Data = merge(DataBadCoords, CentroidData, by.x="sample_id", by.y="sampleid")
    Data = st_as_sf(Data, coords = c("centroid_x", "centroid_y"), dim="XY")
    Data = Data[!is.na(Data$change_at_300m),] # Remove NAs
    
    # Which columns are numeric
    NumCols = c("rowid", "location_id", "sample_id", "bare", "burnt", "crops",
                "fallow_shifting_cultivation", "grassland", "lichen_and_moss", "shrub",
                "snow_and_ice", "tree", "urban_built_up", "water", "wetland_herbaceous",
                "not_sure", "reference_year")
    for (ColName in NumCols)
    {
        Data[[ColName]] = as.numeric(Data[[ColName]])
    }
    st_crs(Data) = 4326
    return(Data)
}

# 2) Extract time series data from the coordinates,
# and cache it in a CSV/GPKG so that we don't need to do that again.
# This is where we select different VIs.
# TODO: Check the case when a point is in more than one UTM zone
# The output is an sf data.frame, with rows being unique points,
# and columns being timesteps, first columns being x, y, and sample_id, last being geometry.
# Data is the value of the vegetation index.
LoadVITS = function(pointlocs, vi="EVI_8d_Int16", sourcedir="/data/users/Public/greatemerald/modis-utm/input-vrt/", prefix="", force=FALSE)
{
    # Cache file. Has location_id, sample_id, x, y, geometry and the values
    VITSFile = paste0("../data/", prefix, vi, "-TS.gpkg")
    if (force || !file.exists(VITSFile))
    {
        print(paste("Cache file", VITSFile, "not found, generating..."))
        # Deduplicate the input. We shouldn't extract from the same point more than once
        # sample_id is unique, but also an option is to dedup on x/y
        pointlocs = pointlocs[!duplicated(pointlocs$sample_id),]
        
        # Input directory that contains all our VRTs is sourcedir,
        # inside we have VIs, and then UTM zones as VRTs
        InputVRTs = list.files(file.path(sourcedir, vi), full.names = TRUE)
        
        OutDF = NULL
        for (UTMfile in InputVRTs)
        {
            print(paste("Processing", UTMfile))
            VIMosaic = brick(UTMfile)
            VIMosaic = setZ(VIMosaic, GetDates(1:nlayers(VIMosaic)))
            names(VIMosaic) = GetDates(1:nlayers(VIMosaic))
            # Reproject all points to this UTM zone
            PointsUTM = st_transform(pointlocs, crs(VIMosaic))
            # Which ones are inside the UTM zone?
            PointsInZone = st_contains(st_as_sfc(st_bbox(VIMosaic)), PointsUTM)[[1]]
            if (length(PointsInZone) <= 0)
                next
            # Extract those
            print(system.time(ChangedVITS <- extract(VIMosaic, PointsUTM[PointsInZone,]))) # method="bilinear" takes too much RAM
            # Put back into sf with orginal coords
            OutDF = rbind(OutDF, cbind(pointlocs[PointsInZone,c("x", "y", "sample_id")], ChangedVITS))
        }
        
        # Finally, save to cache and not do that again
        st_write(OutDF, VITSFile, delete_dsn = TRUE)
    } else OutDF = st_read(VITSFile)
    
    NALocations = apply(GetMatrixFromSF(OutDF), 1, function(x) all(is.na(x)))
    OutDF = OutDF[!NALocations,] # Remove all that are only NAs
    OutDF = OutDF[!duplicated(OutDF$sample_id),] # Deduplicate
    return(OutDF)
}

# Merge the matrix of unique points with the reference data, to get all years info
MergeAllYears = function(df, data)
{
    return(merge(df, as.data.frame(data), "sample_id"))
}
