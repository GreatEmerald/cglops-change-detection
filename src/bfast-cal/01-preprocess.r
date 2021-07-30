library(sf)
library(raster)
library(lubridate)
# Functions for preprocessing data prior to running break detection

# Import the validation CSVs; they contain one entry per year that gives the approximate time of the break in the year.
# Different CSVs are formatted differently and need to be handled differently, so one function per dataset.

# 1) Load the CSV. Ideally we only need an st_read() to get an sf object
LoadReferenceDataAfrica = function(path="../data/ValidationPoints-AfricaPriority.csv")
{
    # Now our centroids are in another file, so we need to do more.
    DataBadCoords = LoadReferenceData(path, check=FALSE)
    CentroidData = read.csv("../data/2019_11_06_training_data_100m.csv")
    Data = merge(DataBadCoords, CentroidData, by.x="sample_id", by.y="sampleid")
    Data = Data[!is.na(Data$change_at_300m),] # Remove NAs
    return(Data)
}

#' Load generic reference data
#' 
#' Assumes that the data is good quality.
#' Otherwise, write your own preprocessing function that calls this one.
#' 
#' @param input If character, path to a CSV file. Otherwise, can be a ready data.frame.
#' @param xname Name of the X column
#' @param yname Name of the Y column
#' @param check Whether to check the result for consistency. NA means check and emit warnings, TRUE means enforce
LoadReferenceData = function(input, xname="centroid_x", yname="centroid_y", check=NA)
{
    if (is.character(input)) {
        Data = st_read(input, stringsAsFactors = FALSE,
            options=c(paste0("X_POSSIBLE_NAMES=", xname), paste0("Y_POSSIBLE_NAMES=", yname)))
    } else Data = input
    NumCols = c("rowid", "location_id", "sample_id", "bare", "burnt", "crops",
                "fallow_shifting_cultivation", "grassland", "lichen_and_moss", "shrub",
                "snow_and_ice", "tree", "urban_built_up", "water", "wetland_herbaceous",
                "not_sure", "reference_year", "x", "y", "centroid_x", "centroid_y", "year_fraction")
    NumCols = c(NumCols, grep(glob2rx("X????.??.??"), names(Data), value = TRUE))
    for (ColName in NumCols)
        if (ColName %in% names(Data))
            suppressWarnings(Data[[ColName]] <- as.numeric(Data[[ColName]]))
    st_crs(Data) = 4326
    
    # Reorder, some functions rely on subsequent steps to indicate change
    Data = Data[order(Data$sample_id, Data$reference_year),]
    # Remove duplicates
    UnessentialCols = c("rowid", "field_1")
    Data = Data[!duplicated(Data[,!names(Data) %in% UnessentialCols]),]
    
    if (is.na(check)) try(CheckReferenceData(Data)) else if (check) CheckReferenceData(Data)
    return(Data)
}

#' Check reference data.frame for inconsistencies
#' 
#' @param data The data.frame to check. This is the full reference data, not just unique locations
CheckReferenceData = function(data)
{
    Consistent = TRUE
    
    # So far there should be 4 years of data, so we should have 4 rows per each entry
    if(!all(table(data$sample_id) == 4))
    {
        print(which(table(data$sample_id) != 4))
        cat("Error: Above samples do not have 4 entries\n")
        Consistent = FALSE
    }
    
    if (!all(data$change_at_300m == "yes" | data$change_at_300m == "no"))
    {
        print(which(data$change_at_300m != "yes" & data$change_at_300m != "no"))
        cat("Error: Change at 300 m is neither yes nor no\n")
        Consistent = FALSE
    }
    
    if (any(is.na(data$year_fraction)))
    {
        print(which(is.na(data$year_fraction)))
        cat("Error: missing year fraction\n")
        Consistent = FALSE
    }
    if (!all(na.omit(data$year_fraction) > 2014.5 & na.omit(data$year_fraction) < 2019.5))
    {
        print(which(na.omit(data$year_fraction) < 2014.5 | na.omit(data$year_fraction) > 2019.5))
        cat("Error: Year fractions out of range\n")
        Consistent = FALSE
    }
    
    if (!all(data[data$reference_year==2015,]$change_at_300m=="no"))
    {
        print(which(data[data$reference_year==2015,]$change_at_300m!="no"))
        cat("Error: Samples labelled as change in 2015\n")
        Consistent = FALSE
    }
    
    YearSteps = data[1:nrow(data)-1,]$reference_year - data[2:nrow(data),]$reference_year
    ConsecutiveSamples = data[1:nrow(data)-1,]$sample_id == data[2:nrow(data),]$sample_id
    #ConsecutiveValidations = data[1:nrow(data)-1,]$validation_id == data[2:nrow(data),]$validation_id
    if (!all(YearSteps == -1 | !ConsecutiveSamples))
    {
        print(which(YearSteps != -1 & ConsecutiveSamples))
        cat("Error: Found years that are out of order! The data is either not sorted properly or not time-filled.\n")
        Consistent = FALSE
    }
    
    if (Consistent) cat("Data file is consistent\n") else stop("Data file is inconsistent, see errors above!")
}

# 2) Extract time series data from the coordinates,
# and cache it in a CSV/GPKG so that we don't need to do that again.
# This is where we select different VIs.
# TODO: Check the case when a point is in more than one UTM zone
# The output is an sf data.frame, with rows being unique points,
# and columns being timesteps, first columns being x, y, and sample_id, last being geometry.
# Data is the value of the vegetation index.
LoadVITS = function(pointlocs, vi="EVI_8d_Int16", sourcedir="/data/users/Public/greatemerald/modis-utm/input-vrt/", prefix="", force=FALSE, ...)
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
        InputVRTs = list.files(file.path(sourcedir, vi), pattern=glob2rx("*.vrt"), full.names = TRUE)
        
        pbi = 1
        pb = txtProgressBar(pbi, length(InputVRTs), 1, style=3)
        OutDF = NULL
        for (UTMfile in InputVRTs)
        {
            print(paste("Processing", UTMfile))
            VIMosaic = brick(UTMfile)
            VIMosaic = setZ(VIMosaic, GetDates(1:nlayers(VIMosaic), ...))
            names(VIMosaic) = GetDates(1:nlayers(VIMosaic), ...)
            # Reproject all points to this UTM zone
            PointsUTM = st_transform(pointlocs, crs(VIMosaic))
            # Which ones are inside the UTM zone?
            PointsInZone = st_contains(st_as_sfc(st_bbox(VIMosaic)), PointsUTM)[[1]]
            if (length(PointsInZone) <= 0)
                next
            # Extract those
            ChangedVITS = extract(VIMosaic, PointsUTM[PointsInZone,])
            # Put back into sf with orginal coords
            OutDF = rbind(OutDF, cbind(pointlocs[PointsInZone,c("x", "y", "sample_id")], ChangedVITS))
            pbi=pbi+1
            setTxtProgressBar(pb, pbi)
        }
        close(pb)
        
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
    Result = merge(df, st_set_geometry(data, NULL), "sample_id")
    return(Result[order(Result$sample_id, Result$reference_year),])
}
