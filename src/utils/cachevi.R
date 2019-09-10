library(sf)
library(raster)

source("utils/dates.r")

# Cache the time series of a vegetation index for easy access
# Returns a data.frame of VI time series for the points indicated
CacheVI = function(Locations, TSStack, VI, CacheFile=paste0("../data/", VI, "TS.csv"), Reset=FALSE, dates=GetDatesFromDir("/data/mep_cg1/MOD_S10/MOD_S10"), ...)
{
    # Cache invalidation is hard.
    # The input files may not exist any more, so we want to make sure to keep the VITSs for plotting.
    # So only update the cache if reset=TRUE.
    # Nevertheless, if the locations requested do not match the cache, also update the cache (append).
    
    if (!Reset && file.exists(CacheFile))
    {
        CurrentDataset = as.matrix(read.csv(CacheFile))#st_read(CacheFile, options=c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y"))
        # TODO: Check if we are missing any points (cache invalidation)
        return(CurrentDataset)
    }

    stopifnot(nrow(Locations) == length(dates))
    VITS = extract(TSStack, Locations)
    colnames(VITS) = as.character(dates)
    VITS = cbind(VITS, X=Locations$sample_x, Y=Locations$sample_y)
    write.csv(VITS, CacheFile, row.names=FALSE)
    
    return(VITS)
}
