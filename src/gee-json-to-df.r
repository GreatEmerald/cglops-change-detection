# Convert GEE JSON output into a data.frame and combine

library(rjson)
library(foreach)
library(iterators)

nullToNA <- function(x) {
    x[sapply(x, is.null)] <- NA
    return(x)
}

# Converts a JSON list in a text file (with the first row as header) to a data.frame
JSONFileToDF = function(filename)
{
    MyList = fromJSON(file=filename, simplify=FALSE)
    MyNA <- lapply(MyList[-1], nullToNA)
    
    Header = unlist(MyList[[1]])
    MyDF = foreach (ListRow = iter(MyNA), .combine=rbind, .multicombine=TRUE) %do%
    {
        Result = data.frame(ListRow, stringsAsFactors=FALSE)
        names(Result) = Header
        Result
    }
    return(MyDF)
}

LoadGEEJSON = function(pattern=glob2rx("breakpoint_validation_ts_landsat_*.json"), qa_keep=c(66,130,322,386))
{
    JSONFiles = list.files("../data/", pattern=pattern, full.names=TRUE)
    DFs = lapply(JSONFiles, JSONFileToDF)
    
    # A weird but working way to make merge() accept a list
    CompleteDF = foreach(DF=iter(DFs), .combine=function(x, y){merge(x, y, all=TRUE)}) %do% DF
    
    # Remove all missing cases
    CompleteDF = CompleteDF[complete.cases(CompleteDF),]
    
    # Convert milliseconds into POSIXct and sort by it
    CompleteDF$time = as.POSIXct(CompleteDF$time/1000, origin="1970-01-01 00:00.00 UTC")
    CompleteDF = CompleteDF[order(CompleteDF$time),]
    
    # Filter clouds
    # Remove all except those marked as 66,130,322,386
    CompleteDF = CompleteDF[CompleteDF$pixel_qa %in% qa_keep,]
    
    # Calculate indices
    CompleteDF$ndvi = (CompleteDF$nir - CompleteDF$red) / (CompleteDF$nir + CompleteDF$red)
    CompleteDF$ndmi = (CompleteDF$nir - CompleteDF$swir) / (CompleteDF$nir + CompleteDF$swir)
    CompleteDF$evi = 2.5 * (CompleteDF$nir*0.0001 - CompleteDF$red*0.0001) / ((CompleteDF$nir*0.0001 + 6*CompleteDF$red*0.0001 - 7.5*CompleteDF$blue*0.0001)+1)
    
    return(CompleteDF)
}
