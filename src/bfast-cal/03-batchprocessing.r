library(foreach)
library(doParallel)
source("../src/bfast-cal/04-validation.r")

# 4c) Wrapper for running 3+4 in one.
# Result is a logical vector saying how many times we accurately predicted a break,
# and how many times not. One value per year.
# VITS is the full VI TS, all years should be there.
TestMODBreakpointDetection = function(VITS, threshold=1, freq=23, quiet=FALSE, ...)
{
    # Output into a new column
    VITS$bfast_guess = NA
    
    if (!quiet) {
        pbi = 0
        pb = txtProgressBar(pbi, length(unique(VITS$sample_id)), style = 3)
    }
    # Detect the breaks in a loop over unique points (don't want to run bfast multiple times)
    for (i in unique(VITS$sample_id))
    {
        SampleMatrix = GetMatrixFromSF(VITS[VITS$sample_id == i,])
        BreakTimes = MODDetectBreaks(GetTS(SampleMatrix[1,], frequency = freq), ..., quiet=TRUE)
        
        if (!quiet) {
        pbi = pbi + 1
        setTxtProgressBar(pb, pbi)
        }
        
        if (length(BreakTimes) < 2 && is.na(BreakTimes))
            next # Already set to NA
        for (year in as.data.frame(VITS)[VITS$sample_id == i,"year_fraction"])
        {
            VITS[VITS$sample_id == i & VITS$year_fraction == year, "bfast_guess"] = IsBreakInTargetYear(BreakTimes, year, threshold)
        }
    }
    if (!quiet)
        close(pb)
    
    return(VITS)
}


#' Function to test statistics on various BFAST0N parameters
#' 
#' @param ParamLists List of parameter lists. Each list will be run on a single core.
#' @param comment Some reference string, e.g. VI name.
#' @param data A data.frame that contains time series information and change information for each year.
#' 
#' @return A data.frame with BFAST guesses vs truth that can be input into FPStats().
TestParams = function(ParamLists, comment="", data)
{
    Result = foreach(ParamList = ParamLists, .combine=rbind, .multicombine = TRUE, .verbose=TRUE) %dopar%
    {
        callstring = paste(names(ParamList), ParamList, collapse=", ")
        ParamList$plot = FALSE
        ParamList$VITS = data
        rm(data)
        gc()
        Result = do.call("TestMODBreakpointDetection", ParamList)
        Result$changed = Result$change_at_300m == "yes"
        Result$comment = as.factor(comment) # Just for reference
        Result$call = callstring
        Result
        # KeepCols = c("sample_id", "year_fraction", "changed", "bfast_guess", "call", "comment")
        # OptionalCols = c("changeclass", "coarsechange")
        # for (OC in OptionalCols)
        #     if (any(names(Result)==OC))
        #         KeepCols = c(KeepCols, OC)
        # Result[,KeepCols]
    }
    return(Result)
}
