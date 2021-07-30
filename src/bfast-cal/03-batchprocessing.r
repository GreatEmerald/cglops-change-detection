library(foreach)
library(doParallel)
library(pbapply)
library(data.table)
source("../src/bfast-cal/04-validation.r")

#' 3c) Wrapper for running 3+4 in one.
#' @return data.frame augmented with a column "bfast_guess" that shows what BFAST guessed.
#' @param VITS Full VI TS data.frame, all years should be there.
#' @param threshold Let for considering whether the break was detected correctly or not.
#' @param NewAccuracy Whether or not to use the new accuracy calculation method (4c)
TestMODBreakpointDetection = function(VITS, threshold=1, freq=23, quiet=FALSE,
    DetectFunction=MODDetectBreaks, NewAccuracy=TRUE, verbose=FALSE, cl=NULL, ...)
{
    # Parse function name
    if (!is.function(DetectFunction) && is.character(DetectFunction))
        DetectFunction = eval(parse(text=DetectFunction))
    # Output into a new column
    if (!NewAccuracy) {
        VITS$bfast_guess = NA
    } else {
        VITS$TP=NA
        VITS$TN=NA
        VITS$FP=NA
        VITS$FN=NA
    }
    
    if (!quiet) {
        pbi = 0
        pb = txtProgressBar(pbi, length(unique(VITS$sample_id)), style = 3)
    }
    # Detect the breaks in a loop over unique points (don't want to run bfast multiple times)
    #for (i in unique(VITS$sample_id))
    #{
    ProcessSample = function(i)
    {
        SampleMatrix = GetMatrixFromSF(VITS[VITS$sample_id == i,])
        if (verbose) {
            print(paste("Processing sample", i))
        }
        BreakTimes = DetectFunction(GetTS(SampleMatrix[1,], frequency = freq), ..., quiet=!verbose)
        
        if (!quiet) {
            pbi = pbi + 1
            setTxtProgressBar(pb, pbi)
        }
        
        SampleChunk = VITS[VITS$sample_id == i, ]
        
        if (length(BreakTimes) < 2 && is.na(BreakTimes))
            return(SampleChunk) # Already set to NA
        if (NewAccuracy) {
            TruthDates = SampleChunk[SampleChunk$change_at_300m == "yes",]$year_fraction
            if (length(BreakTimes) == 1 && !BreakTimes)
                BreakTimes = numeric(0)
            Stats = BreakConfusionStats(BreakTimes, TruthDates, threshold = threshold)
            for (Stat in c("TP", "TN", "FP", "FN"))
                SampleChunk[[Stat]] = Stats[Stat]
        } else {
            for (year in as.data.frame(VITS)[VITS$sample_id == i,"year_fraction"])
            {
                SampleChunk[SampleChunk$year_fraction == year, "bfast_guess"] = IsBreakInTargetYear(BreakTimes, year, threshold)
            }
        }
        return(SampleChunk)
    }
    #}
    ProcessedDF = pblapply(unique(VITS$sample_id), ProcessSample, cl=cl)
    #VITS = do.call(rbind, ProcessedDF)
    VITS = as.data.frame(data.table::rbindlist(ProcessedDF))
    
    if (!quiet)
        close(pb)
    
    # For the new style of accuracy assessment, we don't need repeated years
    if (NewAccuracy)
        VITS = VITS[!duplicated(VITS$sample_id),]
    
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
