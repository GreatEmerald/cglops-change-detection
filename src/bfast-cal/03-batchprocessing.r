# 4c) Wrapper for running 3+4 in one.
# Result is a logical vector saying how many times we accurately predicted a break,
# and how many times not. One value per year.
# VITS is the full VI TS, all years should be there.
TestMODBreakpointDetection = function(VITS, threshold=0.5, TargetYears=TargetYears, ...)
{
    # Output into a new column
    VITS$bfast_guess = NA
    
    pbi = 0
    pb = txtProgressBar(pbi, length(unique(VITS$sample_id)), style = 3)
    # Detect the breaks in a loop over unique points (don't want to run bfast multiple times)
    for (i in unique(VITS$sample_id))
    {
        SampleMatrix = GetMatrixFromSF(VITS[VITS$sample_id == i,])
        BreakTimes = MODDetectBreaks(SampleMatrix[1,], ..., quiet=TRUE)
        
        pbi = pbi + 1
        setTxtProgressBar(pb, pbi)
        
        if (length(BreakTimes) < 2 && is.na(BreakTimes))
            next # Already set to NA
        for (year in as.data.frame(VITS)[VITS$sample_id == i,"year_fraction"])
        {
            VITS[VITS$sample_id == i & VITS$year_fraction == year, "bfast_guess"] = IsBreakInTargetYear(BreakTimes, year, threshold)
        }
    }
    close(pb)
    
    return(VITS)
}


# ParamLists is a list of lists of parameter and value pairs
TestParams = function(ParamLists, vi)
{
    Result = foreach(ParamList = ParamLists, .combine=rbind, .multicombine = TRUE, .verbose=TRUE) %dopar%
    {
        VI = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), vi)
        VI_full = MergeAllYears(VI, RawData)
        rm(VI)
        # Filter out 2015 and 2018, former can't change and latter can't be detected
        VI_full = VI_full[VI_full$reference_year %in% 2016:2017,]
        # Filter out burnt areas
        # IDs to filter out
        BurntIDs = GetBurntIDs(VI_full)
        BurntIDs = BurntIDs[!BurntIDs %in% ChangeIDs_UTM] # Do not exclude IDs that have changed
        VI_full = VI_full[!VI_full$sample_id %in% BurntIDs,] # Exclude burnt areas
        
        callstring = paste(names(ParamList), ParamList, collapse=", ")
        ParamList$plot = FALSE
        ParamList$VITS = VI_full
        rm(VI_full)
        gc()
        Result = do.call("TestMODBreakpointDetection", ParamList)
        Result$changed = Result$change_at_300m == "yes"
        Result$vi = as.factor(vi) # Just for reference
        Result$call = callstring
        Result[,c("sample_id", "year_fraction", "changed", "bfast_guess", "call", "vi")]
    }
    return(Result)
}
