# 4) Did BFAST predict a break in the given year? TRUE/FALSE.
# BreakTimes: all breaks detected by BFAST.
# TargetYear: time at which we want to test.
# Threshold: how much fuzziness is allowed to still consider a break detected.
# This function needs to be run for every year and every point.
# If there was no break for that year, and BFAST predicted one, should return FALSE.
IsBreakInTargetYear = function(BreakTimes, TargetYear, threshold=1)
{
    #if (is.na(TargetYear))
    #    TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
    # TODO: Needs to be updated; previously, lack of break meant that the break time would be set to NA,
    # but now it's no longer the case, it's change_at_300m = FALSE and reference_year=<year>
    return(any(BreakTimes > TargetYear - threshold & BreakTimes < TargetYear+threshold))
}

# 4b) Vectorised version (takes a list of MODDetectBreaks() output and a column of target years)
# Returns a column of whether BFSAT predicted the break at that time or not.
VectorisedIsBreakInTargetYear = function(BreakList, threshold=0.5, TY=TargetYears)
{
    i = 1:length(BreakList)
    return(sapply(i, function(i){return(IsBreakInTargetYear(BreakList[[i]], TY[i], threshold=threshold))}))
}


# 5) Get statistics.
FPStats = function(predictions, truth = NULL, round=3)
{
    # If we get a data.frame, try to automagically determine what is predicted and what is true
    if (is.data.frame(predictions) && is.null(truth))
    {
        if ("change_at_300m" %in% names(predictions)) {
            truth = predictions$change_at_300m == "yes"
        } else if ("changed" %in% names(predictions)) {
            truth = predictions$changed
        } else stop("Could not determine the column with truth, pass truth explicitly")
        predictions = predictions$bfast_guess
    }
    
    # We predicted a break and it was a break
    TruePositiveCount = sum(predictions & truth, na.rm=TRUE)
    # We predicted a break but there wasn't one (overprediction)
    FalsePositiveCount = sum(predictions & !truth, na.rm=TRUE)
    # We predicted no break, and there were none
    TrueNegativeCount = sum(!predictions & !truth, na.rm=TRUE)
    # We predicted no break, but there was one (we missed it)
    FalseNegativeCount = sum(!predictions & truth, na.rm=TRUE)
    # Percent of true positives out of all change
    Sensitivity = TruePositiveCount / (TruePositiveCount + FalseNegativeCount) # AKA Recall, Previously TruePositiveRate
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    # Percent of false positive out of no change
    FalsePositiveRate = FalsePositiveCount / sum(!truth, na.rm=TRUE) # False positive rate or alpha or p-value or Type I Error
    PositiveProportion = TruePositiveCount / FalsePositiveCount
    PositiveLikelihood = Sensitivity / FalsePositiveRate # Likelihood Ratio for Positive Tests
    Precision = TruePositiveCount / (TruePositiveCount + FalsePositiveCount) # AKA positive predictive value
    Accuracy = (TruePositiveCount + TrueNegativeCount) / length(truth)
    F1Score = 2 * (Precision * Sensitivity)/ (Precision + Sensitivity)
    return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
        Sensitivity=round(Sensitivity, round), Specificity=round(Specificity, round),
        Precision=round(Precision, round), F1Score=round(F1Score, round),
        FalsePositiveRate=round(FalsePositiveRate, round),
        PositiveProportion=round(PositiveProportion, round),
        PositiveLikelihood=round(PositiveLikelihood, round), Accuracy=round(Accuracy, round)))
}

#' Utility to run FPStats() on combined dataframes from TestParams()
#' 
#' @param df data.frame that contains information about the call, truth and predictions
#' @return data.frame of FPStats, with one row per call
FPStatsPerParam = function(df)
{
    do.call(rbind, by(df, list(df$call), FPStats))
}
