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

#' 4c) A new version that takes all years as input and returns TP/FP/TN/FN.
#' It is more accurate, as a point may no longer be both a TP and FP.
#' @param PredictionDates Vector with predicted dates of breaks (year fractions)
#' @param TruthDates      Vector with reference dates of breaks (year fractions)
#' @param threshold       How much to allow the predictions to deviate while still considered correct (±years)
#' @param period          Start and end date of the period of interest (without threshold)
#' @return Vector with four values: TP, FP, TN, FN.
BreakConfusionStats = function(PredictionDates, TruthDates, threshold=0.5, period=c(2016, 2019))
{
    # Remove all predictions that fall outside of the period
    if (length(PredictionDates) > 0)
        PredictionDates = PredictionDates[PredictionDates > min(period) - threshold & PredictionDates < max(period) + threshold]
    
    PointDistance = function(point, points) any(abs(point - points) <= threshold+1e-13)
    
    # Algorithm: each category requires one rule to calculate.
    # 1) False positives:
    TruthCloseToPred = if (length(TruthDates) <= 0) {
        # If there is no true break, all predictions are far from truth
        # and all predictions are false positives
        rep(FALSE, length(PredictionDates))
    } else if (length(PredictionDates) <= 0) {
        # If there are no predictions, there are no false positives
        logical(0)
    } else {
        # Is there a true break within threshold of each predicted break?
        sapply(PredictionDates, PointDistance, TruthDates)
    }
    # Any predictions not close to a real break is a false positive
    FP = sum(!TruthCloseToPred)
    
    # 2) False negatives:
    PredCloseToTruth = if (length(PredictionDates) <= 0) {
        # If there is no prediction, all true breaks are far from predictions
        # and all true breaks are false negatives
        rep(FALSE, length(TruthDates))
    } else if (length(TruthDates) <= 0) {
        # If there are no true breaks, there are no false negatives
        logical(0)
    } else {
        # Is there a predicted break within threshold of each true break?
        sapply(TruthDates, PointDistance, PredictionDates)
    }
    # Any true break not close to a prediction is a false negative
    FN = sum(!PredCloseToTruth)
    
    # 3) True positives:
    # Any break that is close to another from both sets
    TP = min(sum(TruthCloseToPred), sum(PredCloseToTruth))
    
    # 4) True negatives:
    # The rest, given the total number of items we should have.
    # e.g. 4 - (sum of above)
    TN = floor(max(period) - min(period) + 2*threshold) - sum(FP, FN, TP)
    return(c(TP=TP, TN=TN, FP=FP, FN=FN))
}

# 5) Get statistics from the prediction result.
# This is generic enough to handle multiple formats.
FPStats = function(predictions, truth = NULL, round=3)
{
    # If we get a data.frame, try to automagically determine what is predicted and what is true
    if (is.data.frame(predictions) && is.null(truth))
    {
        if ("bfast_guess" %in% names(predictions))
        {
            predictions = predictions$bfast_guess
            if ("change_at_300m" %in% names(predictions)) {
                truth = predictions$change_at_300m == "yes"
            } else if ("changed" %in% names(predictions)) {
                truth = predictions$changed
            }
        }
    }
    # New version with stats already in the DF
    if (all(c("TP", "TN", "FP", "TP") %in% names(predictions))) {
        TruePositiveCount = sum(predictions$TP, na.rm=TRUE)
        FalsePositiveCount = sum(predictions$FP, na.rm=TRUE)
        TrueNegativeCount = sum(predictions$TN, na.rm=TRUE)
        FalseNegativeCount = sum(predictions$FN, na.rm=TRUE)
    } else {
        # We predicted a break and it was a break
        TruePositiveCount = sum(predictions & truth, na.rm=TRUE)
        # We predicted a break but there wasn't one (overprediction)
        FalsePositiveCount = sum(predictions & !truth, na.rm=TRUE)
        # We predicted no break, and there were none
        TrueNegativeCount = sum(!predictions & !truth, na.rm=TRUE)
        # We predicted no break, but there was one (we missed it)
        FalseNegativeCount = sum(!predictions & truth, na.rm=TRUE)
    }
    # Percent of true positives out of all change
    Sensitivity = TruePositiveCount / (TruePositiveCount + FalseNegativeCount) # AKA Recall, Previously TruePositiveRate
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    # Percent of false positive out of no change
    FalsePositiveRate = FalsePositiveCount / (TrueNegativeCount + FalsePositiveCount) # False positive rate or alpha or p-value or Type I Error
    PositiveProportion = TruePositiveCount / FalsePositiveCount
    PositiveLikelihood = Sensitivity / FalsePositiveRate # Likelihood Ratio for Positive Tests
    Precision = TruePositiveCount / (TruePositiveCount + FalsePositiveCount) # AKA positive predictive value
    Accuracy = (TruePositiveCount + TrueNegativeCount) /
        sum(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount)
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
