# Script for determining which BFAST/change detection method works the best,
# according to change validation data.

# TODO: Split into functions, move into a subfolder.
# TODO: Make individual runs into an RMarkdown file for ease of reference.

source("bfast-cal/functions.r")

RawData = LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv")
ChangeIDs = as.data.frame(RawData)[RawData$change_at_300m == "yes","sample_id"]
length(ChangeIDs) # 75: total number of (non-unique) change points
length(unique(ChangeIDs)) # 62: unique change locations
plot(RawData[RawData$change_at_300m == "yes","sample_id"]) # All Africa

EVI_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"))
plot(GetMatrixFromSF(EVI_8d_16int)[1,]~GetDates8d(), type="l")

EVI_8d_16int_full = MergeAllYears(EVI_8d_16int, RawData)
ChangeIDs_UTM = as.data.frame(EVI_8d_16int_full)[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]
length(ChangeIDs_UTM) # 50: total number of points in UTM
length(unique(ChangeIDs_UTM)) # 39: unique change locations
length(unique(as.data.frame(EVI_8d_16int_full)[,"sample_id"])) # out of 1900
plot(EVI_8d_16int_full[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]) # UTM zones only

NDMI_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NDMI_8d_Int16")
NDMI_8d_16int_full = MergeAllYears(NDMI_8d_16int, RawData)
plot(GetMatrixFromSF(NDMI_8d_16int)[1,]~GetDates8d(), type="l")

NDMI_16d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NDMI_16d_Int16")
NDMI_16d_16int_full = MergeAllYears(NDMI_16d_16int, RawData)

NIRv_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_8d_Int16")
NIRv_8d_16int_full = MergeAllYears(NIRv_8d_16int, RawData)
plot(GetMatrixFromSF(NIRv_8d_16int)[1,]~GetDates8d(), type="l")

NIRv_16d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_16d_Int16")
NIRv_16d_16int_full = MergeAllYears(NIRv_16d_16int, RawData)

NIRv_16d_byte = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_16d_Byte")
NIRv_16d_byte_full = MergeAllYears(NIRv_16d_byte, RawData)

MODDetectBreaks(GetMatrixFromSF(EVI_8d_16int)[EVI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NDMI_8d_16int)[NDMI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NIRv_16d_16int)[NIRv_16d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NIRv_16d_16int)[NIRv_16d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="BIC", plot=TRUE)

MyBreak = MODDetectBreaks(GetMatrixFromSF(NIRv_8d_16int)[NIRv_8d_16int$sample_id == unique(ChangeIDs_UTM)[2],], breaks="BIC", plot=TRUE)
MyReference = as.data.frame(RawData)[RawData$sample_id==unique(ChangeIDs_UTM)[2] & RawData$change_at_300m=="yes","year_fraction"]
abline(v=MyReference, col="blue")
IsBreakInTargetYear(MyBreak, MyReference)

EVI_8d_16int_defaults = TestMODBreakpointDetection(EVI_8d_16int_full[EVI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
EVI_8d_16int_LWZ = TestMODBreakpointDetection(EVI_8d_16int_full[EVI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)
NDMI_8d_16int_defaults = TestMODBreakpointDetection(NDMI_8d_16int_full[NDMI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
NDMI_8d_16int_LWZ = TestMODBreakpointDetection(NDMI_8d_16int_full[NDMI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)
NIRv_8d_16int_defaults = TestMODBreakpointDetection(NIRv_8d_16int_full[NIRv_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
NIRv_8d_16int_LWZ = TestMODBreakpointDetection(NIRv_8d_16int_full[NIRv_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)

FPStats(EVI_8d_16int_defaults)
FPStats(EVI_8d_16int_LWZ) # Better positive predictive value
FPStats(NDMI_8d_16int_defaults)
FPStats(NDMI_8d_16int_LWZ)
FPStats(NIRv_8d_16int_defaults)
FPStats(NIRv_8d_16int_LWZ)

# This runs with LWZ
VIResult = TestVIs(c(
    "EVI_16d_Int16", "EVI_8d_Int16", "NDMI_16d_Int16", "NDMI_8d_Int16",
    "NIRv_16d_Byte", "NIRv_16d_Int16", "NIRv_8d_Byte", "NIRv_8d_Int16")
    )
FPStats(VIResult[VIResult$vi=="EVI_16d_Int16",]$bfast_guess, VIResult[TestResult$vi=="EVI_16d_Int16",]$changed)

# Get a table of all statistics
tapply(1:nrow(VIResult), VIResult$vi, function(x)FPStats(VIResult[x,]$bfast_guess, VIResult[x,]$changed))

# Test with BIC
BICResult = TestVIs(c("EVI_16d_Int16", "NDMI_16d_Int16", "NIRv_16d_Byte", "NIRv_16d_Int16"), breaks="BIC")
tapply(1:nrow(BICResult), BICResult$vi, function(x)FPStats(BICResult[x,]$bfast_guess, BICResult[x,]$changed))

EVI_unburnt = EVI_8d_16int_full[-GetBurntIndices(EVI_8d_16int_full),]
ChangeIDs_UTM_unburnt = as.data.frame(EVI_unburnt)[EVI_unburnt$change_at_300m == "yes","sample_id"]
length(unique(ChangeIDs_UTM_unburnt)) # 32
IDs_burnt = GetBurntIDs(EVI_8d_16int_full)

BICResult = TestVIs(c("EVI_16d_Int16", "NDMI_16d_Int16", "NIRv_16d_Byte", "NIRv_16d_Int16"), breaks="BIC")
tapply(1:nrow(BICResult), BICResult$vi, function(x)FPStats(BICResult[x,]$bfast_guess, BICResult[x,]$changed))

SCTestTest = TestParams(list(
    list(breaks="LWZ", scsig=0.25), # Sensitivity 0.44, specificity 0.90, a bit better
    list(breaks="LWZ", scsig=0.75), # Sensitivity 0.44, specificity 0.90
    list(breaks="BIC", scsig=0.01), # Sensitivity 0.56, specificity 0.39, a bit better
    list(breaks="BIC", scsig=0.001) # All non-break
    ), "NIRv_16d_Byte")
tapply(1:nrow(SCTestTest), SCTestTest$call, function(x)FPStats(SCTestTest[x,]$bfast_guess, SCTestTest[x,]$changed))

SCTestTest = TestParams(list(
    list(breaks="LWZ", scsig=0.15), # Sensitivity 0.41, specificity 0.90, worse than p<0.25
    list(breaks="BIC", scsig=0.005), # All non-break
    list(breaks="BIC", scsig=0.007), # All non-break
    list(breaks="BIC", scsig=0.003) # All non-break
), "NIRv_16d_Byte")
tapply(1:nrow(SCTestTest), SCTestTest$call, function(x)FPStats(SCTestTest[x,]$bfast_guess, SCTestTest[x,]$changed))

LWZSubset = SCTestTest[SCTestTest$call == "breaks LWZ, scsig 0.15",]
FNSamples = unique(LWZSubset[LWZSubset$changed & !LWZSubset$bfast_guess,]$sample_id)
for (i in 1:length(FNSamples))
{
    MODDetectBreaks(c(GetMatrixFromSF(NIRv_16d_byte[NIRv_16d_byte$sample_id == FNSamples[i],])), scrange=NULL, plot=TRUE)
    abline(v=LWZSubset[LWZSubset$changed & LWZSubset$sample_id == FNSamples[i],]$year_fraction, col="blue")
    title(FNSamples[i])
} # 1363551 is a mismatch of break definition
# A lot of cases like 1363169 will get detected properly in the future

# Try scrange

SCRangeTest = TestParams(list(
    list(breaks="LWZ", scrange=c(2015, 2019)), # Sensitivity 0.29, specificity 0.94, misses a lot
    list(breaks="BIC", scrange=c(2015, 2019))  # Sensitivity 0.44, specificity 0.60, it's a bad version of LWZ
), "NIRv_16d_Byte")

Result = SCRangeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCRangeTest = TestParams(list(
    list(breaks="LWZ", scrange=c(2016, 2018)), # Sensitivity 0.29, specificity 0.95, misses a lot
    list(breaks="BIC", scrange=c(2016, 2018))  # Sensitivity 0.32, specificity 0.66, even worse
), "NIRv_16d_Byte")

Result = SCRangeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

# Try different test types

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="Rec-CUSUM") # Sensitivity 0.38, specificity 0.94, it's even more specific
    ,list(breaks="BIC", sctype="Rec-CUSUM") # Sensitivity 0.47, specificity 0.62, better than limiting the range
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="RE") # Sensitivity 0.38, specificity 0.94, no change
    ,list(breaks="BIC", sctype="RE") # Sensitivity 0.47, specificity 0.62, no change
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="Score-CUSUM") # Sensitivity 0.38, specificity 0.94, no change
    ,list(breaks="BIC", sctype="Score-CUSUM") # Sensitivity 0.47, specificity 0.62, no change
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypePTest = TestParams(list(
    list(breaks="BIC", sctype="Rec-CUSUM", scsig=0.2) # Sensitivity 0.53, specificity 0.46, balanced but still not great, similar to OLS-MOSUM with 0.01
    ,list(breaks="BIC", sctype="Rec-CUSUM", scsig=0.1) # Sensitivity 0.47, specificity 0.55, worse
), "NIRv_16d_Byte")

Result = SCTypePTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))
