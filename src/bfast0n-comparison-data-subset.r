# Select ~300 points based on maximum dissimilarity
library(FactoMineR)
library(sf)
source("../src/bfast-cal/01-preprocess.r")
source("../src/bfast-cal/perclass.r")

GlobalData = LoadReferenceData("../data/Data_Global.csv")
GlobalTS = LoadVITS(GlobalData, "", "/data/cgl_lc_bfast/Modis_VIs_GLOBAL_UTM_FINAL/", "Global")
GlobalData = MergeAllYears(GlobalData, GlobalTS)
CheckReferenceData(GlobalData)
GlobalData$continent = NULL
names(GlobalData)[names(GlobalData)=="field_30"] = "continent"
table(GlobalData$continent)
GlobalData[GlobalData$dominant_lc=="wetland_hebaceous",]$dominant_lc = "wetland_herbaceous"
GlobalData = AddChangeClassCol(GlobalData)
GlobalData = AddChangeProcessCol(GlobalData)
GlobalData = GlobalData[GlobalData$reference_year %in% 2016:2018,] # Discard any changes in 2015, reference does not track that

# Add additional statistics
# Mark all years of a change time series as "changed"
ExtraData = GlobalData
ExtraData[GlobalData$sample_id %in% GlobalData[GlobalData$change_at_300m == "yes",]$sample_id,]$change_at_300m = "yes"
# Number of observations in the time series
ExtraData$nobs = apply(GetMatrix(ExtraData), 1, function(x)sum(!is.na(x)))

# Remove duplicates
ExtraData = ExtraData[!duplicated(ExtraData$sample_id),]

# Remove too cloudy time series
ExtraData = ExtraData[ExtraData$nobs > 50,]

# Remove extra columns
ExtraData = ExtraData[,-c(2:20, 22:27, 29:30)]
ExtraData = ExtraData[,-(ncol(ExtraData)-c(2:3))]

# Get a version with only relevant columns for sampling
SamplingData = ExtraData[,c("change_at_300m", "continent", "nobs")]
SamplingData = st_set_geometry(SamplingData, NULL)

# Convert columns to factors
SamplingData$change_at_300m = factor(SamplingData$change_at_300m)
SamplingData$continent = factor(SamplingData$continent)

# Do PCA sampling, 1% of the data
(n = round(nrow(ExtraData)/100))

PCA1 = FAMD(SamplingData, graph=TRUE)$ind$coord[,1]
ResultRows = sort(PCA1)[ceiling(seq(1, length(PCA1), length.out=n))]
Sample1Perc = ExtraData[rownames(ExtraData) %in% names(ResultRows),]

#Sample1Perc = GlobalData[GlobalData$sample_id %in% sample(unique(GlobalData$sample_id),length(unique(GlobalData$sample_id))/100),]

st_write(Sample1Perc, "../data/bfast0n-benchmark-points-pca.gpkg")
