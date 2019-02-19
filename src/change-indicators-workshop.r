# Stability indicator analysis
library(raster)

# UTM zones: we are interested in 29, 30, 31 north
VRTFiles = "/data/users/Public/greatemerald/probav/change-indicators/"

UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NDMI.vrt"))
AC = LoadTrueChange()

BandNames = c("std. median", "perc change forward median", "perc change backward median", "std. q10", "perc change forward q10", "perc change backward q10", "std. q90", "perc change forward q90", "perc change backward q90")

UTM29Val = extract(UTM29Mosaic, AC)
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NDMI.vrt"))
UTM30Val = extract(UTM30Mosaic, AC)
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NDMI.vrt"))
UTM31Val = extract(UTM31Mosaic, AC)

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NDMIVal=UTM29Val

# NAUC

UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NAUC.vrt"))
UTM29Val = extract(UTM29Mosaic, AC)
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NAUC.vrt"))
UTM30Val = extract(UTM30Mosaic, AC)
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NAUC.vrt"))
UTM31Val = extract(UTM31Mosaic, AC)

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NAUCVal = UTM29Val


UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NIRv.vrt"))
UTM29Val = extract(UTM29Mosaic, AC)
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NIRv.vrt"))
UTM30Val = extract(UTM30Mosaic, AC)
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NIRv.vrt"))
UTM31Val = extract(UTM31Mosaic, AC)

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NIRvVal = UTM29Val

apply(NAUCVal, 2, median, na.rm=TRUE)
apply(NDMIVal, 2, median, na.rm=TRUE)
apply(NIRvVal, 2, median, na.rm=TRUE)

apply(NAUCVal, 2, max, na.rm=TRUE)
apply(NDMIVal, 2, max, na.rm=TRUE)
apply(NIRvVal, 2, max, na.rm=TRUE)

AllChange = LoadAllReference()

UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NIRv.vrt"))
UTM29Val = extract(UTM29Mosaic, AllChange[is.na(AllChange$Change),])
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NIRv.vrt"))
UTM30Val = extract(UTM30Mosaic, AllChange[is.na(AllChange$Change),])
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NIRv.vrt"))
UTM31Val = extract(UTM31Mosaic, AllChange[is.na(AllChange$Change),])

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NIRvAllVal = UTM29Val

UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NAUC.vrt"))
UTM29Val = extract(UTM29Mosaic, AllChange[is.na(AllChange$Change),])
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NAUC.vrt"))
UTM30Val = extract(UTM30Mosaic, AllChange[is.na(AllChange$Change),])
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NAUC.vrt"))
UTM31Val = extract(UTM31Mosaic, AllChange[is.na(AllChange$Change),])

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NAUCAllVal = UTM29Val


UTM29Mosaic = brick(file.path(VRTFiles, "UTM29N-NDMI.vrt"))
UTM29Val = extract(UTM29Mosaic, AllChange[is.na(AllChange$Change),])
colnames(UTM29Val) = BandNames
UTM30Mosaic = brick(file.path(VRTFiles, "UTM30N-NDMI.vrt"))
UTM30Val = extract(UTM30Mosaic, AllChange[is.na(AllChange$Change),])
UTM31Mosaic = brick(file.path(VRTFiles, "UTM31N-NDMI.vrt"))
UTM31Val = extract(UTM31Mosaic, AllChange[is.na(AllChange$Change),])

UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM30Val[UTM29Invalid,]
UTM29Invalid = apply(UTM29Val, 1, function(x){all(is.na(x))})
UTM29Val[UTM29Invalid,] = UTM31Val[UTM29Invalid,]

NDMIAllVal = UTM29Val

apply(NAUCAllVal, 2, median, na.rm=TRUE)
apply(NDMIAllVal, 2, median, na.rm=TRUE)
apply(NIRvAllVal, 2, median, na.rm=TRUE)

apply(NAUCAllVal, 2, max, na.rm=TRUE)
apply(NDMIAllVal, 2, max, na.rm=TRUE)
apply(NIRvAllVal, 2, max, na.rm=TRUE)

apply(NAUCAllVal, 2, quantile, 0.98, na.rm=TRUE)
apply(NDMIAllVal, 2, quantile, 0.98, na.rm=TRUE)
apply(NIRvAllVal, 2, quantile, 0.98, na.rm=TRUE)

# From these, we say that we mask by NAUC: bare==0 & b0 >=115 & b3 >= 145 & b6 >= 160
# NDMI: 120, 105, 170; NIRv: 50, 25, 80
# Currently will use | instead of &
# These will show areas with potential change.
# Bare NDVI threshold is 0.17

# t-test todo: run over Senegal on PV and Landsat; if that's better than 
