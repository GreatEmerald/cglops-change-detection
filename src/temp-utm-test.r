# test UTM stuff

library(raster)
library(sf)
UTMTile = brick("/data/mep_cg1/MOD_S10/PBV_S5/additional_VIs/UTM26/26PQB_ProbaV_UTM_100m_2014-01-01_to_2018-03-31_NDMI.tif")
# -2^15 is fill value
ExampleLayer = UTMTile[[1]]
ExampleLayer[ExampleLayer==-2^15] = NA
plot(ExampleLayer)
hist(ExampleLayer)
EndDate = as.Date("2018.03.31", format="%Y.%m.%d")
click()
MyPoint = st_point(c(780614.1, 1651343))
MyPoint = st_set_geometry(data.frame(a="a"), st_sfc(MyPoint))
MyTS = extract(UTMTile, MyPoint)
MyTS[MyTS == -2^15] = NA
plot(c(MyTS), type="l")
MyTSO = ts(c(MyTS), end = c(2018, 90/5), frequency = 365/5)
StartDate = start(MyTSO)
StartDate[2] = StartDate[2]*5
StartDate = as.Date(paste(StartDate, collapse=" "), format="%Y %j")
plot(MyTSO)

MyDates = seq.Date(from=StartDate+1, by=5, along.with = MyTSO) # Off by one
plot(c(MyTS)~MyDates, type="l")


