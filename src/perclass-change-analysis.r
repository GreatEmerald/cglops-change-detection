library(pbapply)
source("bfast-cal/functions.r")
source("bfast-cal/perclass.r")

# Fix issue where some years are set to 0
RawData[RawData$year_fraction < 2000, "year_fraction"] = RawData[RawData$year_fraction < 2000, "year_fraction"] + RawData[RawData$year_fraction < 2000, "reference_year"]
# Fix precision mismatch
RawData$year_fraction = round(RawData$year_fraction, 3)
MyData$year_fraction = round(MyData$year_fraction, 3)

# Merge the BFAST output with the raw data
LWZ_merged = merge(RawData, MyData[MyData$call=="breaks BIC",], by=c("sample_id", "year_fraction"))

all.equal(LWZ_merged$geometry, LWZ_merged$geom) # TRUE, drop one
LWZ_merged$geom = NULL
LWZ_merged = st_set_geometry(LWZ_merged, "geometry")

# Check output
plot(st_geometry(LWZ_merged))
FPStats(LWZ_merged)

# Run the per-class statistics
LWZ_merged = AddChangeClassCol(LWZ_merged)
sort(table(LWZ_merged$changeclass)) # Most are trees to grass, then water to grass
PerClassStats = FPStatsPerClass(LWZ_merged)
write.csv(PerClassStats, "../output/RoW-change-per-class-BIC.csv")
