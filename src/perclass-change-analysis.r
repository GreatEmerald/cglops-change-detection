library(pbapply)
source("bfast-cal/functions.r")
source("bfast-cal/perclass.r")

# Merge the BFAST output with the raw data (Africa)
LWZ_merged = merge(RawData, MyData, by=c("sample_id", "year_fraction"))

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
write.csv(PerClassStats, "../output/Africa-change-per-class.csv")
