breakval = st_read("../data/valgroup27_exports_20190118_LandCoverChangeDetection_after20181201_withC... (1).csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"), stringsAsFactors=FALSE)
breakval = breakval[-nrow(breakval),] # Remove the NA at the end
st_crs(breakval) = 4326

breakval$subpixel_center_x = as.numeric(breakval$subpixel_center_x) # set the right data types
breakval$subpixel_center_y = as.numeric(breakval$subpixel_center_y)

BVUniq = breakval[!duplicated(breakval$validation_id),]
plot(BVUniq["name.1"])

LCChange = breakval[breakval$name=="CGLOPS Land Cover Change Detection validation - Mark LC change",]
LCCUniq = LCChange[!duplicated(LCChange$validation_id),]
table(LCCUniq$name.1) #130 change points total
table(LCCUniq$confidence) # 11 unsure

# Fix LCChange location IDs by taking the mean of all subpixels
LCChangeAdjusted = by(LCChange, list(LCChange$location_id, LCChange$name.1),
          function(ChangeFrac){
              FirstRow = ChangeFrac[1,]
              FirstRow$sample_x = mean(ChangeFrac$sample_x)
              FirstRow$sample_y = mean(ChangeFrac$sample_y)
              return(FirstRow)
              })
LCChangeAdjusted = LCChangeAdjusted[-which(sapply(LCChangeAdjusted, is.null))] # Remove NULLs
LCChangeAdjusted = Reduce(rbind, LCChangeAdjusted) # rbind all
plot(LCChangeAdjusted)
table(LCChangeAdjusted$name.1) # 132 locations of change (130 + 2 with change from 2 years)

BVUniq[BVUniq$validation_id %in% LCCUniq$validation_id, "changed"] = TRUE
BVUniq$changed[is.na(BVUniq$changed)] = FALSE
plot(BVUniq["changed"])
CN = breakval[breakval$name=="CGLOPS Change Validation (NEW)",]
table(CN$name.1)

######

# Mapping comments to true/false change
Comments = levels(BVUniq$comment)
write.csv(Comments, "comments.csv")

CommentMap = read.csv("../data/comments (1).csv", stringsAsFactors = FALSE)
CommentMap = CommentMap[-1,] # remove empty row
CommentMap[CommentMap == "LC Change"] = "LC change" # Fix capitals
CommentMap[CommentMap == "no sure"] = "not sure"
CommentMap[CommentMap == "maybe dynamic"] = "possible LC dynamics"
CommentMap[CommentMap == "maybe LC dynamic"] = "possible LC dynamics"
CommentMap[CommentMap == "maybe LC dynamic"] = "possible LC dynamics"
CommentMap[CommentMap == "change not sure"] = "possible LC change"
CommentMap[CommentMap == "potential water dynamics"] = "possible water dynamics"
CommentMap[CommentMap == "change in crop cultivation cycle"] = "LC dynamics"
unique(c(CommentMap$X.1, CommentMap$X.2, CommentMap$X.3))

CommentMap$WaterDynamics = CommentMap$X.3 == "possible water dynamics"
CommentMap$ConfidentWaterDynamics = CommentMap$X.1 == "water dynamics" | CommentMap$X.3 == "water dynamics"
WaterDynamicsComments = CommentMap[CommentMap$ConfidentWaterDynamics, "x"]
WaterDynamicsUnsureComments = CommentMap[CommentMap$WaterDynamics, "x"]

WaterDynIDs = BVUniq[BVUniq$comment %in% WaterDynamicsComments, "location_id"]
WaterDynUnsureIDs = BVUniq[BVUniq$comment %in% WaterDynamicsUnsureComments, "location_id"]
cat(as.integer(as.character(WaterDynIDs$location_id)))
cat(as.integer(as.character(WaterDynUnsureIDs$location_id)))

#####

# We have locations of change, now to have locations of "no change"
# This should come from the locations with comments
# 5 categories: LC change, no LC change, LC dynamics (i.e. misc), water dynamics, fire
# Add columns ChangeType (factor) + ChangeYear (integer)
UnmarkedPoints = breakval[!breakval$location_id %in% LCChangeAdjusted$location_id,]
UnmarkedUniq = UnmarkedPoints[!duplicated(UnmarkedPoints$location_id),] # 475
UnmarkedReasons = UnmarkedPoints[UnmarkedPoints$name == "CGLOPS Change Validation (NEW)",] # 277 points with reasons
table(UnmarkedReasons$name.1) # All some kind of change

# Points without reasons but with comments
RemainingUnmarked = UnmarkedPoints[!UnmarkedPoints$location_id %in% UnmarkedReasons$location_id,]
# 198 points that say "no change"
RemUnmSub = RemainingUnmarked[RemainingUnmarked$name == "CGLOPS Land Cover Change Detection validation - Did the LC type at 100m x 100m change?",]
table(RemUnmSub$name.1)
nrow(RemainingUnmarked[!duplicated(RemainingUnmarked$validation_id),]) # out of 198, so it's all of them


# Apply the comment mask to know what is real change and what is not
ChangeVal = rbind(LCChangeAdjusted, UnmarkedReasons, RemUnmSub) # 607 due to two points being different
ChangeVal$ChangeType = "no LC change" # Set all to no change
table(ChangeVal$ChangeType)
ChangeVal[1:nrow(LCChangeAdjusted),"ChangeType"] = "LC change" # Marked changes are change
table(ChangeVal$ChangeType)
# Apply the built-in change info
ChangeVal[ChangeVal$name.1 == "Fires", "ChangeType"] = "Fire dynamics"
table(ChangeVal$ChangeType)
ChangeVal[ChangeVal$name.1 %in% c("Bare->Cropland", "Forest loss/shrub loss-> Grassland",
    "Natural Vegetation-> Bare", "Wetland/grassland -> Cropland",
    "Crop abandonment/fallow", "Not sure LC transition"), "ChangeType"] = "LC change"
table(ChangeVal$ChangeType)
ChangeVal[ChangeVal$name.1 %in% c("Water expansion", "Water reduction", "Wetland degradation"), "ChangeType"] = "Water dynamics"
table(ChangeVal$ChangeType)

# Apply comment info
# Apply uncertainty
CommentMap$confidence = ""
CommentMap[CommentMap$X.1 %in% c("possible LC change", "possible LC dynamics", "possible water dynamics"), "confidence"] = "-1"
CommentMap[CommentMap$X.2 %in% c("possible LC change", "possible LC dynamics", "possible water dynamics"), "confidence"] = "-1"
CommentMap[CommentMap$X.3 %in% c("possible LC change", "possible LC dynamics", "possible water dynamics"), "confidence"] = "-1"
CommentMap$confidence

# These are otherwise no change, so we can just mark all LC dynamics as no change later
CommentMap[CommentMap$X.2 %in% c("possible LC dynamics", "LC dynamics"), "X.1"] = "LC dynamics"
CommentMap[CommentMap$X.3 %in% c("possible LC dynamics", "LC dynamics"), "X.1"] = "LC dynamics"
CommentMap[CommentMap == "water dynamics"] = "Water dynamics"
CommentMap[CommentMap == "fire"] = "Fire dynamics"
CommentMap[CommentMap$X.1 == "possible LC change", "X.1"] = "LC change"
unique(CommentMap$X.1)

for (i in 1:nrow(CommentMap))
{
    Comment = CommentMap[i, "x"]
    ChangeVal[ChangeVal$comment == Comment, "ChangeType"] = CommentMap[i, "X.1"]
    ChangeVal[ChangeVal$comment == Comment & ChangeVal$confidence != "-1", "confidence"] = CommentMap[i, "confidence"]
}
table(ChangeVal$ChangeType)
ChangeVal$ChangeType = as.factor(ChangeVal$ChangeType)

# Dates
ChangeVal$ChangeYear = NA
ChangeVal[grep("2017", ChangeVal$comment), "ChangeYear"] = 2017
ChangeVal[grep("2016", ChangeVal$comment), "ChangeYear"] = 2016
ChangeVal[ChangeVal$name.1 == "LC change 2016", "ChangeYear"] = 2016
ChangeVal[ChangeVal$name.1 == "LC change 2017", "ChangeYear"] = 2017

NoChange = ChangeVal[ChangeVal$ChangeType == "no LC change",] # Mostly interannual variability

barplot(table(ChangeVal$ChangeType), main="Types of change in collected samples")
ConfidentChange = ChangeVal[ChangeVal$confidence == "",]
barplot(table(ConfidentChange$ChangeType), main="Types of change in collected samples (confident)")

# Set sensible data types and export
str(ChangeVal)
ChangeVal$userid = as.integer(ChangeVal$userid)
ChangeVal$validation_id = as.integer(ChangeVal$validation_id)
ChangeVal$expert_rating = NULL
ChangeVal$viewed_ge = NULL
ChangeVal$submission_item_id = as.integer(ChangeVal$submission_item_id)
ChangeVal$enhancement = as.integer(ChangeVal$enhancement)
ChangeVal$skip_reason = as.integer(ChangeVal$skip_reason)
ChangeVal$sample_id = as.integer(ChangeVal$sample_id)
ChangeVal$location_id = as.integer(ChangeVal$location_id)

ExportVal = ChangeVal
class(ExportVal) = "data.frame"
ExportVal$geometry = NULL
write.csv(ExportVal, "../data/ValidationPoints-December.csv", row.names=FALSE)
