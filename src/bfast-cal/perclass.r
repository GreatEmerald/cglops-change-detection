library(pbapply)

# Create a change class column
AddChangeClassCol = function(data)
{
    data$changeclass = NA
    
    fromclasses = data[which(data$change_at_300m=="yes")-1, ][["dominant_lc"]]
    toclasses = data[which(data$change_at_300m=="yes"), ][["dominant_lc"]]
    
    data[which(data$change_at_300m=="yes"),][["changeclass"]] = paste(fromclasses, "to", toclasses)
    return(data)
}

# Get a data frame only with the particular class of change and all the no-changes
FilterChange = function(data, changeclass, column="changeclass") data[data[[column]]==changeclass|is.na(data[[column]]),]

# Run FilterChange+FPStats to get a stat table, needs a changeclass column
# Use column="changeprocess" for the coarse change classification
FPStatsPerClass = function(data, cl=4, column="changeclass")
{
    ChangeClasses = names(sort(table(data[[column]]), decreasing = TRUE))
    PerClassStats = pblapply(ChangeClasses, function(cls) FPStats(FilterChange(data, cls, column=column)), cl=cl)
    PerClassStats = do.call(rbind, PerClassStats)
    rownames(PerClassStats) = ChangeClasses
    return(PerClassStats)
}

#' Add a column that clusters change classes into three:
#' between vegetation, between non-vegetation, and between vegetation and non-vegetation
#' 
#' Depends on an existing changeclass column from AddChangeClassCol
#' 
#' @param data data.frame to be modified
#' @return data.frame with an extra column "changeprocess"
AddChangeProcessCol = function(data)
{
    VegetationClasses = c("tree", "grassland", "shrub", "crops", "wetland_herbaceous")
    
    ChangesMatrix = do.call(rbind, strsplit(data$changeclass, " to "))
    # Get a boolean matrix of whether it is vegetation change
    VegetationChange = matrix(ChangesMatrix %in% VegetationClasses, ncol=2)
    # Define change processes
    ChangeProcesses = rep("between non-vegetation", nrow(data))
    ChangeProcesses[VegetationChange[,1] | VegetationChange[,2]] = "between non-vegetation and vegetation"
    ChangeProcesses[VegetationChange[,1] & VegetationChange[,2]] = "between vegetation"
    ChangeProcesses[is.na(data$changeclass)] = NA
    return(cbind(data, changeprocess = as.factor(ChangeProcesses)))
}
