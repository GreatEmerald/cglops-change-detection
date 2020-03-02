# Create a change class column
AddChangeClassCol = function(data)
{
    data$changeclass = NA
    for (i in which(data$change_at_300m=="yes"))
        data[i, ][["changeclass"]] = paste(data[i-1,][["dominant_lc"]], "to", data[i,][["dominant_lc"]])
    return(data)
}

# Get a data frame only with the particular class of change and all the no-changes
FilterChange = function(data, changeclass) data[data$changeclass==changeclass|is.na(data$changeclass),]

# Run FilterChange+FPStats to get a stat table, needs a changeclass column
FPStatsPerClass = function(data, cl=4)
{
    ChangeClasses = names(sort(table(data$changeclass), decreasing = TRUE))
    PerClassStats = pblapply(ChangeClasses, function(cls) FPStats(FilterChange(data, cls)), cl=cl)
    PerClassStats = do.call(rbind, PerClassStats)
    rownames(PerClassStats) = ChangeClasses
    return(PerClassStats)
}