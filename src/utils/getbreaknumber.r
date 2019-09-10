# Util function: calculate the size of breaks so that the minimum time between them amounts to a year
GetBreakNumber = function(dates)
{
    1/((as.numeric(difftime(max(dates), min(dates), units="weeks")))/52.18)
}
