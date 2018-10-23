GetDatesFromDir = function(dir)
{
    dirnames = list.dirs(dir, FALSE, FALSE)
    dates = parse_date_time(grep(glob2rx("????????"), dirnames, value=TRUE), "ymd")
    return(as.Date(dates))
}

# Util function: calculate the size of breaks so that the minimum time between them amounts to a year
GetBreakNumber = function(dates)
{
    1/((as.numeric(difftime(max(dates), min(dates), units="weeks")))/52.18)
}
