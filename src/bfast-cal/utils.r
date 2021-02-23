# Utilities, small functions to make handing time series easier

# Calculate dates from number of elements in the input
# Returns a ts object
GetTS = function(data, frequency=23, start=2009, years=10)
{
    stopifnot(is.vector(data)) # Univariate only
    # 8-daily frequency is 46, 16-daily is 23, we have 10 years of data
    if (is.null(frequency))
        frequency = length(data)/years
    return(ts(data, start=start, frequency = frequency))
}

# Returns a Date object
GetDates = function(...)
{
    TS = GetTS(...)
    return(as.Date(date_decimal(as.numeric(time(TS)))))
}

GetDates8d = function(...)
{
    return(GetDates(1:460))
}

GetDates16d = function(...)
{
    return(GetDates(1:230))
}

# Utility to extract only the time series table from an input full sf object
GetMatrixFromSF = function(sf)
{
    # Assumes that the pattern of columns is XNNNN.NN.NN
    TSNames = grep(glob2rx("X????.??.??"), names(sf), value = TRUE)
    as.matrix(as.data.frame(sf)[, TSNames])
}

# Same for lists (from apply on rows)
GetMatrixFromList = function(mylist)
{
    # Assumes that the pattern of columns is XNNNN.NN.NN
    TSIdx = grep(glob2rx("X????.??.??"), names(mylist))
    t(as.matrix(unlist(mylist[TSIdx])))
}

# Generic
GetMatrix = function(obj)
{
    if ("sf" %in% class(obj))
        return(GetMatrixFromSF(obj))
    if (is.list(obj))
        return(GetMatrixFromList(obj))
}
