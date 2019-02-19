library(raster)
library(bfast)
library(lubridate)

md = brick("/data/mep_cg1/MOD_S10/PBV_S5/additional_VIs/UTM28/28PFU_ProbaV_UTM_100m_2014-01-01_to_2018-03-31_NDMI.tif")

BFMBreaks = function(pixel)
{
    pixel[pixel==-2^15] = NA # Mask out values always
    endyear = year(today(tzone="Europe/Brussels"))
    Results = numeric(endyear-startyear)
    Results[] = NA
    if (all(is.na(pixel)))
        return(Results)
    
    # This is optimised for the dates that we have
    PixelTS = ts(pixel, c(2014, 6), frequency=round(365.25/5)) # 5-daily, start from date_decimal(2014.068) ~2014-01-25
    i = 1
    for (year in startyear:(endyear-1))
    {
        if (year+1 == endyear)
        {
            ShortenedTS = PixelTS
        } else {
            ShortenedTS = window(PixelTS, end=year+1)
        }
        result = tryCatch(bfastmonitor(ShortenedTS, year),
                          error = function(e){print(e); traceback(e); cat(c("Note: pixel values were: ", toString(pixel), "\n")); return(NA)})
        if (all(!is.na(result)))
            Results[i] = yday(date_decimal(result$breakpoint)) # Fractional year, so convert to day of year
        i = i + 1
    }
    return(Results)
}

mdv = getValues(md, nrow(md)/2, 1)
plot(mdv[1,], type="l")
startyear=2016
pixel = mdv[10,]
BFMBreaks(mdv[1,])
