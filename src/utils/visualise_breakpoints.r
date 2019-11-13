library(strucchange)
library(bfast)
source("utils/getbreaknumber.r")

MyPlotBfast = function(MyBF, PlotData=TRUE, ...)
{
    MyIter = length(MyBF$output)
    MyFit = MyBF$output[[MyIter]]$St + MyBF$output[[MyIter]]$Tt
    MyTrend = MyBF$output[[MyIter]]$Tt
    MyTimes = as.numeric(time(MyBF$output[[MyIter]]$Tt))
    MyTimes = MyTimes[!is.na(MyBF$output[[MyIter]]$Tt)]
    if (any(!is.na(MyBF$output[[MyIter]]$bp.Vt)))
    {
        MyBreaks = MyTimes[MyBF$output[[MyIter]]$bp.Vt$breakpoints]
    } else {
        MyBreaks = NULL
    }
    if (PlotData)
    {
        plot(MyBF$Yt, ...)
        lines(MyFit, col="blue")
    } else {
        plot(MyFit, col="blue", ...)
    }
    lines(MyTrend, col="green")
    abline(v=MyBreaks, col="red")
}

# Plot a given time series (pixel values) and return a breakpoint object
# 'bfast' controls whether the whole BFAST is calculated, or just the breakpoints
visualise_breakpoints = function(pixel, order=3, timestep="10-day", dates=NULL,
                                 t0=as.Date("2014-03-12"), bfast=FALSE, approxna=FALSE, breaks=NULL, ...)
{
    if (approxna)
        pixel = na.approx(pixel)
        
    if (timestep == "regular")
    {
        z = zoo(pixel,dates) # make a zoo
        bfts = as.ts(z)
        print(bfts)
        if (is.null(dates))
            dates = as.Date(as.numeric(time(bfts)))
        bpp = .bfastpp.full(bfts, order=order) # Requires a redefined .bfastpp.full
    } else {
        stopifnot(!is.null(dates))
        bfts = bfastts(pixel, dates, type = timestep)
        bpp = bfastpp(bfts, order=order)
    }
    
    #plot(read.zoo(bpp))
    if (!bfast)
    {
        bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates), breaks=breaks)
        #print(paste("Breakpoints:", toString(bf$breakpoints)))
        #print(paste("Times:", toString(bpp$time[[max(bf$breakpoints)]])))
        #print(as.integer(as.Date(bpp$time[[max(bf$breakpoints)]]) - t0)) ## Gives NA for irregular because dates do not correspond to bpp$trend in this case
        plot(response/10000~time, data=bpp, type="l", ...)#col=ifelse(getSceneinfo(names(timeseries))$sensor == "OLI", "blue", "green"), ...)#main=POI$comment[point_id], ...)
        if (is.finite(bf$breakpoints[1]))
        {
            for (i in 1:length(bf$breakpoints)){
                #print(i)
                #print(bf$breakpoints)
                Break = bpp$time[[bf$breakpoints[[i]]]]
                print(Break)
                abline(v=Break, col="red")
            }
        }
    } else {
        bfo = bfast(bfts, GetBreakNumber(dates), season="harmonic", max.iter=1)
        MyPlotBfast(bfo, ...)
        bf = bfo$breakpoints
    }
    return(bf)
}
