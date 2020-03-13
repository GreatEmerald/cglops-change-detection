#' Plot the time series and results of bfast0n
#' 
#' @param bp     breakpoints object from strucchange::breakpoints()
#' @param bpp    data.frame output from bfastpp()
#' @param breaks Number of breaks or optimal break selection method, see strucchange::breakpoints()
#' @param bpMag  Break magnitudes, as returned by magnitude(bp)
#' @param ...    Other parameters to pass to plot()
plot.bfast0n = function(bp, bpp, breaks, bpMag=NULL, ...)
{
    # Plot the original time series
    plot(response~time, data=bpp, type="l", ...)
    # Plot the fitted model in green
    lines(fitted(bp, breaks=breaks)~bpp$time, col="green")
    
    # Get the requested breaks
    bpOptim = breakpoints(bp, breaks=breaks)
    
    if (length(bpOptim$breakpoints) > 0 && !all(is.na(bpOptim$breakpoints))) {
        bpTimes = bpp[bpOptim$breakpoints, "time"]
        bpY = bpp[bpOptim$breakpoints, "response"]
        abline(v=bpTimes, col="blue") # Detected breakpoints in blue
        # If magnitudes requested, plot whiskers denoting magnitude
        arrows(bpTimes, bpY-bpMag[,"RMSD"], bpTimes, bpY+bpMag[,"RMSD"], length=0.05, angle=90, code=3, col="blue")
    }
}
