library(strucchange)
library(bfast)
data = c(0.2, 0.4, 0.8, 1, 0.8, 0.4, 0.2, 0.2, 0.4, 0.8, 1, 0.8, 0.4, 0.2, 0.2, 0.4, 0.8, 1, 0.72, 0.32, 0.14, 0.12, 0.2, 0.32, 0.3, 0.16, 0.04, 0.02, 0.02, 0.04, 0.08, 0.1, 0.08, 0.04, 0.02)
data = data + c(rnorm(ceiling(length(data)/2), sd=0.1), rnorm(floor(length(data)/2), sd=0.025))
dates = seq(2012, 2017, by=1/7)
dates = dates[1:length(data)]
t_ts = ts(data, 2012, frequency=7)
plot(t_ts)

bpp = bfastpp(t_ts, order=1)
plot(read.zoo(bpp))
bfr = breakpoints(response ~ harmon+trend, data=bpp, h=7) # Min 7 observations per segment
bfr

## This is weird: there are more breakpoints detected than there should be with a full model (harmon+trend), and only response~1 works correctly

plot(t_ts, ylab="NDVI")
for (i in 1:length(bfr$breakpoints))
    abline(v=dates[bpp$trend[[bfr$breakpoints[[i]]]]], col="red")
