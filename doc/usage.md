% Land cover change detection algorithm theoretical basis
% Dainius Masiliūnas
% 2019-08-12

# Land cover change detection

## Overview and algorithm choice

For producing yearly map updates, a change detection algorithm is needed to assess whether a given pixel is temporally stable or not.
If it is stable, then it is desirable to keep the previous land cover classification, in order to avoid spurious changes that comes from minor changes in the input data for the classifier.

There are several algorithms available for land cover change detection.
We evaluated several of those: BFAST, BFAST Monitor and t-test.
For consolidated map production, we chose [BFAST](http://dx.doi.org/10.1016/j.rse.2009.08.014), since the algorithm is able to detect all breaks in a given time series, is more precise with the estimation of breakpoint timing, and gives an uncertainty measure about the time of the breaks.
This allows for potential optimisation of the classifier, since knowing the time of the breaks allows avoiding training the classifier with metrics produced from unstable periods of change.
In addition, knowing all the breaks in the time series, rather than just the most significant one, allows for more efficient consolidation, as newly incoming data is reprocessed and thus contributes to a more precise estimation of breaks in all years simultaneously.
This allows for updating (consolidating) previously-released land cover maps as part of the map updating process for a new year.

In contrast, BFAST Monitor was designed for detecting change at the end of a time series, and thus reports only whether there was a break or not and the timing of the first break since the start of the monitoring period.
This break is not necessarily the largest break, and only a single break is reported.
Similarly, a t-test run over the history and compared with a monitoring period is a simple way to test whether there has been a break in the time series, however, it only indicates whether the mean of the monitoring period is different from the mean of the history period, with no indication about the time of the break.

## BFAST algorithm details

BFAST is used for detecting breaks in the time series of vegetation indices.
The algorithm works by first decomposing the input (cloud-free) time series into a seasonal (harmonic: sine and cosine functions, with a defined frequency order of a year or half-year), a trend and a remainder component.
Next, the time series is segmented into segments via piecewise linear regression, with a minimal segment size that is controlled by the parameter *h*.
The number of breaks in the time series is optimised by iteratively minimising the Bayesian Information Criterion (BIC).
This procedure is only run in case a structural change test indicates with *p < 0.05* that there is at least one break in the time series.
The algorithm has several configurable parameters: vegetation index that it is run on, minimal segment size *h*, components on which the breaks are being detected, and the harmonic order of the seasonal component(s).

We set the parameter *h* to equal the maximum number of timesteps within a single year of the time series.
This means that at most one break can be detected in any given year.
Also, this implies that the first break that can be detected is at least one year away from the start or end of the time series.
As such, the algorithm is not suited for detecting changes at the end of the time series, i.e. for NRT updating; in that case, BFAST Monitor is used instead.

We detected breaks on both seasonal and trend components.
The harmonic order of the seasonal component was set to 3, i.e. the mean, annual and semiannual frequencies were included in the model.
This was determined by testing the algorithm output over Sahel when using order 3 vs order 2 harmonics and comparing the change detection accuracy.

The algorithm was run on time series from MODIS 250 m Terra+Aqua reflectance product from 2009 until 2018, after calculating a number of vegetation indices from it.
The algorithm was run on three vegetation indices: NDMI, EVI and NIRv, as they result in different output that may be used complementarily.

Once the breaks are determined, the output of the algorithm is then used in further processing, in combination with the classifier output and expert rules based on extra stability metrics to determine which pixels are likely to have changed.

## Algorithm limitations

Any algorithm based on time series analysis requires a *history* period which defines how a stable time series should look like.
Deviations from such a stable state result in a reported break.
As such, a long time series history is needed for the algorithm to function properly.
Hence the algorithm is run over MODIS data, which provides a much longer history period than Proba-V.
The downside is that the spatial resolution of MODIS data is coarser (250 m instead of 100 m).
As such, areas of change that are smaller than the 250 m MMU cannot be detected; likewise, areas of no change within an otherwise changing MMU will also not be reported.
To deal with this issue, it is planned to transition from MODIS to Landsat data (30 m) as the provider of the history period.

Since BFAST requires a minimum segment size *h* (i.e. time between breaks) to be specified, no breaks can be detected for *h* time steps from the end of the time series.
In our case, setting *h* to a year means that we can only detect changes one year after they happen at the earliest.
To work around this limitation, BFAST Monitor will be used to detect breaks at the end of the time series for NRT maps, and BFAST will be used for map consolidation after a year has passed.

In addition, BFAST tends to overestimate the number of breaks in the time series, as it is sensitive to relatively small variations that are caused by seasonal variability.
To work around this limitation, we constrain the change further by combining it with classifier output (e.g. even if a break is specified, the class before and after the break may be the same, therefore the change is not considered land cover change). In addition, we combine the output with extra stability metrics and expert rules that define what class transitions are possible and which ones are not.

Furthermore, while the algorithm provides uncertainly about the timing of each break, it does not provide uncertainly about the existence of the breaks.
As an alternative, the BIC score can be used as a per-pixel metric for how well the modelled harmonics fit the actual data after it is segmented.
However, it is a metric for all the breaks in the time series, rather than each one.
<!--It is possible to further extract goodness-of-fit statistics for each segment, however, the breaks are defined as the transition points between the segments.-->

Lastly, the BFAST algorithm is computationally intensive, much more so than BFAST Monitor or t-test.
As such, the Terrascope computing cluster is used for computations.
In the future, as more validation data is gathered, the optimal vegetation index can be determined and it may no longer be necessary to run the algorithm over multiple vegetation indices.

## Implementation

The algorithm is implemented in the R programming language and is designed to run on the Terrascope platform.
To make use of the parallelisation capabilities provided by the Spark cluster set up on the Terrascope platform, the algorithm makes use of the `SparkR` package.
To make the calculations efficient when scaled across a cluster, the input data for the algorithm needs to be split into chunks.
After processing is complete, the chunks are reassembled (mosaicked) back into an output raster.

The main script that runs break detection can be found at `src/detect-breaks.R`.
The script is partially self-contained, in that it is sent to the allocated Spark driver and is then executed there.
The script for pushing R scripts to driver can be found at `src/spark-submit.sh` and is a convenience wrapper around the `spark-submit` script provided by the Spark 2 installation.
Further, the convenience wrapper `src/process-tile.sh` ensures the correct order of processing: first, input files are split into several thousand chunks locally, since it is a disk-intensive rather than a CPU-intensive task; and then the script is pushed to the Spark driver to process the chunks.

The main processing script `detect-breaks.R` takes several arguments as input.
First, option `-v` describes the vegetation index to use, one of EVI, NDMI and NIRV.
Next, option `-o` defines the harmonic order (usually set to 3).
Lastly, option `-t` defines the tile ID to process.
Tile IDs are in Proba-V grid, e.g. X16Y06.
Optional arguments include `-f` for the break detection method selection (`bfast` vs `bfastmonitor`), `-m` for the multithreading method (`spark`, `none` and `foreach`, the latter for doing local multithreading), and `--crop-only` to perform cropping but not change detection (used when doing cropping locally).
The input and output directories are preset in the script.

The script also depends on packages that implement BFAST in R.
In order to avoid every executor downloading the same packages for processing every chunk, and to make sure the package versions are well-known, the packages were preset and saved into a publicly-accessible directory within Terrascope.
The script first loads packages from that directory, if found.
The package that provides BFAST functionality is a [development version of BFAST](https://github.com/GreatEmerald/bfast) that is optimised for speed by moving parts of the R code into C++.

For convenience, the `process-tile.sh` script is further wrapped by `src/process-tiles.sh`, which runs `process-tile.sh` over all Proba-V tile IDs over Africa.
This wrapper also skips any tiles that already have been processed.
