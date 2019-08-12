% Land cover change detection algorithm theoretical basis
% Dainius Masiliūnas
% 2019-08-12

# Land cover change detection

## Overview

For producing yearly map updates, a change detection algorithm is needed to assess whether a given pixel is temporally stable or not.
If it is stable, then it is desirable to keep the previous land cover classification, in order to avoid spurious changes that comes from minor changes in the input data for the classifier.

In order to do so, we make use of the [BFAST algorithm](http://dx.doi.org/10.1016/j.rse.2009.08.014) to detect breaks in the time series of vegetation indices.
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
Optional arguments include `-m` for the method selection (BFAST vs BFAST Monitor), and `--crop-only` to perform cropping but not change detection (used when doing cropping locally).
The input and output directories are preset in the script.

The script also depends on packages that implement BFAST in R.
In order to avoid every executor downloading the same packages for processing every chunk, and to make sure the package versions are well-known, the packages were preset and saved into a publicly-accessible directory within Terrascope.
The script first loads packages from that directory, if found.
The package that provides BFAST functionality is a [development version of BFAST](https://github.com/GreatEmerald/bfast) that is optimised for speed by moving parts of the R code into C++.

For convenience, the `process-tile.sh` script is further wrapped by `src/process-tiles.sh`, which runs `process-tile.sh` over all Proba-V tile IDs over Africa.
This wrapper also skips any tiles that already have been processed.
