# rawrr_mrm_extraction
Repository for a small R/rawrr script to extract MRM data directly from Thermo RAW files

The main script is rawrr_mrm_extraction.R, which can be run after the appropiate R packages and extras have been installed.

To locally install R/R packages/extras, follow these instructions:

After installing R, install the rawrr package via the Bioconductor commands (on R shell):
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(version = "3.17")
```

The rawrr package requires specific DLLs and an exe to extract information from the binary Thermo RAW files. Execute the following in R:

```
if (isFALSE(rawrr::.checkDllInMonoPath())){
  rawrr::installRawFileReaderDLLs()
}
```
It will ask you to accept Terms and Conditions by entering "y" and the last lines of the output should look like

```
##                   ThermoFisher.CommonCore.Data.dll 
##                                                  0 
## ThermoFisher.CommonCore.MassPrecisionEstimator.dll 
##                                                  0 
##          ThermoFisher.CommonCore.RawFileReader.dll 
##                                                  0
```
Then, install the required exe application using 
```
rawrr::installRawrrExe()
```
Additional packages can be installed using:
```
install.packages(tidyverse)
```
or your preferred method.


Change the folder path in the R script accordingly before using it. The time it takes to process a batch of files depends on the number of files and the number of scans. The functions in the rawrr package are built thinking on more complex Orbitrap data, but small modifications can be performed to extract any kind of MS1/MS2 experiments from Thermo RAW format. Here we used it for MRM experiments. 

The output, which is a master table (I called it hypertable), contains the scan data for each file and each precursor/fragments defined in the MS method. For each scan, the counts and the time are available, therefore the desired peak curve is given by joining the scan intensities with a line over the whole acquisition time.

We haven't used any sophiscated method to select an specific time window to define the peak (peak picking), we just considered the whole acquisition time as our window of interest and integrate (area under the curve) the peak over that interval.

Before the integration, we performed baseline correction using a common algorithm, Asymmetric Least Squares. The implementation I found is in Python, so I create a Python script for it, that also performs the integration part at the end using a trapezoid rule. The input of the Python script is the master hypertable from the R script.

The Python script requires the [pybaselines](https://pybaselines.readthedocs.io/en/latest/) package. You can install it using pip or conda.

