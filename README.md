# Model selection with ABC

Instructions, softwares and files for demographic model selection through Approximate Bayesian Computation for the manuscript [Late Pleistocene climate change shapes population divergence of an Atlantic Forest passerine: a model-based phylogeographic hypothesis test](https://link.springer.com/article/10.1007%2Fs10336-019-01650-1)

## Steps

1. Calculate summary statistics for the observed dataset
2. Simulate 1 million genetic datasets for each demographic scenario to be tested.
3. Calculate summary statistics for each simulated dataset in each demographic scenario.
4. Use provided R script to implement model selection using abc package (Csillery et al., 2010). For more information, check [this vignette](https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf)

## Data provided

* Arlsumstat files for Linux (as used in the manuscript) - file *arlsumstat_linux.zip*
* Scripts for the simulations in fastsimcoal of the 7 scenarios tested in the manuscript, for the three molecular markers used (cytb, g3pdh and tgfb2) - file *scen_fsimcoal_scripts.zip*
* Summary statistics calculated for each scenario and each marker - file *simulations_arlsumstat_results.zip* stored in Dropbox (click [here]() to access)
* Summary statistics calculated for each marker and the mean, for observed data (file *observed.zip*)
* ssdefs file defining the summary statistics to be calculated by arlsumstat (file *ssdefs.txt*)
* Script in R to perform the organization of the data outputted by arlsumstat, the calculation of means, the PCA pre-evaluation and the model selection in the abc package (file *script_abc.R*). More info [here](https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf).

## A few notes

* Simulations were made on fastsimcoal separately for each scenario and each molecular marker. .tpl and .est files for each scenario and each marker are in the file *scen_fsimcoal_scripts.zip*. For more information on how to create your scenarios, check the [fastSimcoal manual](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf).
* [Arlsumstat](http://cmpg.unibe.ch/software/arlequin35/Arl35Downloads.html) was used to calculate summary statistics. Since this software works on .arp Arlequin files, observed alignments for each molecular marker was loaded into DnaSP and converted to .arp files. fastSimcoal output files are already in .arp format.
* The final summary statistic (SS) used (both for observed and simulated values) was the mean of the SSs for the molecular markers (i.e., we had three molecular markers and we calculated the mean value of each SS for each marker).
* The final data format used in the abc R package is: 1) a vector with all the SS calculated for the observed data (the mean of each SS for each marker); 2) a data frame with simulated SS calculated from the simulated scenarios (the mean value of each SS for each marker, in each simulated dataset). In this data frame, columns must represent SS and rows must represent each simulated dataset; 3) a vector with the name of the scenario for each row of the data frame of simulated SSs. This vector must have the same number of rows of the data frame of simulated SSs (in our case, 1 million datasets for each of 7 scenarios, a total of 7 million datasets). 
