# Morphometrics practical
This repository is used to give an example of how to perform a Procrustes analysis (PA) with the R package `geomorph`.

The data set used is a subset of the [data](http://datadryad.org/resource/doi:10.5061/dryad.nr210) provided by [Ivanović A, Arntzen JW](https://academic.oup.com/biolinnean/article-lookup/doi/10.1111/bij.12314), which consists of one specimen of each species randomly sampled. This subset can be found in the file *"Triturus_and_Calotriton_lmk_reduced.csv"*

The R script includes a series of commands to:
* Prepare the environment
* Load and prepare the data set *"Triturus_and_Calotriton_lmk_reduced.csv"*
* Perform the PA with the R function `geomorph::gpagen`
* Plot the superimposition of landmark points before and after the PA.

A detailed explanation of the usage of PA with morphological data can be found at 
