# shinySeqBrowser: An R shiny web application for exploring read data across a genome

`shinySeqBrowser` is similar to genome browser tools (such as IGV), mainly built around the R package `Gviz`. I created this project to learn how to build R shiny apps, how to use the Gviz package, all while providing something potentially useful for my lab. The main purpose is to be able to explore RNA-seq data across specific transcripts and easily export these plots for publications. This tool is still under development and feedback/issue reporting is welcome. 

## Features
* Analyze read coverage of your BAM files
* Search transcripts by name, select regions by position/chromosome
* Customize styling and export plots for publication

## Requirements
In order to use this app you need to first install `R` from [cran-r](https://cran.r-project.org/). I recommend that you also install [RStudio](https://www.rstudio.com/) to make your life easier if you plan on using R long-term. 

Install the following dependencies from within R:
```
# Shiny
install.packages("shiny")
install.packages("shinythemes")

# Bioconductor tools
source("http://bioconductor.org/biocLite.R")
biocLite("Gviz")
biocLite("Rsamtools")
biocLite("GenomicFeatures")
biocLite("GenomicRanges")
```

## Running the app
Then the easiest way to run the app is to run this from within R/RStudio:
```
library(shiny)
runGitHub("shinySeqBrowser", "biokcb")
```

## Usage
TODO: Add some examples and images here once it's up and running well enough. 

## Contributing
I do not have contributing guidelines for this repository yet, but if you would like to contribute to the code feel free to fork the repository and create a new branch for your feature, bug fix, etc then create a pull request. Collaborative draft pull requests also welome. Create an issue if you have run into a bug or have a suggestion. Thanks!

## License
Due to R shiny being licensed GPLv3, I believe this app should also be licensed that way, but it is unclear to me.. clarifications welcome. :) See [LICENSE](LICENSE)
