# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(Gviz)

options(ucscChromosomeNames=FALSE)

setwd('~/Documents/GitHub/shinySeqBrowser/')
# Implement an example case in wild type C. elegans
cel_n2_bam <- "./data/cel_example.bam"
indexBam("./data/cel_example.bam")
cel_ws230_gff <- "./data/cel_ws230_example.gff3"

cel_ws230_db <- GenomicFeatures::makeTxDbFromGFF(file = cel_ws230_gff, format = "gff3")

cel_gtrack <- Gviz::GeneRegionTrack(range=cel_ws230_db, chromosome = "CHROMOSOME_I", 
                      start = 9989752, end = 9994025, stacking = "dense",
                      geneSymbols = TRUE)
cel_altrack <- Gviz::AlignmentsTrack(range = cel_n2_bam, chromosome = "CHROMOSOME_I", 
                                     start = 9989752, end = 9994025, stacking = "dense",
                                     isPaired = TRUE)
  
shinyServer(function(input, output) {
  
  # Get inputs
  bam_upload <- reactive({
    bam1 <- scanBam(input$bam_file1)
    return(bam1)
  })
  
  gff_upload <- reactive({
    gff <- read.delim(input$gff_file)
    return(gff)
  })
  
  output$read_gviz_plot <- renderPlot({
    
      Gviz::plotTracks(c(cel_altrack, cel_gtrack), 
                       from = 9989752, to = 9994025, stacking = "dense")
    
  })
  
})
