# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(Gviz)

options(ucscChromosomeNames=FALSE)
options(shiny.maxRequestSize=500*1024^2)

# Always will have the generic axis, no need to re-render.
axisTrack <- GenomeAxisTrack()

shinyServer(function(input, output) {
  # Get input alignment file to plot (bam format)
  # Needs to be indexed to run
  bam_process <- reactive({
    req(input$bam_file1)
    indexBam(input$bam_file1$datapath)
  })
  
  # Create data track with bam alignments
  bam_track <- reactive({
    req(input$bam_file1)
    Gviz::DataTrack(range = input$bam_file1$datapath, stream=TRUE, legend=TRUE, 
                  isPaired = input$ispaired)
  })
  
  # Get input annotations (gtf/gff format)
  txDb <- reactive({
    req(input$gff_file)
    GenomicFeatures::makeTxDbFromGFF(file = input$gff_file$datapath, format = "auto")
  })
  
  # Create annotation track from gff
  annot_track <- reactive({
    Gviz::GeneRegionTrack(txDb(), collapseTranscripts = "meta")
  }) 
  
  })
  
  # Update the main plot window
  output$read_gviz_plot <- renderPlot({
      bam_process()
      Gviz::plotTracks(list(axisTrack, annot_track(), bam_track()), type='hist',
                       chromosome = "CHROMOSOME_I",
                       from = 9989752, to = 9994025)
  })
  
})
