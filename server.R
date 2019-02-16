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
  
  # TODO add reactive for ispaired checkbox, I don't think it actually updates when 
  # checked off?
  
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
  
  # Get chromosome options from gff
  chroms <- reactive({
    req(input$gff_file)
    seqlevels(txDb())
  })
  
    # TODO add parsing the txDb for a list of features, ranges
    # TODO select a random one within range
  })
  
  # Update the chromosome select widget with chromosomes
  output$chromosomes <- renderUI({
    selectInput("chrom", label = h5("Select chromosome"),
                choices = chroms())
  })
  
  # TODO Update the search transcript widget with searchable list
  # TODO add Typeahead text input from shinysky
  
  # Update the main plot window
  output$read_gviz_plot <- renderPlot({
      bam_process()
      Gviz::plotTracks(list(axisTrack, annot_track(), bam_track()), type='hist',
                       from = 9989752, to = 9994025)
                       chromosome = input$chrom,
  })
  
  # TODO add ability to export plot 
  
})
