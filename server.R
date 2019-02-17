# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(Gviz)
library(Cairo)

options(ucscChromosomeNames=FALSE)
options(shiny.maxRequestSize=500*1024^2)

shinyServer(function(input, output) {
  # Get input alignment file & index
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
    seqlevels(txDb())
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
                       chromosome = input$chrom,
                       from = input$start, to = input$end)
      Gviz::plotTracks(list(GenomeAxisTrack(), annot_track(), bam_track()), type='hist',
  })
  
  # Download the plot as a PDF
  output$save <- downloadHandler(
    filename = function(){
      paste(Sys.Date(), "_read_plot", ".pdf", sep="")
    },
    content = function(filename){
      cairo_pdf(filename)
      bam_process()
      Gviz::plotTracks(list(GenomeAxisTrack(), annot_track(), bam_track()), type='hist',
                       chromosome = input$chrom, from = input$start, to = input$end,
                       background.panel = input$bgcol, background.title = input$pancol)
      dev.off()
    }
  )
  
})
