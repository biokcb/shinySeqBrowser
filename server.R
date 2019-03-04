# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(Gviz)
library(Cairo)

options(ucscChromosomeNames=FALSE)
options(shiny.maxRequestSize=5000*1024^2)

shinyServer(function(input, output, session) {
  # Get input alignment file & index
  bam_process <- reactive({
    req(input$bam_file1)
    indexBam(input$bam_file1$datapath)
  })

  # Create data track with bam alignments
  bam_track <- reactive({
    req(input$bam_file1)
    bam_process()
    Gviz::DataTrack(range = input$bam_file1$datapath, stream=TRUE, legend=TRUE, 
                  isPaired =  input$ispaired, fill = input$readcol, col = input$readcol)
  })
  
  # Get input annotations (gtf/gff format)
  txDb <- reactive({
    req(input$gff_file)
    GenomicFeatures::makeTxDbFromGFF(file = input$gff_file$datapath, format = "auto")
  })
  
  # Create annotation track from gff
  annot_track <- reactive({
    Gviz::GeneRegionTrack(txDb(), collapseTranscripts = 'meta', fill = input$annotcol, showId=TRUE)
  }) 
  
  # Get chromosome options from gff
  chroms <- reactive({
    seqlevels(txDb())
  })
  
  # Create a list of genes from txDb
  gene_list <- reactive({
    req(input$gff_file)
    names(genes(txDb()))
  })
  
  # Get the chromosome of selected gene
  gene_chrom <- reactive({
    seqnames(genes(txDb())[gene_list()==input$gene])
  })
  
  # Update the start & end positions based on gene selection
  observe({
    gene_pos <- ranges(genes(txDb())[gene_list()==input$gene])
    gene_start <- start(gene_pos) - 500
    gene_end <- end(gene_pos) + 500
    updateNumericInput(session, "start", value=gene_start)
    updateNumericInput(session, "end", value=gene_end)
  })
  
  # Update the chromosome select widget with chromosome options + selected gene chromosome
  output$chromosomes <- renderUI({
    selectInput("chrom", label = h5("Select chromosome"),
                choices = chroms(), selected = gene_chrom())
  })

  # Generate the gene search bar from list
  output$gene_search <- renderUI({
    selectizeInput("gene", label = h5("Search transcripts"), selected = NULL, multiple=FALSE, 
                   options=list(placeholder = "Enter gene name...", maxItems = 1, maxOptions = 10),
                   choices = gene_list())
  })

  # Update the main plot window
  output$read_gviz_plot <- renderPlot({
      Gviz::plotTracks(list(GenomeAxisTrack(), annot_track(), bam_track()), type='hist',
                       chromosome = input$chrom, from = input$start, to = input$end,
                       background.panel = input$bgcol, background.title = input$pancol,
                       window = -1, col.histogram = input$readcol)
  })
  
  # Download the plot as a PDF
  output$save <- downloadHandler(
    filename = function(){
      paste(Sys.Date(), "_read_plot", ".pdf", sep="")
    },
    content = function(filename){
      cairo_pdf(filename)
      Gviz::plotTracks(list(GenomeAxisTrack(), bam_track(), annot_track()), type='hist',
                       chromosome = input$chrom, from = input$start, to = input$end,
                       background.panel = input$bgcol, background.title = input$pancol,
                       window = -1, col.histogram = input$readcol)
      dev.off()
    }
  )
  
})
