# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 
library(shiny)
library(shinycssloaders)
library(colourpicker)

shinyUI(fluidPage(
  ### Basic app elements
  titlePanel("shinySeqBrowser"),
  
  ### User-modifiable options
  sidebarLayout(
    sidebarPanel(
      h4("File inputs"),
      
      # Upload the bam file
      fileInput("bam_file1", "Choose a BAM File",
                multiple = FALSE,
                accept = c(".bam")),
      
      # indicate paired-end
      checkboxInput("ispaired", label="Paired-end data", value=FALSE),
        
      # Upload the GFF3 file
      fileInput("gff_file", "Choose a GFF3 File",
                 multiple = FALSE,
                accept = c(".gff3")),
        
      # Add horizontal line to split from upload section
      hr(),
      h4("Search parameters"),
        
      # Gene input
      textInput("gene", label = h5("Search transcripts"), placeholder = "Enter gene name..." ),
        
      # Chromosome input
      uiOutput("chromosomes"),
        
      # position input
      numericInput("start", label = h5("Enter start position"), value=1),
      numericInput("end", label = h5("Enter end position"), value=100000),

      hr(),
      downloadButton("save", "Download plot as PDF")
    ),
    # Data visualizations here!
    mainPanel(
      # The main plot contains a distribution of aligned reads
      withSpinner(plotOutput("read_gviz_plot")),
      h4("Plot Options"),
      hr(),
        # Color palette choices
        colourInput("bgcol", "Select background color", "white"),
        colourInput("pancol", "Select panel color", "#51B9CF"),
        colourInput("annotcol", "Select gene track color", "#FDC010"),
        colourInput("readcol", "Select read color", "#333333")
    )
  )
))
