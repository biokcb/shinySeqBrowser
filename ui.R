# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  
  # Basic app elements
  theme = shinytheme("slate"),
  titlePanel("Browse RNA-seq reads aligned across transcripts"),
  
  # User-modifiable options
  sidebarLayout(
    sidebarPanel(
        h4("File inputs"),
        # Upload the bam file
        fileInput("bam_file1", "Choose a BAM File",
                  multiple = FALSE,
                  accept = c(".bam")),
        
        checkboxInput("ispaired", label="Paired-end data", value=FALSE),
        
        # Upload the GFF3 file
        fileInput("gff_file", "Choose a GFF3 File",
                  multiple = FALSE,
                  accept = c(".gff3")),
        
        # Add horizontal line to split from upload section
        hr(),
        h4("Search inputs"),
        
        # Gene input
        textInput("gene", label = h5("Search transcripts"), placeholder = "Enter gene name..." ),
        
        # Chromosome input
        uiOutput("chromosomes"),
        
        # position input
        numericInput("start", label = h5("Enter start position"), value=0),
        numericInput("end", label = h5("Enter end position"), value=0),
        # Add horizontal line to split from search section
        hr(),
        h4("Plot Options"),
        
        # Color palette choices
        selectInput("colormap", label = h5("Select color palette"), 
                    choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
                    selected = 1),
        
        # Slider input to adjust number of bins
        sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 50,
                   value = 30)
    ),
    
    # Data visualizations here!
    mainPanel(
       # The main plot contains a distribution of aligned reads
       plotOutput("read_gviz_plot")
    )
  )
))
