# This app will allow a user to explore their mapped RNA seq reads across transcripts
# or positions of interest. 

library(shiny)
library(Rsamtools)

shinyServer(function(input, output) {
  
  # Get inputs
  #req(input$bam_file1)
  #req(input$gff_file)
  
  #bam1 <- scanBam(input$bam_file1)
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
})
