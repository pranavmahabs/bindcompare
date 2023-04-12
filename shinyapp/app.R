#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(shinythemes)

freqpoly <- function(x1, x2, binwidth = 0.1, xlim = c(-3, 3)) {
  df <- data.frame(
    x = c(x1, x2),
    g = c(rep("x1", length(x1)), rep("x2", length(x2)))
  )
  
  ggplot(df, aes(x, colour = g)) +
    geom_freqpoly(binwidth = binwidth, size = 1) +
    coord_cartesian(xlim = xlim)
}

t_test <- function(x1, x2) {
  test <- t.test(x1, x2)
  
  # use sprintf() to format t.test() results compactly
  sprintf(
    "p value: %0.3f\n[%0.2f, %0.2f]",
    test$p.value, test$conf.int[1], test$conf.int[2]
  )
}

ui <- shinyUI(
  navbarPage(
    title = "BindCompare",
    theme = shinytheme("flatly"), # Use a pre-defined theme to style the app
    tabPanel(
      "Compare Two BED Files",
      icon = icon("chart-column"), # Add an icon to the tab label
      fluidPage(
        # Insert input elements for Tab 1 here
        sidebarLayout(
          sidebarPanel(
            
            # First tab input options
          tabPanel("Tab 1",
                     textInput("base_bed_path", "Base Bed File Path"),
                     textInput("olay_bed_path", "Overlayed Bed File Path"),
                     numericInput("scope", "Scope", value = 1000),
                     textInput("sample_name", "Sample Name"),
                     textInput("out_dir", "Output Directory"),
                     textInput("genes_gtf_filepath", "Genes GTF File Path"),
                     textInput("genome_fa_filepath", "Genome fasta File Path"))
          ),
          
          # Main panel for displaying outputs
          mainPanel(
            
            # Placeholder text
            h3("This is where the output will be displayed.")
            
          )
          
        )
      )
    ),
    tabPanel(
      "Compare Unqiue Experiments",
      icon = icon("table"), # Add a different icon to the second tab label
      fluidPage(
        # Insert input elements for Tab 2 here
      )
    ),
    footer = list(
      # Customize the footer with any relevant information or links
      p("  Copyright © 2023"),
      a("  Privacy Policy", href = "https://example.com/privacy"),
      a("  Terms of Service", href = "https://example.com/terms")
    )
  )
)

server <- function(input, output, session) {
}

# Run the application 
shinyApp(ui = ui, server = server)
