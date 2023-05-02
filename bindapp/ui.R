library(shiny)
library(ggplot2)
library(shinythemes)
library(RVenn)

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
                   textInput("shell_run", "BindCompare Directory Filepath"),
                   textInput("base_bed_path", "Base Bed File Path"),
                   textInput("olay_bed_path", "Overlayed Bed File Path"),
                   numericInput("scope", "Scope", value = 1000),
                   textInput("sample_name", "Sample Name"),
                   textInput("out_dir", "Output Folder (Must Already Exist!)"),
                   textInput("genes_gtf_filepath", "Genes GTF File Path"),
                   textInput("genome_fa_filepath", "Genome fasta File Path")),
          actionButton("run_command_btn", strong("Run BindCompare!"), class="btn btn-primary"),
          ),
          
          # Main panel for displaying outputs
          mainPanel(
            # h1, h2, h3, h4, h5, h6
            # Placeholder text
            helpText("Please enter the full filepaths into the boxes to run BindCompare!",
                     "A full explanation can be found on the GitHub Repo found below."),
            htmlOutput("status"), br(),
            textOutput("overlap"),
            imageOutput("overlapprofile"), br(), br(), br(), br(), br(),
            imageOutput("bargraph"), br(), br(), br(),
            htmlOutput("output"), br(),
          )
          
        )
      )
    ),
    tabPanel(
      "Compare Unqiue Experiments",
      icon = icon("table"), # Add a different icon to the second tab label
      fluidPage(
        # Insert input elements for Tab 2 here
        titlePanel("Compare Gene Lists from Two BindCompare Experiments"),
        sidebarLayout(
          sidebarPanel(
            textInput("gene_list_1", "Gene List 1", ""),
            textInput("gene_list_2", "Gene List 2", ""),
            actionButton("submit_button", "Submit")
          ),
          mainPanel(
            helpText("The below UI will fill when you submit the Gene Lists.",
                     "Please copy and paste the Gene Lists printed in the summary file!"),
            h3("Jaccard Similarity Results"),
            verbatimTextOutput("jaccard_output"),
            h3("Size-Biased Venn Diagram"),
            plotOutput("venn_plot"),
            h3("Gene Lists"),
            h4("Counts in Each Category"),
            htmlOutput("gene_list_output"),
            h4("Genes Only in List 1"),
            verbatimTextOutput("only_list_1"),
            h4("Genes Only in List 2"),
            verbatimTextOutput("only_list_2"),
            h4("Genes in Both Lists"),
            verbatimTextOutput("both_lists"),
          )
        )
      )
    ),
     tabPanel(
      "How to BindCompare",
      icon = icon("table"), # Add a different icon to the second tab label
      fluidPage(
        # Insert input elements for Tab 2 here
        titlePanel("Help Manual for BindCompare"),
        mainPanel(
            helpText("For a full how-to, please see the GitHub README!"),
        )
      )
    ),
    footer = list(
      # Customize the footer with any relevant information or links
      p("Copyright © 2023"),
      a("Github README!", href = "https://github.com/pranavmahabs/bindcompare/blob/main/README.md"),
      a("Created by Pranav Mahableshwarkar and the Larschan Lab", href = "larschanlab.com")
    )
  )
)