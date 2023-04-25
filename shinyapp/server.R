library(shiny)
library(ggplot2)
library(shinythemes)
library(eulerr)

processFile <- function(filepath) {
  con = file(filepath, "r")
  outlines <- c()
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    outlines <- append(outlines, line)
  }
  close(con)
  return(outlines)
}

server <- function(input, output, session) {
    observeEvent(input$run_command_btn, {
      # Define the command and its arguments as separate variables
      shellfile <- reactive({input$shell_run})
      base_bed <- reactive({input$base_bed_path})
      olay_bed <- reactive({input$olay_bed_path})
      scope <- reactive({input$scope})
      name <- reactive({input$sample_name})
      outdir <- reactive({input$out_dir})
      gtf_path <- reactive(input$genes_gtf_filepath)
      genome_path <- reactive({input$genome_fa_filepath})
      
      # Run the command and capture both the output and the exit status
      bindcompare <- paste0(shellfile(), "/bindcompare.sh")
      cmd <- paste(c(bindcompare, base_bed(), olay_bed(), scope(), name(), outdir(), gtf_path(), genome_path()), collapse = " ")
      command_output <- system(cmd, intern = TRUE)
      summary_path <- paste0(outdir(), "/", name(), "_summary.txt")
      
      # Display the output and exit status in the app
      output$status <- renderUI(HTML("<b>Outputs Have Been Saved! Here is some of the Script Output (Full Results can be found in Outdir):</b>"))
      output$output <- renderUI({
        str2 <- paste(processFile(summary_path), collapse = '<br/>')
        str3 <- paste(command_output, collapse = '<br/>')
        HTML(paste(str2, str3, sep = '<br/>'))
      })
      output$overlap <- renderText("1) Binding Overlap Profile generated from your BED Files! 2) Number of Overlaps Found in Comparison to Total Number of Binding Sites.")
      overlaps <- paste0(outdir(), "/", name(), "_overlaps.png")
      output$overlapprofile <- renderImage({list(src=overlaps)}, deleteFile=FALSE)
      bargraph <- paste0(outdir(), "/", name(), "_bartotals.png")
      output$bargraph <- renderImage({list(src=bargraph)}, deleteFile=FALSE)
    })
  
  observeEvent(input$submit_button, {
    # Parse the input gene lists
    # Parse the input gene lists
    gene_list_1 <- strsplit(input$gene_list_1, ",")[[1]]
    gene_list_2 <- strsplit(input$gene_list_2, ",")[[1]]
    
    # Calculate the Jaccard similarity
    intersection <- length(intersect(gene_list_1, gene_list_2))
    union <- length(union(gene_list_1, gene_list_2))
    jaccard_similarity <- intersection / union
    
    # Output the Jaccard similarity and p-value results
    output$jaccard_output <- renderPrint({
      paste0("Jaccard similarity between Gene List 1 and Gene List 2: ", jaccard_similarity)
    })
    
   venn_list <- list(
      gene_list_1 = gene_list_1,
      gene_list_2 = gene_list_2
    )
    genes.venn <- euler(venn_list)
    
    # Render the Venn diagram
    output$venn_plot <- renderPlot({
      plot(genes.venn, quantities = TRUE)
    })

    # Categorize genes into lists
    only_list_1 <- setdiff(gene_list_1, gene_list_2)
    only_list_2 <- setdiff(gene_list_2, gene_list_1)
    both_lists <- intersect(gene_list_1, gene_list_2)
    all_outs <- c(paste0("<b>Number of genes only in Gene List 1: </b>", length(setdiff(gene_list_1, gene_list_2)),
                         "<b> Number of genes only in Gene List 2: </b>", length(setdiff(gene_list_2, gene_list_1)), sep = ' '),
          paste0("<b>Number of genes in both lists: </b>", length(intersect(gene_list_1, gene_list_2)),
                 "<b> Total number of unique genes: </b>", length(unique(c(gene_list_1, gene_list_2))), sep = ' '))
    output$gene_list_output <- renderUI({
        str1 <- paste(all_outs, collapse = '<br/>')
        HTML(str1)
    })
    output$only_list_1 <- renderPrint(paste(only_list_1, collapse = ' '))
    output$only_list_2 <- renderPrint(paste(only_list_2, collapse = ' '))
    output$both_lists <- renderPrint(paste(both_lists, collapse = ' '))

  })
}