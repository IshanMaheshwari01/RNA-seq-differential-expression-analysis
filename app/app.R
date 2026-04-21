library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(plotly)

# Load data
de_results <- read.csv("../results/tables/DE_results_full.csv")
sig_genes <- read.csv("../results/tables/DE_results_significant.csv")
go_up <- read.csv("../results/tables/GO_enrichment_upregulated.csv")
go_down <- read.csv("../results/tables/GO_enrichment_downregulated.csv")

# UI
ui <- dashboardPage(
  skin = "purple",
  
  dashboardHeader(title = "RNA-Seq DE Explorer", titleWidth = 300),
  
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Gene Explorer", tabName = "explorer", icon = icon("search")),
      menuItem("Volcano", tabName = "volcano", icon = icon("chart")),
      menuItem("Pathways", tabName = "pathways", icon = icon("project-diagram")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .skin-purple .main-header .logo {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        }
        .skin-purple .main-header .navbar {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        }
        .content-wrapper { background-color: #f4f6f9; }
        .box { border-top: 3px solid #667eea; }
      "))
    ),
    
    tabItems(
      
      # OVERVIEW
      tabItem(tabName = "overview",
              h2("RNA-Seq Differential Expression Analysis"),
              p("Interactive explorer for tumor vs normal tissue comparison"),
              
              fluidRow(
                valueBox(format(nrow(de_results), big.mark = ","), "Total Genes", color = "purple", width = 3),
                valueBox(format(nrow(sig_genes), big.mark = ","), "DE Genes", color = "red", width = 3),
                valueBox(format(sum(sig_genes$direction == "Upregulated"), big.mark = ","), "Upregulated", color = "orange", width = 3),
                valueBox(format(sum(sig_genes$direction == "Downregulated"), big.mark = ","), "Downregulated", color = "blue", width = 3)
              ),
              
              fluidRow(
                box(title = "Key Findings", status = "primary", solidHeader = TRUE, width = 6,
                    h4("3,296 Differentially Expressed Genes"),
                    tags$ul(
                      tags$li("Fold changes: 0.04x to 28.1x"),
                      tags$li("Min p-value: 1.36×10⁻¹⁴"),
                      tags$li("Up: Cell cycle, DNA replication"),
                      tags$li("Down: Immune response, defense")
                    )
                ),
                
                box(title = "Cancer Hallmarks", status = "danger", solidHeader = TRUE, width = 6,
                    tags$ul(
                      tags$li("Sustained proliferation"),
                      tags$li("Immune evasion"),
                      tags$li("Evading growth suppressors"),
                      tags$li("Replicative immortality")
                    )
                )
              )
      ),
      
      # EXPLORER
      tabItem(tabName = "explorer",
              h2("Gene Explorer"),
              
              fluidRow(
                box(title = "Filter Options", status = "primary", solidHeader = TRUE, width = 12,
                    column(3, sliderInput("fc_cutoff", "Log2 FC", min = 0, max = 5, value = 1)),
                    column(3, sliderInput("padj_cutoff", "-Log10(P-adj)", min = 0, max = 15, value = 1.3)),
                    column(3, selectInput("direction", "Direction", c("All", "Upregulated", "Downregulated"))),
                    column(3, textInput("gene_name", "Search Gene"))
                )
              ),
              
              fluidRow(
                box(title = "Results", status = "success", solidHeader = TRUE, width = 12,
                    DTOutput("gene_table")
                )
              )
      ),
      
      # VOLCANO
      tabItem(tabName = "volcano",
              h2("Volcano Plot"),
              
              fluidRow(
                box(title = "Controls", status = "primary", solidHeader = TRUE, width = 12,
                    column(4, sliderInput("vol_fc", "FC Threshold", min = 0.5, max = 3, value = 1)),
                    column(4, sliderInput("vol_p", "-Log10(P-adj)", min = 1, max = 15, value = 1.3)),
                    column(4, p("Hover over points for gene names"))
                )
              ),
              
              fluidRow(
                box(title = "Interactive Plot", status = "warning", solidHeader = TRUE, width = 12,
                    plotlyOutput("volcano_plot", height = "600px")
                )
              )
      ),
      
      # PATHWAYS
      tabItem(tabName = "pathways",
              h2("Pathway Enrichment"),
              
              fluidRow(
                box(title = "Upregulated", status = "danger", solidHeader = TRUE, width = 6,
                    DTOutput("pathways_up")
                ),
                box(title = "Downregulated", status = "info", solidHeader = TRUE, width = 6,
                    DTOutput("pathways_down")
                )
              )
      ),
      
      # DOWNLOAD
      tabItem(tabName = "download",
              h2("Download Data"),
              
              fluidRow(
                box(title = "Results", status = "primary", solidHeader = TRUE, width = 12,
                    h4("Download CSV Files:"),
                    downloadButton("dl_all", "All Genes"),
                    downloadButton("dl_sig", "Significant Only"),
                    downloadButton("dl_up", "Upregulated"),
                    downloadButton("dl_down", "Downregulated"),
                    
                    br(), br(),
                    h4("View on GitHub:"),
                    tags$a(href = "https://github.com/IshanMaheshwari01/RNA-seq-differential-expression-analysis", 
                           target = "_blank", class = "btn btn-success", "GitHub Repository")
                )
              )
      ),
      
      # ABOUT
      tabItem(tabName = "about",
              h2("About"),
              
              fluidRow(
                box(title = "Project", status = "primary", solidHeader = TRUE, width = 6,
                    h4("RNA-Seq Analysis"),
                    p("Tumor vs Normal tissue comparison"),
                    h4("Features:"),
                    tags$ul(
                      tags$li("Interactive exploration"),
                      tags$li("Dynamic plots"),
                      tags$li("Data downloads"),
                      tags$li("Pathway analysis")
                    )
                ),
                
                box(title = "Author", status = "success", solidHeader = TRUE, width = 6,
                    h4("Ishan Maheshwari"),
                    p("Bioinformatics"),
                    p("Biostatistics"),
                    p("Computational Biology"),
                    p("Health Data Science"),
                    p("Genetics Data Science"),
                    hr(),
                    p(icon("envelope"), " ishanmaheshwari02@gmail.com"),
                    p(icon("linkedin"), tags$a(href = "https://www.linkedin.com/in/ishanmaheshwari2001", target = "_blank", "LinkedIn")),
                    p(icon("github"), tags$a(href = "https://github.com/IshanMaheshwari01", target = "_blank", "GitHub"))
                )
              )
      )
    )
  )
)

# SERVER
server <- function(input, output, session) {
  
  filtered_genes <- reactive({
    data <- de_results %>%
      filter(!is.na(padj)) %>%
      mutate(abs_lfc = abs(log2FoldChange), neg_log_padj = -log10(padj))
    
    data <- data %>%
      filter(abs_lfc >= input$fc_cutoff, neg_log_padj >= input$padj_cutoff)
    
    if (input$direction != "All") {
      data <- data %>% filter(direction == input$direction)
    }
    
    if (input$gene_name != "") {
      data <- data %>% filter(grepl(input$gene_name, gene, ignore.case = TRUE))
    }
    
    data
  })
  
  output$gene_table <- renderDT({
    filtered_genes() %>%
      select(gene, log2FoldChange, padj, direction) %>%
      mutate(padj = format(padj, scientific = TRUE, digits = 2)) %>%
      datatable(options = list(pageLength = 25), colnames = c("Gene", "Log2 FC", "P-adj", "Direction"))
  })
  
  output$volcano_plot <- renderPlotly({
    plot_data <- de_results %>%
      filter(!is.na(padj)) %>%
      mutate(neg_log_p = -log10(padj),
             color = case_when(
               abs(log2FoldChange) >= input$vol_fc & neg_log_p >= input$vol_p & log2FoldChange > 0 ~ "Up",
               abs(log2FoldChange) >= input$vol_fc & neg_log_p >= input$vol_p & log2FoldChange < 0 ~ "Down",
               TRUE ~ "NS"))
    
    plot_ly(plot_data, x = ~log2FoldChange, y = ~neg_log_p, color = ~color,
            colors = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#95a5a6"),
            text = ~gene, type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.6)) %>%
      layout(title = "Volcano Plot", xaxis = list(title = "Log2 FC"), 
             yaxis = list(title = "-Log10(P-adj)"), hovermode = "closest")
  })
  
  output$pathways_up <- renderDT({
    go_up %>% select(Description, GeneRatio, p.adjust, Count) %>%
      mutate(p.adjust = format(p.adjust, scientific = TRUE, digits = 2)) %>%
      datatable(options = list(pageLength = 10), colnames = c("Pathway", "Ratio", "P-adj", "Count"))
  })
  
  output$pathways_down <- renderDT({
    go_down %>% select(Description, GeneRatio, p.adjust, Count) %>%
      mutate(p.adjust = format(p.adjust, scientific = TRUE, digits = 2)) %>%
      datatable(options = list(pageLength = 10), colnames = c("Pathway", "Ratio", "P-adj", "Count"))
  })
  
  output$dl_all <- downloadHandler(
    filename = "DE_all.csv",
    content = function(file) write.csv(de_results, file, row.names = FALSE)
  )
  
  output$dl_sig <- downloadHandler(
    filename = "DE_significant.csv",
    content = function(file) write.csv(sig_genes, file, row.names = FALSE)
  )
  
  output$dl_up <- downloadHandler(
    filename = "DE_upregulated.csv",
    content = function(file) write.csv(sig_genes %>% filter(direction == "Upregulated"), file, row.names = FALSE)
  )
  
  output$dl_down <- downloadHandler(
    filename = "DE_downregulated.csv",
    content = function(file) write.csv(sig_genes %>% filter(direction == "Downregulated"), file, row.names = FALSE)
  )
}

shinyApp(ui = ui, server = server)