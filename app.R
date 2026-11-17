# Interactive version of the differential expression results.
# Run from the repo root: Rscript -e 'shiny::runApp()'
# Needs shiny and ggplot2: install.packages(c("shiny", "ggplot2"))

library(shiny)
library(ggplot2)

de <- read.csv("results/tables/DE_full_results.csv")
de$neg_log_p <- -log10(de$padj)
de$direction <- with(de, ifelse(
    is.na(padj) | padj >= 0.05 | abs(log2FoldChange) < 1,
    "not significant",
    ifelse(log2FoldChange > 0, "up in tumor", "down in tumor")))

top_up <- read.csv("results/tables/DE_top20_upregulated.csv")
top_dn <- read.csv("results/tables/DE_top20_downregulated.csv")
go_up <- head(read.csv("results/tables/enrichment_up_GO_BP.csv"), 15)
go_dn <- head(read.csv("results/tables/enrichment_down_GO_BP.csv"), 15)
pca <- read.csv("results/tables/pca_coordinates.csv")

# The top tables use gene symbols; the CSVs hold them as plain ordered lists.
up_symbols <- read.csv("results/tables/up_gene_symbols.csv")$symbol
dn_symbols <- read.csv("results/tables/down_gene_symbols.csv")$symbol

plot_volcano <- function(show_sig_only) {
    d <- de
    if (show_sig_only) d <- d[d$direction != "not significant", ]
    d$direction <- factor(d$direction,
                          levels = c("not significant", "up in tumor", "down in tumor"))
    ggplot(d, aes(x = log2FoldChange, y = neg_log_p, colour = direction)) +
        geom_point(size = 0.6, alpha = 0.5) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
        scale_colour_manual(values = c("grey70", "#d6604d", "#4a7bb7"), drop = FALSE) +
        labs(x = "log2 fold change (tumor vs normal)", y = "-log10 adjusted p",
             colour = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top")
}

plot_enrichment <- function(df) {
    df$term <- substr(df$Term, 1, 60)
    ggplot(df, aes(x = -log10(`Adjusted.P.value`), y = reorder(term, -log10(`Adjusted.P.value`)))) +
        geom_col(fill = "#4C72B0") +
        labs(x = "-log10 adjusted p", y = NULL) +
        theme_minimal(base_size = 11)
}

ui <- fluidPage(
    titlePanel("TCGA-BRCA tumor vs normal"),
    p("Differential expression on 25 matched tumor/normal pairs (pydeseq2, paired design).",
      "All numbers come from the CSVs in results/tables."),
    checkboxInput("sig_only", "Show only significant genes", FALSE),
    plotOutput("volcano", height = "520px"),
    fluidRow(
        column(6, h4("Top 20 up in tumor"), tableOutput("top_up")),
        column(6, h4("Top 20 down in tumor"), tableOutput("top_dn"))
    ),
    h4("GO Biological Process enrichment (top 15)"),
    fluidRow(
        column(6, plotOutput("go_up", height = "500px")),
        column(6, plotOutput("go_dn", height = "500px"))
    ),
    h4("PCA, tumor vs normal"),
    plotOutput("pca", height = "450px")
)

server <- function(input, output) {
    output$volcano <- renderPlot(plot_volcano(input$sig_only))

    fmt_top <- function(df, sym) {
        data.frame(
            gene = sym[seq_len(nrow(df))],
            log2FC = sprintf("%.2f", df$log2FoldChange),
            padj = format(df$padj, digits = 2),
            base_mean = format(round(df$baseMean), big.mark = ",")
        )
    }

    output$top_up <- renderTable(fmt_top(top_up, up_symbols), rownames = FALSE)
    output$top_dn <- renderTable(fmt_top(top_dn, dn_symbols), rownames = FALSE)
    output$go_up <- renderPlot(plot_enrichment(go_up))
    output$go_dn <- renderPlot(plot_enrichment(go_dn))

    output$pca <- renderPlot({
        ggplot(pca, aes(x = PC1, y = PC2, colour = condition)) +
            geom_point(size = 3) +
            geom_line(aes(group = patient_id), colour = "grey80", linewidth = 0.4) +
            scale_colour_manual(values = c(Normal = "#4a7bb7", Tumor = "#d6604d")) +
            labs(colour = NULL) +
            theme_minimal(base_size = 13) +
            theme(legend.position = "top")
    })
}

shinyApp(ui, server)
