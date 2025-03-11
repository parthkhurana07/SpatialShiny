# Libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)   # for combining ggplots
library(DT)          # for interactive tables
library(dplyr)

setwd("/Users/parthkhurana/Desktop/Projects/SpatialShiny")
data_dir <- '/Users/parthkhurana/Desktop/Projects/SpatialShiny/data'

spatial_seurat <- Load10X_Spatial(data.dir = data_dir)

spatial_seurat <- SCTransform(spatial_seurat, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:10)

# Create a vector of gene choices for selection.
gene_choices <- rownames(spatial_seurat)

# Get cluster assignments from metadata (if available)
cluster_choices <- sort(unique(as.character(spatial_seurat$seurat_clusters)))

ui <- navbarPage("Spatial Transcriptomics Dashboard",
                 tabPanel("Spatial Feature Plot",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("gene", "Select Gene:",
                                          choices = gene_choices, selected = gene_choices[1]),
                              sliderInput("pt.size", "Point Size:", min = 0.5, max = 5, value = 2, step = 0.5),
                              sliderInput("alpha", "Transparency:", min = 0.1, max = 1, value = 0.8, step = 0.1)
                            ),
                            mainPanel(
                              plotOutput("spatialPlot", height = "600px")
                            )
                          )
                 ),
                 tabPanel("Cluster Plot",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("clusterSelect", "Select Cluster to Highlight:",
                                          choices = c("All", cluster_choices), selected = "All")
                            ),
                            mainPanel(
                              plotOutput("clusterPlot", height = "600px")
                            )
                          )
                 ),
                 tabPanel("Gene Expression Heatmap",
                          sidebarLayout(
                            sidebarPanel(
                              selectizeInput("heatGenes", "Select Genes (max 20):",
                                             choices = gene_choices, selected = gene_choices[1:5],
                                             multiple = TRUE,
                                             options = list(maxItems = 20))
                            ),
                            mainPanel(
                              plotOutput("heatmapPlot", height = "600px")
                            )
                          )
                 ),
                 tabPanel("Data Summary",
                          fluidPage(
                            h3("Spatial Seurat Object Summary"),
                            verbatimTextOutput("seuratSummary"),
                            h3("Spot-level Metadata"),
                            DTOutput("metaTable")
                          )
                 )
)

server <- function(input, output, session) {
  
  # Reactive expression to generate a spatial feature plot for a selected gene.
  output$spatialPlot <- renderPlot({
    req(input$gene)
    # Generate spatial plot using SpatialFeaturePlot
    p <- SpatialFeaturePlot(spatial_seurat, features = input$gene, pt.size.factor = input$pt.size, alpha = input$alpha) +
      ggtitle(paste("Expression of", input$gene)) +
      theme_minimal()
    print(p)
  })
  
  # Reactive expression to generate a cluster plot.
  output$clusterPlot <- renderPlot({
    # Default: show clusters as assigned by Seurat on the spatial image.
    # If a particular cluster is selected, highlight its spots.
    p <- SpatialDimPlot(spatial_seurat, label = TRUE, label.size = 3) +
      ggtitle("Spatial Distribution of Clusters") +
      theme_minimal()
    
    if (input$clusterSelect != "All") {
      # Create a new metadata column to highlight selected cluster.
      spatial_seurat$highlight <- ifelse(spatial_seurat$seurat_clusters == input$clusterSelect, "Selected", "Other")
      p <- SpatialFeaturePlot(spatial_seurat, features = "highlight", pt.size.factor = 1.6, alpha = input$alpha) +
        scale_color_manual(values = c("Selected" = "red", "Other" = "grey")) +
        ggtitle(paste("Highlighting Cluster", input$clusterSelect)) +
        theme_minimal()
    }
    print(p)
  })
  
  # Generate a heatmap of selected genes using Seurat's DoHeatmap.
  output$heatmapPlot <- renderPlot({
    req(input$heatGenes)
    # For heatmaps, we can use DoHeatmap. We need to subset our spatial data.
    # To speed up plotting, we use a random subset if there are many spots.
    cells_to_use <- colnames(spatial_seurat)
    if(length(cells_to_use) > 1000) {
      set.seed(123)
      cells_to_use <- sample(cells_to_use, 1000)
    }
    # Create the heatmap; scale data if needed.
    p <- DoHeatmap(spatial_seurat[, cells_to_use],
                   features = input$heatGenes,
                   size = 3) +
      ggtitle("Heatmap of Selected Genes") +
      theme_minimal()
    print(p)
  })
  
  # Output summary of the Seurat object
  output$seuratSummary <- renderPrint({
    spatial_seurat
  })
  
  # Output metadata table as an interactive datatable
  output$metaTable <- renderDT({
    meta_df <- spatial_seurat@meta.data
    datatable(meta_df, options = list(pageLength = 10, scrollX = TRUE))
  })
}

# ---- 5. Run the Shiny App ----
shinyApp(ui = ui, server = server)
