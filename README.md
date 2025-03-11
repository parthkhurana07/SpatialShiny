# SpatialShiny: Interactive Spatial Transcriptomics Dashboard

SpatialShiny is an interactive dashboard built with R Shiny that processes and visualizes spatial transcriptomics data from 10x Genomics Visium platforms. By leveraging powerful R packages such as **Seurat**, **ggplot2**, **patchwork**, and **DT**, the dashboard offers dynamic visualizations—spatial feature plots, cluster overlays, and gene expression heatmaps—to facilitate intuitive exploration of tissue-level gene expression.

## Table of Contents
- [Introduction](#introduction)
- [Dataset](#dataset)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Running the App](#running-the-app)
  - [Exploring the Dashboard](#exploring-the-dashboard)
- [Project Structure](#project-structure)
- [License](#license)
- [References](#references)

## Introduction

Spatial transcriptomics maps gene expression data directly onto tissue sections, enabling the investigation of cellular organization and tissue architecture. SpatialShiny simplifies the process of visualizing and analyzing such data. The dashboard preprocesses spatial data using Seurat (applying SCTransform, PCA, clustering, and UMAP), and provides interactive modules to:
- Visualize gene expression patterns over the tissue image.
- Explore clustering results and highlight specific clusters.
- Generate heatmaps for selected genes.
- Inspect spot-level metadata and summaries.

## Dataset

The dataset used in this project is:
- **CytAssist_Fresh_Frozen_Sagittal_Mouse_Brain_filtered_feature_bc_matrix.h5**

This dataset, generated using the 10x Genomics Visium platform, contains spatial transcriptomics data from a fresh frozen, sagittal mouse brain sample. Ensure that you have the H5 file and its accompanying spatial directory (which includes images and spot coordinates).

## Requirements

- **R** (version 4.0 or above)
- **RStudio** (optional, but recommended)
- Required R packages:
  - `shiny`
  - `Seurat`
  - `ggplot2`
  - `patchwork`
  - `DT`
  - `dplyr`
- Data: Place the dataset (and associated spatial files) in an organized folder (e.g., `data/CytAssist_Mouse_Brain`).

## Installation

1. **Clone or download the project files.**
2. **Install required R packages:**

   ```r
   install.packages(c("shiny", "ggplot2", "patchwork", "DT", "dplyr"))
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("Seurat")
