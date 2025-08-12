# Load the required packages

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ghibli)

# Standard Seurat Clustering Workflow 
cluster_seurat<-function(seurat){
  seurat <- NormalizeData(object = seurat)
  seurat <- FindVariableFeatures(object = seurat,selection.method = "vst", nfeatures = 2000)
  seurat <- ScaleData(object = seurat)
  seurat <- RunPCA(object = seurat)
  seurat <- FindNeighbors(object = seurat, dims = 1:40)
  seurat <- FindClusters(object = seurat, resolution = 1.2)
  seurat <- RunUMAP(object = seurat, dims = 1:40)
  plot<-DimPlot(object = seurat, reduction = "umap")
  result<-list(seurat=seurat,
               DimPlot=plot)
  return(result)
}



#Abundance Analysis on subtypes
Abudance_analysis<- function(data, column,disease,log_prop,covariates) {
  
  # Check if the specified column exists in the data frame
  if (!(column %in% names(data))) {
    stop("Specified column not found in the data frame.")
  }
  
  # Perform Wilcoxon test using compare_means
  wilcox_test <- wilcox.test(data[[column]] ~ data[[disease]])
  
  # Create missing column
  data$cell_missing <- as.numeric(is.na(data[[column]]))
  
  # Group by Cognitive Status and cell_missing
  AD_cell <- data %>%
    group_by(!!sym(disease), cell_missing) %>%
    summarise(count = n())
  
  # Get in correct format for Chi Squared tests
  chi_format <- pivot_wider(AD_cell, id_cols = cell_missing, values_from = count, names_from = !!sym(disease)) %>%
    dplyr::select(-1) # Exclude the first column
  
  # Perform Chi-squared test
  chisq_AD <- chisq.test(chi_format)
  
  #Perform Multivarite model
  if (log_prop==TRUE){
  multivar<-glm(as.formula(paste0("log(",column,")",covariates)),data=data)
  } else{
  multivar<-glm(as.formula(paste0(column,"~",covariates)),data=data)
  }
  # Output the plots and test results
  output <- list(
    wilcox_test = wilcox_test,
    chi_squared_test = chisq_AD,
    model=multivar
  )
  
  return(output)
}

get_abundance <- function(cluster_column_name, df, covariates,var_of_interest) {
  tryCatch({ 
  # Step 1: Mutate and create a new column that normalizes by total_cells
  new_column <- paste0(cluster_column_name, "_norm")
  df[[new_column]] <- df[[cluster_column_name]] / df$total_cells
  
  # Step 2: Create the formula for the model using the new normalized column and covariates
  formula_str <- paste0(new_column, " ~ ", covariates)
  
  # Step 3: Fit the generalized linear model
  multivar <- glm(as.formula(formula_str), data = df)
  
  # Step 4: Extract the p-value for the 3rd variable (Pathologic_diagnosis_of_ADyes)
  # Check if the variable exists in the model coefficients
  if (var_of_interest %in% rownames(summary(multivar)$coefficients)) {
    p_value <- summary(multivar)$coefficients["Pathologic_diagnosis_of_ADyes", "Pr(>|t|)"]
    coefficient <- summary(multivar)$coefficients["Pathologic_diagnosis_of_ADyes", "Estimate"]
    direction <- ifelse(coefficient > 0, "positive", "negative")
  } else {
    p_value <- NA  # Return NA if the variable is not present
    direction <- NA
  }
  
  # Return the p-value
  return(list(p_value = p_value, direction = direction))
  }, error = function(e) {
    # If an error occurs, return NA
    return(list(p_value = NA, direction = NA))
  })
}

# Main function to apply get_abundance to relevant columns
apply_abundance <- function(df,tag,covars,var_of_interest) {
  
  # Get the column names starting with the tag
  tag_cols <- grep(paste0("^",tag), names(df), value = TRUE)
  
  # Create a dataframe to store the results
  result_df <- data.frame(Cluster = character(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Apply the get_abundance function to each tag
  for (col in tag_cols) {
    # Get the p-value for each cluster column
    result <- get_abundance(col, df, covars,var_of_interest)
    
    # Append the result to the dataframe
    result_df <- rbind(result_df, data.frame(Cluster = col, P_Value = result$p_value, Direction = result$direction, stringsAsFactors = FALSE))
  }
  return(result_df)
}