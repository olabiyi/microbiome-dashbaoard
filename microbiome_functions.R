library(exactRankTests)
library(nlme)
library(pavian)
library(shiny)
library(DT)
library(tidyverse)
library(pheatmap)
#library(heatmaply)
library(here)
library(plotly)
library(ggpubr)
library(tools)
library(DescTools)
library(broom)
library(phyloseq)
library(vegan)
library(compositions)
library(car)
library(multcomp)
library(multcompView) # Package needed to display compact-letter display
#library(randomForest)
#library(plyr) # for the "arrange" function
#library(rfUtilities) # to test model significance
#library(caret) # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
library(gplots)
library(venn)
#library(d3heatmap) # for interactive heatmaps
#library(ancom.R)
library(patchwork)
library(factoextra)
library(NbClust)
library(pavian) # remotes::install_github("fbreitwieser/pavian")
library(qiime2R)
library(ALDEx2)
library(lme4) # Package for fitting linear and generalized linear mixed-effects models.
library(lmerTest)
library(emmeans) # Package to conduct mutliple comparisons on lm and lmem models
library(fmsb)
library(performance)
library(qiimer)
library(corrplot)
library(pairwiseAdonis) # #devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(furrr)
plan(strategy = multisession, workers=8)
library(mice) # replaceing missing values using predictive modelling


publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=ggplot2::margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=ggplot2::margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))


custom_palette <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                    "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                    "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033",
                    "#999933","#FF9933","#FFFAFAFF",colors()) 
# remove white colors
custom_palette <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow|azure|gray|#FFFAFAFF|aliceblue",
                                         x = custom_palette, ignore.case = TRUE))]


# including visualization
get_number_of_clusters <- function(mat, clustering_function=factoextra::hcut, ...){
  
  # A function to return the optimumn number of clusters the silhouette method
  # mat  is a matrix with the elements you are interested in determining the number of clusters a rows e.g genes or taxons
  # clustering_function - a clustering function to use. options are c(kmeans, cluster::pam, factoextra::hcut)
  # ... other arguments that can be passed to ?factoextra::fviz_nbclust
  
  # save and print cluster plot 
  invisible(print(cluster_plot <- fviz_nbclust(mat, hcut, method = 'silhouette', ...)))
  
  
  optimum_number_of_clusters <- cluster_plot$data[which.max(cluster_plot$data$y), 'clusters']
  
  return(optimum_number_of_clusters)
}

#  A more robust way to determine the optimum number of clusters using the NbClust package
# Without viusualization

number_of_clusters <- function(mat,...){
  
  # mat  is a matrix with the elements you are interested in determining the number of clusters a rows e.g genes or taxons
  # ... any arguement that can be passed to ?NbClust::NbClust
  # Example:
  # number_of_clusters(mat2plot, diss=NULL, distance = "euclidean", 
  #                min.nc=2, max.nc=10, 
  #               method = "complete", index = "silhouette")
  
  res <- NbClust::NbClust(mat,...)
  
  res$Best.nc['Number_clusters']
}


extract_clusters_from_pheatmap <- function(heatmap_result, mat){
  # heatmap_result is the output generate from running pheatmap
  # mat is the matrix used when plotting the heatmap
  
  number_of_clusters <- get_number_of_clusters(mat)
  mat.clust <- cbind(cluster = cutree(heatmap_result$tree_row, 
                                      k = number_of_clusters), mat)
  
  return(mat.clust)
}


# A set of functions to process the tables of different methods use for metagenomics

methods_stats <- function(taxon_level,unique_assignments,
                          analyses, taxon_tables){
  
  res_matrix <- map_dfr(.x = unique_assignments[[taxon_level]],
                        .f = function(taxa_name){
                          # Return a logical vector for whether a taxa_name was 
                          # detected in the different analyses methods
                          taxa_name <- taxa_name %>% 
                            # this is to ensure that we match the work precisely
                            str_replace_all('(.+)', '^\\1$') %>%  
                            # Escape meta characters
                            str_replace_all('\\.', '\\\\.')  %>% 
                            str_replace_all('\\[', '\\\\[') %>% 
                            str_replace_all('\\]', '\\\\]') %>% 
                            str_replace_all('\\(', '\\\\(') %>% 
                            str_replace_all('\\)', '\\\\)') %>% 
                            str_replace_all('\\-', '\\\\-') %>% 
                            str_replace_all('\\+', '\\\\+') %>% 
                            str_replace_all('\\*', '\\\\*')
                          
                          map_lgl(.x = analyses, 
                                  .f = function(analysis){ 
                                    
                                    # Get all the taxon assignments for the current analysis
                                    taxa_names <- taxon_tables[[analysis]][[taxon_level]] %>% 
                                      rownames() %>% trimws()
                                    # test if taxa_name was found in the taxa_names for the current analysis
                                    grepl(pattern = taxa_name, x = taxa_names) %>% 
                                      any()
                                    
                                    
                                    
                                  })
                          
                        }) %>% 
    mutate(!!taxon_level:=unique_assignments[[taxon_level]]) %>% 
    dplyr::select(!!taxon_level, everything())  %>% 
    arrange(!!taxon_level)
  
  
  # Select only relevant columns  
  res_matrix <- res_matrix  %>% 
    dplyr::select(1,(colSums(res_matrix[,2:ncol(res_matrix)]) > 0) %>% 
             which() + 1) %>% 
    mutate(percent=(rowMeans(across(where(is.logical))) * 100) %>% round,
           prevalence=glue::glue("{rowSums(across(where(is.logical)))}/{ncol(.)-1}"))
  
  
  return(res_matrix)
  
}



get_high_confidence_matrix <- function(taxon_level, analysis, taxon_tables,
                                       sub_comparison, REPLACE_NA = FALSE){
  
  taxon_table <- taxon_tables[[analysis]][[taxon_level]]
  high_confidence_taxa <- sub_comparison[[taxon_level]]
  taxa_pattern <- high_confidence_taxa %>% 
    # Make sure you match the exact word
    str_replace_all('(.+)', '^\\1$') %>% 
    str_c(collapse = "|") %>%
    # Escape meta characters
    str_replace_all('\\.', '\\\\.')  %>% 
    str_replace_all('\\[', '\\\\[') %>% 
    str_replace_all('\\]', '\\\\]') %>% 
    str_replace_all('\\(', '\\\\(') %>% 
    str_replace_all('\\)', '\\\\)') %>% 
    str_replace_all('\\-', '\\\\-') %>% 
    str_replace_all('\\+', '\\\\+') %>% 
    str_replace_all('\\*', '\\\\*')
  all_taxa <- rownames(taxon_table) %>% trimws()
  high_confidence_taxa_indices <- grep(pattern = taxa_pattern,
                                       x = all_taxa)
  # Return an empty matrix if the taxon_level does not exit for the analysis
  if(is.null(taxon_table)){
    samples <- taxon_table %>% colnames()
    ideal_matrix <- matrix(data = NA, 
                           nrow = length(high_confidence_taxa),
                           ncol = length(samples))
    
    rownames(ideal_matrix) <- high_confidence_taxa
    colnames(ideal_matrix) <- samples
    return(ideal_matrix)
  }
  
  # Return a matrix of the taxa counts for the high confidence taxa
  if(REPLACE_NA){
    # Retrieve the reads count matrix of the high confidence taxonomy assignments
    # then replaces NA rows (rows for which the taxa was not assigned for the analysis)
    # with zeros
    ideal_matrix <- taxon_table[high_confidence_taxa_indices,, drop=FALSE] %>% 
      mutate_if(is.numeric, .f=replace_na, replace=0)
    rownames(ideal_matrix) <-  all_taxa[high_confidence_taxa_indices]
  }else{
    
    # # Retrieve the reads count matrix of the high confidence taxonomy assignments
    # then drop NA rows (rows for which the taxa was not assigned for the analysis)
    ideal_matrix <- taxon_table[high_confidence_taxa_indices,, drop=FALSE] %>% 
      as.data.frame() %>% 
      drop_na() %>% 
      as.matrix()
    
  }
  
  return(ideal_matrix)
}


combine_high_confidence_matrices <- function(taxon_level, analyses,
                                             high_confidence_matrices){
  
  base_analysis <- analyses[1] 
  
  # initialize the df to be joined
  df <- high_confidence_matrices[[base_analysis]][[taxon_level]] %>% 
    as.data.frame(check.names=FALSE) %>% 
    rownames_to_column("taxa")
  
  for(analysis in analyses){
    
    if(analysis != base_analysis){
      
      # Join the high confidences metrics 
      df <- df   %>% 
        full_join(high_confidence_matrices[[analysis]][[taxon_level]] %>% 
                    as.data.frame(check.names=FALSE) %>% 
                    rownames_to_column("taxa"), by='taxa')
      
    }
    
  }
  rownames(df) <- df[,"taxa"] 
  taxa_column <- which(colnames(df) == "taxa")
  df <- df[,-taxa_column]
  return(df)
  
}


##########################################################
# Functions to annotate a plot (box plot) by adding a vertxivcal line
# to split groups and top section titles

get_metadata_value <- function(search, metadata, group_searched, group_returned){
  
  # group1 is 
  # USAGE:
  # get_metadata_value('E.G_CAR', metadata, "Group", "Biofilter_media")
  # get only the index of the first value just incase there a duplicate
  # values in the metadata for search vales in group searched
  first_index = which(metadata[,group_searched] == search)[1]
  
  return(as.character(metadata[,group_returned])[first_index])
}


mutate_df <- function(data, subet_from, subset_value,
                      search_group,
                      group1= "Biofilter_media",
                      group1_mapping=mapping_expr,
                      group2="Biofilter_configuration",
                      group3="Event"){
  # A function to add additional columns to a dataframe
  # in order to map specific values to a "search_group"
  # column of the dataframe data
  df <- data %>% 
    filter(!!sym(subet_from) == subset_value) %>% 
    # convert the search group variable to character
    dplyr::mutate(!!sym(search_group) := as.character(!!sym(search_group))) %>% 
    dplyr::mutate(
      !!sym(group1) := map_chr(.x = !!sym(search_group),
                               .f =  get_metadata_value,
                               metadata=metadata, 
                               group_searched=search_group,
                               group_returned=group1)
    ) %>% 
    dplyr::mutate(
      !!sym(group1) := case_when(!!!group1_mapping)
    ) %>%
    dplyr::mutate(
      !!sym(group2):= map_chr(.x = !!sym(search_group),
                              .f =  get_metadata_value,
                              metadata=metadata, 
                              group_searched=search_group,
                              group_returned=group2)
    )  %>% 
    dplyr::mutate(
      !!sym(group3) := map_chr(.x = !!sym(search_group),
                               .f =  get_metadata_value,
                               metadata=metadata, 
                               group_searched=search_group,
                               group_returned=group3)
    ) 
  
  
  return(df)
}

annotate_box_plot <- function(ylabel, df, x="Group",
                              y="value", y_interval=NULL, xlabel=NULL,
                              color="Biofilter_configuration",
                              color_vect= annotation_colors[["Biofilter_configuration"]],
                              text_lab1="GAC", text_lab2="SAND",
                              text_size=14, custom_palette=publication_format,
                              xlabel_order=NULL, set_x_order=TRUE,
                              adjust_x_angle=FALSE, x_text_angle=90,
                              add_title=FALSE, draw_vertical_line=TRUE){
  
  # Calculate points for adding vertical and horizontal lines
  mid_point <- length(xlabel_order)/2
  vertical_line <- mid_point + 0.5
  x_text_lab1 <- (mid_point/2) + 0.5
  x_text_lab2 <- (mid_point/2) + mid_point + 0.5
  max_y <- ifelse(test=round(max(df[[y]])) == 0,
                  yes = round(max(df[[y]]), 4),
                  no= max(df[[y]])
                  )
  min_y <- ifelse(test=floor(min(df[[y]])) == 0,
                  yes=round(min(df[[y]]), 4),
                  no=min(df[[y]])
  )
  
  
  
  
  # Create box plot and add colors
  gg <- ggplot(data = df, mapping = aes_string(x=x, y = y,color=color)) +
    geom_boxplot() + 
    geom_point() +
    scale_color_manual(values = color_vect)
  
  if(!is.null(y_interval)){
    # get the position for the topmost horizontal line
    to_add <- y_interval/2
    #to_minus <- y_interval/4
    #y_min <- min_y - to_minus
    y_title <- (max_y + to_add) + (to_add/2)
    y_title_lab <- (max_y + to_add) + (to_add)
    
    # set the y breaks interval
    if(add_title){
      gg <- gg + scale_y_continuous(limits = c(min_y, y_title_lab), 
                                    breaks=seq(min_y, max_y, by = y_interval))
    }else{
      
      gg <- gg + scale_y_continuous(limits = c(min_y,(max_y + to_add)), 
                                    breaks=seq(min_y, max_y, by = y_interval))
    }
    
    
    # draw a vertical line switch geom_vline otherwise draw with geom_segment
    if(draw_vertical_line){
      gg <- gg + geom_vline(xintercept = vertical_line)
    }else{
      
      gg <- gg +  annotate(geom = "segment", x = vertical_line,
                           y=min_y, xend = vertical_line, 
                           yend = y_title) 
    }
    
    gg <- gg +
      # Draw horizontal line
      geom_hline(yintercept = max_y + y_interval/4) +
      # custome theme
      custom_palette +
      # Annotate plot
      # Left-seide
      annotate(geom = "text", x = x_text_lab1, y= max_y + to_add,
               label=text_lab1, size=text_size) +
      # right-side
      annotate(geom = "text", x = x_text_lab2, y=max_y + to_add,
               label=text_lab2, size=text_size)
  }else{
    # Add vertical line
    gg <- gg +
      geom_vline(xintercept = vertical_line) +
      geom_hline(yintercept = max_y+0.2) +
      custom_palette + 
      # Annotate plot
      annotate(geom = "text", x = x_text_lab1, y=max_y+0.4,
               label=text_lab1, size=text_size) +
      annotate(geom = "text", x = x_text_lab2, y=max_y+0.4,
               label=text_lab2, size=text_size)
  }
  
  
  # Set x-text angle
  if(adjust_x_angle){
    gg <- gg + theme(axis.text.x = element_text(angle = x_text_angle,
                                                hjust = 0.5, vjust = 0.5 ))
  }
  
  # Set x-text to the x labels 
  if(set_x_order && (!is.null(xlabel_order))){
    gg <- gg + scale_x_discrete(labels=xlabel_order)
  }
  
  gg <- gg + labs(x=xlabel,y=ylabel%>% str_replace("\\.", " "),
                  color= color %>% str_replace("_", " ") %>% toTitleCase) 
  
  if(add_title){ 
    
    
    #gg <- gg + ggtitle(ylabel %>% toTitleCase)
    gg <- gg + 
      geom_hline(yintercept= y_title) +
      annotate(geom = "text", x = mid_point,
               y=y_title_lab,
               label=ylabel %>% toTitleCase, 
               size=text_size)
  }
  
  return(gg)
  
}

# Boxplot function specific for on the biofilter experiment
ancom_boxplot <- function(gg,microbes,intervals, Group_order, ncols=3){
  
  
  gg_plots <- map2(.x = microbes, .y = intervals,
                   .f = function(x,y){
                     
                     # Add mapping columns
                     df <- mutate_df(data=gg$data, subet_from="OTU_name",
                                     subset_value = x, search_group="Group")
                     # Enforce the grouping order
                     df$Group <- factor(x = df$Group, levels = Group_order)
                     # Make and annotate boxplot  
                     annotate_box_plot(y = "OTU", ylabel = sprintf("%s\n(Log abundance)", x), df=df,
                                       xlabel_order = xlabel_order,
                                       y_interval = y,
                                       text_size = 8)
                   })
  p <- wrap_plots(gg_plots, ncol = ncols, guides = 'collect') 
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size=20,face = "bold"))
  
  return(p)
}


##########################################################



#######################################################################

# A Function to process the taxonomy derived from protein coding sequences 
# produced by MetaErg
get_taxa_names <- function(taxon_level, taxonomy_vector){
  
  # taxon_level STR - one of  c('Kingdom', "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # taxonomy_vector - A vector of taxonomy assigments to be parsed
  # For example:  
  
  patterns <- c('Kingdom'='d',"Phylum"='p', "Class"='c', 
                "Order"='o', "Family"='f', "Genus"='g',
                "Species"="s")
  
  taxa_names <- strcapture(pattern = paste0(patterns[taxon_level],"__([-A-Za-z0-9 _]+);?"),
                           x=taxonomy_vector,
                           proto = list(taxon_level = character()))
  
  df <- data.frame(taxa_level=taxa_names)
  colnames(df) <- taxon_level
  
  return(df)
  
  
}


process_kraken <- function(file_path, delimiter="\t", 
                           drop_cols=c("taxRank", "taxID", "Max", "lineage"), kingdom=NULL, sample_order=NULL){
  
  # file_path {String] - path to parsed taxonlevel file derived after running pavian
  # delimiter [String] - the delimiter sepating sample counts in the pavian parsed count table
  # drop_cols [String vector] - Containing the name of columns to be dropped from the mutated dataframe
  # kingdom [String] - A regex pattern to return on assignments belong to one or more specific kingdoms e.g
  #                   "Bacteria" or "Virus|Archaea" to return only bacteria or a combination of only archae and 
  #                   viruses, respectively.
  #
  # df <- read.table(file=file_path, header=TRUE, 
  #                  sep="\t", check.names=FALSE, 
  #                  stringsAsFactors=FALSE )
  # 
  # taxonomy_table <- df[,"Lineage", drop=FALSE] %>% 
  #   mutate(Lineage=str_replace_all(Lineage,"&nbsp;","")) %>% 
  #   mutate(Lineage=str_replace_all(Lineage,"cellularorganisms>","")) %>%
  #   separate(col=Lineage, 
  #              into=c("Kingdom", "Super_Phylum", "Phylum", "Class", "Order", "Family", "Genus", 
  #                     sep=">"))
  
  # Read-in the count table
  df <- read_delim(file_path, delim = delimiter)
  
  non_microbial_indices <- grep(pattern = "Metazoa|Chordata|Nematoda|Arthropoda|Annelida|Brachiopoda|Mollusca|Cnidaria|Streptophyta|Primate|Chloroplast|Mitochondria",
                                x = df$lineage)
  
  if(!is_empty(non_microbial_indices)){
    
    df <- df[-non_microbial_indices,, drop=FALSE]
  }
  
  if(!is.null(kingdom)){
    
    # look in every row to find the indices of the kindom of interest
    kingdom_indices <- which(apply(df,MARGIN = 1,
                                   FUN = function(row) any(grepl(KINGDOM,row))
                                   )
                             )
    
    if(!is_empty(kingdom_indices)){
    
      df <- df[kingdom_indices,]
      
        }
    
  }
  
  df <- df  %>%
    dplyr::select(-!!drop_cols ) %>% 
    map_df(replace_na,0) %>% as.data.frame(check.names=FALSE)
  
  df <- aggregate(.~name, data=df, FUN=sum)
  rownames(df) <- df[,1]
  df <- df[,-1]
  
  clean_col_names <- gsub(pattern = "\\.cladeReads",replacement = "", x = colnames(df))
  colnames(df) <- clean_col_names
  
  if(!is.null(sample_order)){
    
    df <- df[,sample_order]
    
  }
  return(df)
  
}

process_mOtu_table <- function(file_path){
  
  # df <- read.table(file=file_path, 
  #                  header=TRUE, sep="\t", 
  #                  row.names=1, 
  #                  check.names=FALSE, 
  #                  stringsAsFactors=FALSE,
  #                  comment.char = "", 
  #                  skip=2)
  
  df <- read_delim(file=file_path, delim="\t",skip=2) %>% as.data.frame(check.names=FALSE)
  colnames(df)[1] <- "taxonomy"
  if(!grepl(pattern = "mOTU", x = file_path)){
    df$taxonomy <- gsub(pattern = ".+\\[(.+)?\\]", replacement = "\\1", x = df$taxonomy)
  }
  df <- aggregate(. ~ taxonomy, data = df, FUN = sum)
  rownames(df) <- df[,1]
  minus_one_row <- grep(pattern = "-1", x = rownames(df))
  df <- df[-minus_one_row,]
  df <- df[,-1]
  # Drop taxa without reads assigned to them
  df[rowSums(df) > 0, ]
}


# Generate count table for detecting proteins of interest
generate_count_table <- function(file_name, COLNAMES, gene_names) {
  
  if (file.size(file_name) > 0){
    df <- read_delim(file = file_name, delim = "\t", col_names = COLNAMES)
    df %>% count(stitle) %>% mutate(short_name=case_when(
      stitle == gene_names[1] ~ names(gene_names)[1],
      stitle == gene_names[2] ~ names(gene_names)[2],
      stitle == gene_names[3] ~ names(gene_names)[3],
      stitle == gene_names[4] ~ names(gene_names)[4],
      stitle == gene_names[5] ~ names(gene_names)[5],
      stitle == gene_names[6] ~ names(gene_names)[6],
      stitle == gene_names[7] ~ names(gene_names)[7],
      stitle == gene_names[8] ~ names(gene_names)[8],
      stitle == gene_names[9] ~ names(gene_names)[9],
      stitle == gene_names[10] ~ names(gene_names)[10],
      stitle == gene_names[11] ~ names(gene_names)[11],
      stitle == gene_names[12] ~ names(gene_names)[12],
      stitle == gene_names[13] ~ names(gene_names)[13],
      stitle == gene_names[14] ~ names(gene_names)[14],
      stitle == gene_names[15] ~ names(gene_names)[15],
      stitle == gene_names[16] ~ names(gene_names)[16],
      TRUE ~ stitle)
    ) %>% dplyr::select(short_name, protein=stitle, n)
  }
}

make_gene_count_matrix <- function(count_tables,gene_names,samples){
  # create a matrix of gene counts for all the sample
  genes_table <- matrix(data = 0, nrow = length(gene_names), ncol = length(samples))
  rownames(genes_table) <- names(gene_names)
  colnames(genes_table) <- samples
  
  for(sample in samples){
    
    for(gene in names(gene_names)){
      row <- which(count_tables[[sample]]$short_name == gene)
      if(!is_empty(row)){
        genes_table[gene,sample] <- count_tables[[sample]]$n[row]
      }
    }
  }
  return(genes_table)
}


map_sample_to_category <- function(group_col, sample_vector, metadata, sample_col="Sample.ID"){
  
  # Example: map_sample_to_category("Biofilter_media", rarefaction_table$sample, no_control_metadata)

  mapping <- map(.x = sample_vector, 
                 function(name) as.character(metadata[[group_col]])[metadata[[sample_col]] == name]) %>%
    unlist
  
  df <- data.frame(group=mapping)
  colnames(df) <- group_col
  
  return(df)
}


combine_mean_se_and_letters <- function(mean_se_df, compact_letters_df, colname, sep=" ± "){
  
  # mean_se_df - a dataframe with two columns to be pasted together using the sep provided
  #      The rownames of this dataframe must correspond to the rownames of compact_letters_df
  # compact_letters_df - a one column dataframe with a column of cattters to be added to the newly
  #             combined dataframe
  # sep  - a string specifying the seprator to be used to paste together the two columns of
  #     of mean_se_df
  # Colname - the column name to be used for the newly generated combined dataframe
  
  mean_se_df <- mean_se_df[rownames(compact_letters_df),]
  pasted_columns <- paste(mean_se_df[,1], mean_se_df[,2], sep = sep)
  column_and_letters <- paste0(pasted_columns,compact_letters_df[,1])
  combined_df <- data.frame(combined=column_and_letters)
  colnames(combined_df) <- colname
  rownames(combined_df) <- rownames(mean_se_df)
  return(combined_df)
  
}


combine_mean_and_se <- function(mean_se_df, colname, sep=" ± "){
  
  # mean_se_df - a dataframe with two columns to be pasted together using the sep provided
  # sep  - a string specifying the seprator to be used to paste together the two columns of
  #     of mean_se_df
  # Colname - the column name to be used for the newly generated combined dataframe
  
  pasted_columns <- paste(mean_se_df[,1], mean_se_df[,2], sep = sep)
  combined_df <- data.frame(combined=pasted_columns)
  colnames(combined_df) <- colname
  rownames(combined_df) <- rownames(mean_se_df)
  return(combined_df)
  
}


# A funcion to convert a named vector to a dataframe sorted in alphabetical order 
vector2dataframe <- function(named_vector, column_name){
  # named_vector - A named vector of compacter letters with the groups as the name and letters as the value
  
  sorted_group_names <- sort(names(named_vector))
  df <- data.frame(group=named_vector[sorted_group_names])
  colnames(df) <- column_name 
  rownames(df) <- sorted_group_names
  return(df)
  
}

# Get a list of species detected in a sample/treatment used to make venn diagrams
treatment_taxon_list <-function(treatment, taxon_table,taxa_are_rows=FALSE){
  # Function to get the taxon/species that were detected in a specified treatment from a taxon table
  #USAGE:treatments <- rownames(collapsed_tables$taxon_table)
  #treatment_taxon_list(treatments[2], collapsed_tables$taxon_table)
  if(taxa_are_rows){ 
    taxon_table <- t(taxon_table)
  }
  colnames(taxon_table)[taxon_table[treatment,] > 0]
  
}



create_distance_object <- function(ps, method = "PCoA", distance = "bray"){
  # A function to create an object similar to that produce by read_qza
  # so that it can be used with my pcoa_combination function
  # ps is a phloseq object with sample metadata and featuxre table
  
  # Perform pcoa
  pcoa <- ordinate(ps,method = method, distance = distance)
  
  # Intialize a 
  bray_res2 <- vector(mode = "list", length = 3)
  
  names(bray_res2) <- c("Eigvals", "ProportionExplained", "Vectors")
  
  # Eigen values
  bray_res2$Eigvals <- pcoa$values[,"Eigenvalues"]
  names(bray_res2$Eigvals) <- paste0("PC", 1:length(bray_res2$Eigvals))
  
  # Proportion of variance explained
  Cumul_eig_column <- ncol(pcoa$values) - 1
  pc2Tolast <- c(pcoa$values[,Cumul_eig_column][-1],0) - pcoa$values[,Cumul_eig_column]
  pc2Tolast <- pc2Tolast[-length(pc2Tolast)]
  bray_res2$ProportionExplained <- c(pcoa$values[,Cumul_eig_column][1], pc2Tolast)
  names(bray_res2$ProportionExplained) <- paste0("PC", 1:length(bray_res2$ProportionExplained)) 
  
  # Euhin vectors
  bray_res2$Vectors <- pcoa$vectors
  colnames(bray_res2$Vectors) <- paste0("PC", 1:ncol(bray_res2$Vectors))
  bray_res2$Vectors <- bray_res2$Vectors %>% 
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% dplyr::select("SampleID", everything()) 
  
  
  return(bray_res2)
}



combine_plots <- function(graphs,Title,Xlab,Ylab=NULL,Legend='right', ncols=NULL){
  
  # Get none null graphs and remove the legends and x labs
  graphs <- map(graphs[lengths(graphs)>0],
                .f = function(x){ x + theme(legend.position = Legend , axis.title.x = element_blank())} )
  
  if(!is.null(ncols)){
    
    p <- wrap_plots(graphs, ncol = ncols,
                    guides = 'collect')
  }else{
    
    p <- wrap_plots(graphs, ncol = if_else(condition = length(graphs) < 3,
                                           true = 1,false = 2),
                    guides = 'collect')
  }
  if(!is.null(Ylab)){
    
    p <- (wrap_elements(grid::textGrob(Ylab, rot = 90, gp=grid::gpar(fontsize=16,fontface='bold'))) + p) + 
      plot_layout(widths = c(1,30))
    
  }
  
  # combine plots
  p <- wrap_elements(grid::textGrob(Title, x=unit(0.53,"npc"), gp=grid::gpar(fontsize=18,fontface='bold'))) / 
    p / 
    wrap_elements(grid::textGrob(Xlab, x=unit(0.53,"npc"), gp=grid::gpar(fontsize=16,fontface='bold'))) + 
    plot_layout(heights = c(1,30,1))
  
  
  
  p
  
}

combine_enrichement_dfs <- function(enrich_dfs,GO_df){
  
  # Get only dataframes with enrichment results    
  dfs <- Filter(function(x) nrow(x) > 0, enrich_dfs)
  
  # Add comparison and Go category columns and re-arrange the columns
  sig_enrich_df <- map2_dfr(.x = dfs ,
                            .y = names(dfs), 
                            .f = function(x,y){
                              
                              if(!is_empty(x) && nrow(x) > 0){
                                x[,'GO_Category'] <- GO_df[GO_df$go_id %in% x$ID,"ontology"]
                                x[,'Comparison'] <- y
                                x %>% dplyr::select(Comparison,ID,GO_Category,everything())
                              }
                              
                              
                            })
  
  return(sig_enrich_df )
}

create_dt <- function(table2show, caption=NULL){
  DT::datatable(table2show,
                rownames = FALSE, # remove row numbers
                filter = "top", # add filter on top of columns
                extensions = "Buttons", # add download buttons
                options = list(
                  autoWidth = TRUE,
                  dom = "Blfrtip", # location of the download buttons
                  buttons = c("copy", "csv", "excel", "pdf", "print"), # download buttons
                  pageLength = 5, # show first 5 entries, default is 10
                  order = list(0, "asc") # order the title column by ascending order
                ),
                escape = FALSE # make URLs clickable) 
  )
  
}


# Function to perform single or multiway anova with assumptions testing and posthoc test with TukeyHSD
factorial_anova <- function(formula,Data,y,is_one_predictor=TRUE,letters=TRUE,...){
  
  # 1. formula is a formula
  # 2. Data is a dataframe containing the variables to be analysed
  # 3. y is a string specifying the response variable
  # 4. is_one_predictor a boolean specifying if you have only one predictor variable or not
  # 4. ... any parameter that can be passed to aov or lm functions
  
  # #############    Usage ##################
  # # For one predictor variable
  # result<- factorial_anova(formula = reformulate(termlabels = "growth_medium*sterility", response = "Acinetobacter"),
  #                          Data = new_data, y = "Acinetobacter", is_one_predictor = TRUE)
  # 
  # # For multiple predictor variables
  # result<- factorial_anova(formula = reformulate(termlabels = "growth_medium+sterility+growth_medium*sterility",
  #                                                response = "Acinetobacter"),
  #                          Data = new_data, y = "Acinetobacter", is_one_predictor = FALSE)
  # # Get the results
  # result$anova_test
  # result$assumptions_test
  
  
  
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    exp_aov <-aov(formula=formula,data=Data,...)
    tidy_aov <- broom::tidy(exp_aov)
    exp_tukey <- TukeyHSD(exp_aov)
    tidy_tukey <- exp_tukey %>% broom::tidy()
    if(letters){
    compact_letters <- multcompView::multcompLetters4(exp_aov,exp_tukey)
    }
    # Assumption testing for normality
    ifelse(test = length(unique(residuals(exp_aov))) > (length(residuals(exp_aov))/2),
           yes = shapiro_pval <- (shapiro.test(residuals(exp_aov)))$p.value, no = shapiro_pval<-0)
    # Assumption testing for outliers - NA occurs when p > 1 meaning there is no indication of outliers
    # in the dataset
    exp_outlier <- car::outlierTest(exp_aov)$bonf.p
    
  }else{
    shapiro_pval <- NA
    Levene_pval <- NA
  }
  
  # Levene_pval<- DescTools::LeveneTest(exp_aov)$'Pr(>F)'[1]
  # assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
  #                                Levene_pval = Levene_pval )
  # 
  
    Levene_pval<- DescTools::LeveneTest(lm(formula = formula, data = Data))$'Pr(>F)'[1]
    assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
                                   Levene_pval = Levene_pval,
                                   bonf_outlier_pval=exp_outlier)

  if(letters){
    results <- list(exp_aov,tidy_aov,assumptions_test,tidy_tukey,compact_letters)
    names(results) <- c("fit","anova_test","assumptions_test", "Tukey_posHoc_test", "compact_letters")
    }else{
  results <- list(exp_aov,tidy_aov,assumptions_test,tidy_tukey)
  names(results) <- c("fit","anova_test","assumptions_test", "Tukey_posHoc_test")
  
  }
  return(results)
}



ordinal_regression <- function(Formula, comparison, Data){
  
  library(ordinal)
  library(car)
  library(RVAideMemoire)
  library(rcompanion)
  library(broom)
  library(emmeans)
  library(multcomp)
  
  # Your dependent variable (y) must be a factor
  # 1. Formula is a formula
  # 2.comparison is a formula for comparing means see emeans::emeans() e.g Time|Location for Time by Location
  # 3. Data is a dataframe containing the variables to be analysed
  #comparison <- "Instructor|Question|extra"
  
  # Example:
  # res <- ordinal_regression(Formula = Likert.f ~ Instructor + Question + Instructor:Question, comparison = "Instructor|Question", Data = Data) 
  
  # Run ordinal regression
  model = clm(Formula,
              data = Data,
              threshold="symmetric")
  # Get the ANOVA table
  anova_table <- Anova(model, type = "II") %>% tidy
  
  # Get p-value for the model and pseudo R-squared
  model_test <- nagelkerke(fit = model)
  
  R_squared <- model_test$Pseudo.R.squared.for.model.vs.null[3,"Pseudo.R.squared"]
  p_value <-  model_test$Likelihood.ratio.test[1,"p.value"]
  
  model_significance <- data.frame(R_squared = R_squared, p_value=p_value)
  
  
  model.emm <- emmeans(model, reformulate(comparison))
  
  compact_letters <- multcomp::cld(model.emm, alpha = 0.05, 
                                   Letters = letters, reversed = TRUE) %>% 
    as.data.frame()
  
  # Check model assumptions
  nominal_effect_test <- nominal_test(model) %>%
    tidy() %>% 
    rename_all(~str_c("nominal_", .x)) %>%  
    dplyr::rename(term=nominal_term)
  
  scale_effect_test <- scale_test(model) %>% 
    tidy() %>% 
    rename_all(~str_c("scale_", .x)) %>%  
    dplyr::rename(term=scale_term)
  
  assumptions_test <- nominal_effect_test %>% inner_join(scale_effect_test,"term")
  
  
  results <- list(model,anova_table,assumptions_test,model_significance,compact_letters)
  names(results) <- c("fit","anova_test","assumptions_test", "model_significance", "compact_letters")
  
  
  return(results)
  
  
}




# Function to perform mixed model anova (aka within group anova) with assumptions testing and posthoc test with emmeans
mixed_model_anova <- function(formula,Data,y,comparison,...){
  library(lmerTest)
  library(emmeans)
  # 1. formula is a formula
  # 2. Data is a dataframe containing the variables to be analysed
  # 3. y is a string specifying the response variable
  # 4. comparison is a formula for comparing means see emeans::emeans() e.g Time|Location for Time by Location
  # 4. ... any parameter that can be passed to lmer functions
  
  # #############    Usage ##################
  # # random_factor = "Site"
  # result<- mixed_model_anova(formula =reformulate(c(str_c(c("Location", "Time"),collapse = "*"), glue::glue("(1|{random_factor})")), "value"),
  #                          Data = diversity_long.df %>% dplyr::filter(measures == "Observed"), 
  #                           y = "value", comaprison= "Time|Location")
  # 
  # # Get the results
  # result$anova_test
  # result$assumptions_test
  
  
  
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    lmer_fit <-lmer(formula=formula,data=Data,...)
    # performance::check_model(lmer_fit)
    # performance::check_homogeneity(lmer_fit)
    # performance::check_normality(lmer_fit)
    # performance::check_collinearity(lmer_fit)
    # performance::check_outliers(lmer_fit)
    lmer_aov <- anova(lmer_fit)
    tidy_aov <- broom::tidy(lmer_aov)
    summ_fit <- summary(lmer_fit)
    # GEt nmarginal means and standard errors
    lmer.emm <- emmeans(lmer_fit, reformulate(comparison))
    compact_letters <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE) %>% 
      as.data.frame()
    
    # Assumption testing for normality
    ifelse(test = length(unique(residuals(lmer_fit))) > (length(residuals(lmer_fit))/2),
           yes = shapiro_pval <- (shapiro.test(residuals(lmer_fit)))$p.value, no = shapiro_pval<-0)
    # Assumption testing for outliers - NA occurs when p > 1 meaning there is no indication of outliers
    # in the dataset
    exp_outlier <- car::outlierTest(lmer_fit)$bonf.p
    
  }else{
    shapiro_pval <- NA
    Levene_pval <- NA
  }
  
  
  
  #Levene_pval<- DescTools::LeveneTest(lm(formula = formula, data = Data))$'Pr(>F)'[1]
  factors <- str_split(comparison, "\\|") %>% unlist
  names(factors) <- factors
  Levene_pval <- map_dfc(.x = factors, .f = function(group){ 
    
    df <- data.frame(LeveneTest(reformulate(group, "residuals(lmer_fit)"), 
                                data = Data)$'Pr(>F)'[1])
    
    colnames(df) <- glue::glue("{group}_levene_pval")
    return(df)
  })
  
  assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
                                 bonf_outlier_pval=exp_outlier,
                                 Levene_pval)
  
  results <- list(lmer_fit,tidy_aov,assumptions_test,summ_fit,compact_letters)
  names(results) <- c("fit","anova_test","assumptions_test", "summary", "compact_letters")
  
  
  return(results)
}

# repeated measures anova
repeated_measures_anova <- function(formula,Data,y,comparison,...){
  
  # 1. formula is a formula
  # 2. Data is a dataframe containing the variables to be analysed
  # 3. y is a string specifying the response variable
  # 4. comparison is a formula for comparing means see emeans::emeans() e.g Time|Location for Time by Location
  # 4. ... any parameter that can be passed to lmer functions
  # Example
  # aov( value ~ Time*Location + Error(Site), data = diversity_long.df %>% dplyr::filter(measures == "Observed")) %>% tidy()
  
  
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    exp_aov <-aov(formula=formula,data=Data,...)
    tidy_aov <- broom::tidy(exp_aov)
    exp_tukey <- TukeyHSD(exp_aov)
    tidy_tukey <- exp_tukey %>% broom::tidy()
    if(letters){
      compact_letters <- multcompView::multcompLetters4(exp_aov,exp_tukey)
    }
    # Assumption testing for normality
    ifelse(test = length(unique(residuals(exp_aov))) > (length(residuals(exp_aov))/2),
           yes = shapiro_pval <- (shapiro.test(residuals(exp_aov)))$p.value, no = shapiro_pval<-0)
    # Assumption testing for outliers - NA occurs when p > 1 meaning there is no indication of outliers
    # in the dataset
    exp_outlier <- car::outlierTest(exp_aov)$bonf.p
    
  }else{
    shapiro_pval <- NA
    Levene_pval <- NA
  }
  
  # Levene_pval<- DescTools::LeveneTest(exp_aov)$'Pr(>F)'[1]
  # assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
  #                                Levene_pval = Levene_pval )
  # 
  
  Levene_pval<- DescTools::LeveneTest(lm(formula = formula, data = Data))$'Pr(>F)'[1]
  assumptions_test <- data.frame(shapiro_pval = shapiro_pval,
                                 Levene_pval = Levene_pval,
                                 bonf_outlier_pval=exp_outlier)
  
  if(letters){
    results <- list(exp_aov,tidy_aov,assumptions_test,tidy_tukey,compact_letters)
    names(results) <- c("fit","anova_test","assumptions_test", "Tukey_posHoc_test", "compact_letters")
  }else{
    results <- list(exp_aov,tidy_aov,assumptions_test,tidy_tukey)
    names(results) <- c("fit","anova_test","assumptions_test", "Tukey_posHoc_test")
    
  }
  return(results)
  
  
}



SE <- function(measure){sd(measure,  na.rm=TRUE)/sqrt(n())}

run_random_forest <- function(feature_table, metadata, group, remove_rare=TRUE, threshold=10/100,
                              transform_method="scale", transform_function=NULL,
                              ntree=501, top=10000, model_significance=TRUE, ...){
  # feature_table - A feature matrix with samples / observations as column names and features as row names
  # metadata - A metadata data frame mapping the samples to groups for classification / regression. The sample 
  #            names should be the row names and must match with the column names of the feature table.
  # group - A string specifying the name of the column in the metadata fro classification
  # remove_rare - A boolean set TRUE OR FALSE is rare features should be removed
  # threshold - A fraction specifying the threshold for removing rare features. e,g 10/100 to keep
  #            only features that are present in at least 10% of the samples
  # transform_method - A string specifying which transformation that should be applied
  #    on the samples before training the model
  # transform_function - A custom function to be applied on the samples before training the model.
  # ntree - the random of random trees tries to be used when running random forest
  # ... - other arguements that can be passed to the randomForest function. Please see  help(randomForest)
  
  # USAGE:
  # Remove rare features, normalize by scaling to Z-scores and then run random forest with 501 trees
  #  run_random_forest(taxon_tables$Strain, no_control_metadata, 'Biofilter_media', TRUE)
  
  # Set see for reproducibility of results
  set.seed(151) 
  
  taxon_table <- feature_table[,rownames(metadata)]
  
  if(remove_rare){
    # Remove rare feature - in this example, features not present in atleast 10% of our samples
    rare_removed_taxon_table <- remove_rare_features(feature_table = taxon_table,
                                                     cut_off_percent = threshold,
                                                     top=top)
  }else{
    
    rare_removed_taxon_table <- taxon_table
    
  }
  
  # convert to relative abundance
  rel_abund_taxon_table <- apply(X = rare_removed_taxon_table, MARGIN = 2, FUN = function(x) x/sum(x))
  
  # Transforming your data 
  if( transform_method == "none" ){
    taxon_table_norm <- rare_removed_taxon_table
  }else if(transform_method == "scale"){
    # Here we scale the data i.e convert the sample abundances to Z-scores 
    # (subtracting each sample's mean (center) and then dividing by the sample's standard deviation (scale))
    taxon_table_norm <- scale(rel_abund_taxon_table , center = TRUE, scale = TRUE)  
  } else if(transform_method == "asinh"){
    # Here we transform using the inverse hyperbolic sine and then mean center by sample
    taxon_table_norm<- scale(asinh(rel_abund_taxon_table), center=TRUE, scale=FALSE)  
  } else if(transform_method == "clr"){
    # centered log ratio transformation - clr 
    # Replace zero values with an estimate probability of zero
    d.n0 <- zCompositions::cmultRepl(rare_removed_taxon_table,  label=0, method="CZM", output="p-counts")
    # Apply clr transformation on every sample
    taxon_table_norm <- apply(X = d.n0  %>% as.matrix(), MARGIN = 2, FUN = compositions::clr)
  }else{
    
    if(!is.null(transform_function)){
      
      taxon_table_norm <- apply(X = rare_removed_taxon_table, MARGIN = 2, FUN = transform_function)
      
    }else{
      
      stop("You must provide a transformation mothod or transformation function to be applied to xyour samples")
    }
    
  }
  
  # Running the Model -  note that the metadata column will be the last column of each input dataframe
  taxon_table_norm_df <- data.frame(t(taxon_table_norm), check.names = FALSE) 
  
  # 2 parameters for a RF model are the number of trees in the forest (ntree) and the number of features randomly sampled at each node in a tree (mtry). The more trees you run in your forest the better the model will converge.
  #  choose odd numbers of trees to ensure there are never any ties for binary classification models.
  # use the default value of mtry
  taxon_table_norm_df[,group] <- metadata[rownames(taxon_table_norm_df), group]
  
  features <- taxon_table_norm_df[,1:(ncol(taxon_table_norm_df)-1)]
  y <- taxon_table_norm_df[ , ncol(taxon_table_norm_df)] 
  
  RF_fit <- randomForest::randomForest( x= features, y= y, ntree=ntree, importance=TRUE, proximities=TRUE, ...)
  #Assessing the model fit
  # RF_fit # Assessing accuracy by Out-of-bag error
  if(model_significance){
    # Permutation Test
    # Testing whether your model's performance metric is more extreme than expected by chance based on a permutation test is one method to assess model significance.
    # To run these significance tests with 1000 permutations:
    RF_fit_sig <- rfUtilities::rf.significance( x= RF_fit ,  xdata= features, nperm=1000 , ntree=ntree) 
  }
  # Accuracy Estimated by Cross-validation using the caret package
  # The first step is defining the parameters we want to use for training:
  fit_control <- caret::trainControl( method = "LOOCV" )
  # The command below will run the leave-one-out cross-validation, note the same ntree and mtry parameters are being used as above:
  RF_fit_loocv <- caret::train(features, y=y , method="rf", ntree=RF_fit$ntree,
                               tuneGrid=data.frame( mtry=RF_fit$mtry ) , trControl=fit_control )
  
  # Identifying Important Features
  RF_fit_imp <- as.data.frame( RF_fit$importance )
  RF_fit_imp$features <- rownames( RF_fit_imp )
  if(RF_fit$type == "classification"){
    RF_fit_imp_sorted <- plyr::arrange( RF_fit_imp  , MeanDecreaseAccuracy)
  }else{
    RF_fit_imp_sorted <- plyr::arrange( RF_fit_imp  , `%IncMSE`)  
  }
  
  RF_fit_imp_sorted$features <- factor(RF_fit_imp_sorted$features, levels = RF_fit_imp_sorted$features)
  twenty_most_important_features_indices <- (nrow(RF_fit_imp_sorted)-20):nrow(RF_fit_imp_sorted) # the last 20 rows
  
  if( RF_fit$type == "classification"){
    p <- ggplot(data = RF_fit_imp_sorted[twenty_most_important_features_indices,],
                mapping = aes(x=features,
                              y=MeanDecreaseAccuracy))
  }else{
    p <- ggplot(data = RF_fit_imp_sorted[twenty_most_important_features_indices,],
                mapping = aes(x=features,
                              y= `%IncMSE`))
  }
  
  p <- p +
    geom_col() + 
    coord_flip() + 
    labs(x=NULL, y=ifelse(test = RF_fit$type == "classification",
                          yes = "Mean Decrease in Accuracy (Variable Importance)",
                          no = "% Increase in Mean Squared Error (Variable Importance)")) +
    publication_format
  
  if(model_significance){
    
    final <- list(fit=RF_fit, 
                  significance=RF_fit_sig,
                  cross_validation=RF_fit_loocv,
                  feature_importance=RF_fit_imp_sorted,
                  plot=p) 
  }else{
    final <- list(fit=RF_fit, 
                  cross_validation=RF_fit_loocv,
                  feature_importance=RF_fit_imp_sorted,
                  plot=p) 
  }
  return(final)
}


make_venn_list <- function(...){
  # A function to generate a venn list for plotting from a taxon table that will be collapsed
  # by a group specified in a metadata
  # ... Arguements that can be passed to the function collapse_samples
  # Input taxon table should have sample names as rownames
  # USAGE: 
  # make_venn_list(taxon_table = t(taxon_tables$Genus), metadata = no_control_metadata, group = "Biofilter_media", convertToRelativeAbundance = TRUE)
  collapsed_tables <- collapse_samples(...)
  treatments <- rownames(collapsed_tables$taxon_table)
  names(treatments) <-  treatments
  VENN.LIST <- map(.x = treatments, 
                   .f =  treatment_taxon_list, taxon_table=collapsed_tables$taxon_table)
  
  # Get intesections
  a <- gplots::venn(VENN.LIST, show.plot=FALSE)
  
  # Store the intersections in a new object named inters
  inters <- attr(a,"intersections")
  
  final <- list(VENN.LIST=VENN.LIST, v.table=a, intersections=inters, treatments=treatments)
  return(final)
}


run_ANCOM_I <- function(feature_table, metadata, group, remove_rare=TRUE, threshold=10/100, ...){
  # feature_table - A feature matrix with samples / observations as column names and features as row names
  # metadata - A metadata data frame mapping the samples to groups for classification / regression. The sample 
  #            names should be the row names and must match with the column names of the feature table.
  # group - A string specifying the name of the column in the metadata fro classification
  # remove_rare - A boolean set TRUE OR FALSE is rare features should be removed
  # threshold - A fraction specifying the threshold for removing rare features. e,g 10/100 to keep
  #            only features that are present in at least 10% of the samples
  # ... - other arguments that can be passed to the ancom.R::ANCOM function. Please see  help(ancom.R::ANCOM)
  # USAGE: 
  # run_ANCOM_I(taxon_tables$Strain, no_control_metadata, 'Biofilter_media', TRUE, sig = 0.05, multcorr = 2, repeated = FALSE)
  
  if(remove_rare){
    # Remove rare features by on relative abundance
    feature_table <- apply(X = feature_table, MARGIN = 2, FUN = function(column) column/sum(column))
    abundant_features <- apply(X = feature_table, 
                               MARGIN = 1, 
                               FUN = function(x) max(x) >= threshold)
    feature_table <- feature_table[abundant_features, ,drop=FALSE]
    feature_table <- remove_rare_features(feature_table = feature_table,
                                          cut_off_percent = threshold,
                                          by_sample = FALSE)
  }
  # make the sample rows and then convert to dataframe
  feature_df <- feature_table %>% t %>% as.data.frame(check.names=FALSE)
  feature_df <- feature_df[rownames(metadata), ,drop=FALSE]
  
  # handling a case where their is only one feature to be analyzed
  if(ncol(feature_df) == 1) { 
    feature <- colnames(feature_df)
    dummy_feature <- glue::glue("{feature}2")
    feature_df[,dummy_feature] <- feature_df[,feature]
    }
  feature_df[, group] <- metadata[,group, drop=FALSE]
  ancom.out <- ancom.R::ANCOM( OTUdat = feature_df, ...)
  
  # ggplotly(gg_ancom)
  # cat('The significantly differntially abundant fetaures detected by ancom are: \n', ancom.out$detected)
  
  W_statistics_table <- data.frame(feature=colnames(ancom.out$dframe)[-ncol(ancom.out$dframe)],
                                   W=ancom.out$W) %>% 
    arrange(desc(W))
  
  if(length(ancom.out$detected) > 0){
    gg_ancom <- ancom.R::plot_ancom(ancom.out)
    final <- list(result=ancom.out, plot=gg_ancom, stat_table=W_statistics_table)
  }else{
    final <- list(result=ancom.out, stat_table=W_statistics_table)
  }
  
  
  return(final)
}

remove_rare_features <- function(feature_table, cut_off_percent=3/4, top=100, by_sample=TRUE){
  # feature_table - feature table matrix with samples as columns and features as rows
  #cut_off_percent - cut-off fraction  or decimal between 0.001 to 1 
  #of the total number of samples for determining the most abundant features
  # by default removes features that are not present in 3/4 of the total number of samples
  # top - an integer specifying the maximum number of abundant features to retian after filtering
  # by_sample - A boolean set to TRU or FALSE whether features should be removed based on
  #  there presences in less the then cut_off_perecent. If set to FALSE, features will be dropped
  # based on relative abundance at the set threshold /cut_off_percent
  # all also returns a maximum of top (by default 100) abundant features
  # remove rare OTUs/ASVs
  
  if(by_sample){
    # Filter by occurrence in a fraction of samples
    # Define a cut-off for determining what's rare
    cut_off <- cut_off_percent * ncol(feature_table)
    # Get the occurrence for each feature
    feature_occurence <- rowSums(feature_table > 0)
    # Get the names of the abundant features
    abund_features <- names(feature_occurence[feature_occurence >= cut_off])
  }else{
    # filter based on relative abundance
    abund_features <- apply(X = feature_table, MARGIN = 1, FUN = function(x) max(x) > cut_off_percent)
  }
  
  
  abun_features.m <- feature_table[abund_features,]
  # drop excess features
  if(nrow(abun_features.m) > top){
    abund_features <- names(abun_features.m %>% rowMeans() %>% sort(decreasing = TRUE))[1:top]
    # This is to ensure that my edge correlation function does not break
    # As it assumes that the ASV are in order
    abund_features <- DescTools::SortMixed(abund_features)
    abun_features.m <- feature_table[abund_features,]
  }
  return(abun_features.m)
}

process_metaphlan_taxonomy <- function(taxonomy, prefix='\\w__') {
  #function to process a metaphlan2 taxonopmy assigment table
  #1. ~ file_path is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
  
  #taxon_levels <- c("kingdom","phylum","class","order","family","genus","species", "strain")
  
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  
  #taxonomy[,'species'] <- paste(taxonomy[,'genus'],taxonomy[,'species'])
  # replace NAa with Other and delete the D_num__ prefix from the taxonomy names
  for (rank in colnames(taxonomy)) {
    #delete the taxonomy prefix
    taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],replacement = '')
    indices <- which(is.na(taxonomy[,rank]))
    taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
    #replace empty cell
    indices <- which(taxonomy[,rank] == "")
    taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
  }
  taxonomy <- apply(X = taxonomy,MARGIN = 2,
                    FUN =  gsub,pattern = "_",replacement = " ") %>% 
    as.data.frame(stringAsfactor=F)
  #taxonomy <- as.data.frame(taxonomy,stringsAsFactors=F)
  return(taxonomy)
}


#function for format a taxonomy assignment table by appending sufix to a known name
format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other", suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    #indices <- which(taxonomy[,taxa_index] == stringToReplace)
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}


fix_names<- function(taxonomy,stringToReplace,suffix){
  #1~ taxonomy is a taxonomy dataframe with taxonomy ranks as column names
  #2~ stringToReplace is a vector of regex strings specifying what to replace
  #3~ suffix is a string specifying the replacement value
  
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}


replace_duplicate_cells <- function(df, col1_index, col2_index, replacement="Other"){
  
  # Get the rows that are similar
  row_indices <- which(df[,col1_index] == df[,col2_index])
  
  # replace the values in the cells from the index of colum2 to the ncols of the dataframe 
  
  df[row_indices,col2_index:ncol(df)] <- replacement
  return(df)
}


combine_previous_and_current_cell <- function(df,current_col_index=8,suffixes, suffixes_indices){
  
  stopifnot(length(suffixes) == length(suffixes_indices), is.numeric(current_col_index))
  
  df[suffixes_indices,current_col_index] <-  paste0(df[suffixes_indices,current_col_index-1],suffixes)
  
  return(df)
  
}

# A function to generate taxon level count matrix based on a taxonomy table and an existing feature table
make_feature_table <- function(count_matrix,taxonomy,taxon_level, samples2keep=NULL){
  
  # EAMPLE:
  # make_feature_table(count_matrix = feature_counts_matrix,taxonomy = taxonomy_table,taxon_level = "Phylum")
  
  feature_counts_df <- data.frame(taxon_level=taxonomy[,taxon_level],
                                  count_matrix, check.names = FALSE, stringsAsFactors = FALSE)
  
  feature_counts_df <- aggregate(.~taxon_level,data = feature_counts_df,FUN = sum)
  rownames(feature_counts_df) <- feature_counts_df[,"taxon_level"]
  feature_table <- feature_counts_df[,-1]
  # Retain only taxa found in at least one sample
  taxa2keep <- rowSums(feature_table) > 0
  feature_table <- feature_table[taxa2keep,]
  
  if(!is.null(samples2keep)){
    feature_table <- feature_table[,samples2keep]
    # Retain only taxa found in at least one sample
    taxa2keep <- rowSums(feature_table) > 0
    feature_table <- feature_table[taxa2keep,]
  }
  
  return(feature_table)
}

get_rarefaction_df <- function(sample_index, rare_curve_list, sample=NULL) {
  # rare_curve_list - phloseq rarecurve list generated after runnig alpha_diversity function
  x <- attr(rare_curve_list[[sample_index]], which = "Subsample")
  sample_name <- names(x)[length(x)]
  if(is.null(sample_name) && !is.null(sample)){
    sample_name <- sample
  }
  data.frame(sample=sample_name, sequences=x, observed_species=rare_curve_list[[sample_index]])
  
}

get_last_assignment <- function(taxonomy_string, split_by =';', remove_prefix = NULL){
  # A function to get the last taxonomy assignment from a taxonomy string 
  split_names <- strsplit(x =  taxonomy_string , split = split_by) %>% 
    unlist()
  
  level_name <- split_names[[length(split_names)]]
  
  if(level_name == "_"){
    return(taxonomy_string)
  }
  
  if(!is.null(remove_prefix)){
    level_name <- gsub(pattern = remove_prefix, replacement = '', x = level_name)
  }
  
  return(level_name)
}

mutate_taxonomy <- function(df, taxonomy_column="taxonomy"){
  
  # make sure that the taxonomy column is always named taxonomy
  col_index <- which(colnames(df) == taxonomy_column)
  colnames(df)[col_index] <- 'taxonomy'
  df <- df %>% 
    pavian::zero_if_na() %>% 
    dplyr::mutate(taxonomy=map_chr(taxonomy,.f = function(taxon_name=.x){
      last_assignment <- get_last_assignment(taxon_name) 
      last_assignment  <- gsub(pattern = "\\[|\\]|'", replacement = '',x = last_assignment)
      trimws(last_assignment, which = "both")
    })) %>% 
    as.data.frame(check.names=FALSE, StringAsFactor=FASLE)
  # Ensure the taxonomy names are unique by aggregating duplicates
  df <- aggregate(.~taxonomy,data = df, FUN = sum)
  return(df)
}


process_kaiju_table <- function(file_path, long=TRUE, kingdom=NULL){
  
  kaija_table <- read_delim(file = file_path,
                            delim = "\t",
                            col_names = TRUE)
  # Remove non-microbial and unclassified assignments in this case Metazoa for animal assignments
  non_microbial_indices <- grep(pattern = "unclassified|assigned|Metazoa|Chordata|Nematoda|Arthropoda|Annelida|Brachiopoda|Mollusca|Cnidaria|Streptophyta",
                                x = kaija_table$taxon_path)
  
  if(!is_empty(non_microbial_indices)){
    kaija_table <- kaija_table[-non_microbial_indices,]
  }
  
  if(!is.null(kingdom)){
    kingdom_indices <- grep(pattern = kingdom ,
                                  x = kaija_table$taxon_path)
    if(!is_empty(kingdom_indices)){
    kaija_table <- kaija_table[kingdom_indices,]
    }
  }
  
  if(long){
    abs_abun_df <- pivot_wider(data = kaija_table %>% dplyr::select(sample,reads,taxonomy=taxon_path), 
                               names_from = "sample", values_from = "reads",
                               names_sort = TRUE) %>% mutate_taxonomy
    
    rel_abun_df <- pivot_wider(data = kaija_table %>% dplyr::select(sample,percent,taxonomy=taxon_path), 
                               names_from = "sample", values_from = "percent",
                               names_sort = TRUE) %>% mutate_taxonomy
    
    
    
  }else{
    
    abs_abun_df <-  kaija_table %>% dplyr::select(taxon_path, starts_with('read')) %>% 
      setNames(gsub('reads\\s+', '', names(.))) %>% 
      dplyr::rename(taxonomy=taxon_path) %>% 
      mutate_taxonomy
    
    rel_abun_df <- kaija_table %>% dplyr::select(taxon_path, starts_with('percent')) %>% 
      setNames(gsub('percent\\s+', '', names(.))) %>% 
      dplyr::rename(taxonomy=taxon_path) %>% 
      mutate_taxonomy
  }
  
  
  # Set the taxon names as row names, drop the taxonomy column and convert to a matrix
  rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
  rownames(rel_abun_df) <- rel_abun_df[,"taxonomy"]
  
  abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
  rel_abun_df <- rel_abun_df[,-(which(colnames(rel_abun_df) == "taxonomy"))]
  
  abs_abun_matrix <- as.matrix(abs_abun_df)
  rel_abun_matrix <- as.matrix(rel_abun_df)
  
  final_tables <- list("relative_table"=rel_abun_matrix,
                       "abundance_table"=abs_abun_matrix)
  return(final_tables)
  
}

process_phyloFlash_table <- function(file_path){
  
  
  phyloseq_table <- read_delim(file = file_path,
                               delim = "\t",
                               col_names = c('taxonomy', 'sample', 'reads'))
  
  abs_abun_df <- pivot_wider(data = phyloseq_table, 
                             names_from = "sample", values_from = "reads",
                             names_sort = TRUE) %>% 
    pavian::zero_if_na() %>% 
    as.data.frame(check.names=FALSE, StringAsFactor=FASLE)
  
  # convert to a matrix
  rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
  abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
  abs_abun_matrix <- as.matrix(abs_abun_df)
  
  return(abs_abun_matrix)
  
}

process_phyloFlash_table_species_table <- function(file_path, 
                                                   replacement="uncultured$|unidentified$|metagenome|uncultured bacterium|unclassified"){
  
  
  phyloflash_table <- read_delim(file = file_path,
                                 delim = "\t",
                                 col_names = c('taxonomy', 'sample', 'reads'))
  
  
  
  abs_abun_df <- pivot_wider(data = phyloflash_table, 
                             names_from = "sample", values_from = "reads",
                             names_sort = TRUE) %>% 
    pavian::zero_if_na() %>% 
    as.data.frame(check.names=FALSE, StringAsFactor=FASLE)
  # Remove plant and plastid assignments
  plant_indices <- grep(pattern = "Chloroplast;|Mitochondria;|Embryophyta|Metazoa|Chordata|Nematoda|Arthropoda|Annelida|Brachiopoda|Mollusca|Cnidaria|Streptophyta|Primate",x =  abs_abun_df$taxonomy, ignore.case = T)
  if(!is_empty(plant_indices)){
    # remove plant sequences
    abs_abun_df <- abs_abun_df[-plant_indices,]
  }
  abs_abun_df$taxonomy <- gsub(pattern = "\\[|\\]",replacement = '', x =  abs_abun_df$taxonomy, ignore.case = T)
  
  taxon_levels <- c("Kingdom", "Phylum", "Class",  "Order", "Family", "Genus", "Species")
  names(taxon_levels) <- taxon_levels
  
  taxonomy_table <- abs_abun_df %>% 
    dplyr::select(taxonomy) %>% 
    separate(taxonomy,into = taxon_levels, sep = ";")
  
  taxonomy_table <- process_metaphlan_taxonomy(taxonomy = taxonomy_table, prefix = "\\(.+\\)")
  
  
  stringToReplace <- c(grep(pattern = replacement,
                            x = taxonomy_table$Species, value = TRUE) %>% unique,
                       "Other")
  suffix <- rep(x = ";_", times=length(stringToReplace))
  
  taxonomy_table <- fix_names(taxonomy=taxonomy_table,stringToReplace=stringToReplace,suffix=suffix)
  
  
  feature_counts_matrix <- abs_abun_df[,-1] %>% as.matrix
  
  
  taxon_tables <- map(.x = taxon_levels,
                      .f = make_feature_table,
                      count_matrix = feature_counts_matrix,
                      taxonomy = taxonomy_table)
  
  
  return(taxon_tables)
}



# run permanova on a matrix using a choice distance metric, euclidean by default.
permanova <- function(mat, y_vect=NULL, dist_method="euclidean", formula=NULL, metadata=NULL){
  # EXAMPLE: permanova(vsd.df,metadata$treatment)
  
  if(!is.null(formula) && !is.null(metadata)){
    d <- vegan::vegdist(mat, method = dist_method)
    adonis2(formula,data = metadata)
  }else{
    d <- vegan::vegdist(mat, method = dist_method)
    adonis(d  ~ y_vect)
  }
  
}

# Make pcoa from distance object for an experiment with only one Treatment factor
make_pcoa <- function(distance_res, metadata,color_factor,
                      annotation_colors,
                      sample_column = "Sample.ID", 
                      custom_theme=publication_format,
                      axis1="PC1", axis2="PC2"){
  
  
  vectors_df <- distance_res$Vectors
  rownames(vectors_df) <- vectors_df[,"SampleID"]
  vectors_df <- vectors_df[rownames(metadata),
                           c("PC1", "PC2", "PC3")] %>% 
    bind_cols(metadata)
  
  All_p <- ggplot(vectors_df,aes_string(x=axis1, y=axis2, color=color_factor)) +
    geom_point(size=5) +
    scale_color_manual(values = annotation_colors[[color_factor]])
  
  All_p <- All_p + xlab(sprintf("%s: %1.2f%% of variance", axis1,
                                distance_res$ProportionExplained[axis1]*100)) + 
    ylab(sprintf("%s: %1.2f%% of variance", axis2,
                 distance_res$ProportionExplained[axis2]*100)) + 
    labs(color=color_factor) +
    coord_fixed() + custom_theme
  return(All_p)    
}

# axis1 <- "PC1"
# axis2 <- "PC2"

# Make PCoA plot combinations from QIIME2 output
pcoa_combination  <- function(distance_res, metadata,combinations,
                              annotation_colors,
                              sample_column = "Sample.ID", 
                              custom_theme=publication_format,
                              plot_dir=NULL, save_plots=TRUE, axis1="PC1", axis2="PC2",
                              shapes=NULL){
  # distance_res - A list generated after running read_qza(pcoa_results)
  
  vectors_df <- distance_res$Vectors
  rownames(vectors_df) <- vectors_df[,"SampleID"]
  vectors_df <- vectors_df[rownames(metadata), c("PC1", "PC2", "PC3")] %>%
    bind_cols(metadata)
  
  All_plots = list()
  for(column in 1:ncol(combinations)){
    
    color_factor <- combinations[1,column]
    shape_factor <- combinations[2,column]
    plot_title <- str_c(color_factor,' Vs ',shape_factor) %>% 
      str_replace_all('_',' ') %>% 
      toTitleCase()
    
    All_p <- ggplot(vectors_df,aes_string(x=axis1, y=axis2, color=color_factor, shape=shape_factor)) +
      geom_point(size=5) +
      scale_color_manual(values = annotation_colors[[color_factor]])
    
    if(!is.null(shapes)){
    All_p <- All_p +
      scale_shape_manual(values = shapes)
    }
    
    All_p <- All_p + xlab(sprintf("%s: %1.2f%% of variance", axis1, distance_res$ProportionExplained[axis1]*100)) + 
      ylab(sprintf("%s: %1.2f%% of variance", axis2, distance_res$ProportionExplained[axis2]*100)) + 
      labs(color=color_factor, shape=shape_factor) +
      coord_fixed() + custom_theme
    
    All_plots[[plot_title]] <- All_p
    
    if(save_plots && (!is.null(plot_dir))){
      
      ggsave(filename = paste0(plot_dir,'/', color_factor,'-',shape_factor,'-_pca.png'), 
             plot = All_p, device = 'png',
             width = 14, height = 8, dpi = 600 )
      
      
    }
    
    
  }
  
  return(All_plots)
}

# Plot 2D PCA
plot_pca <- function(mat, metadata,  axis1 = "PC1", axis2 = "PC2", color_factor,
                     shape_factor = NULL, sample_column='Sample.ID'){
  # mat - matrix with samples as rows and column as features
  
  common.ids <- intersect(rownames(metadata), rownames(mat))
  metadata <- metadata[common.ids,]
  samples <- rownames(metadata)
  mat <- mat[common.ids,]
  
  # samples <- rownames(metadata)
  # mat <- mat[samples,]
  
  mat_pca <- prcomp(mat)
  mat_pca_percentVar <- mat_pca$sdev^2/sum(mat_pca$sdev^2)
  names(mat_pca_percentVar) <- colnames(mat_pca$x)
  
  if(!is.null(shape_factor)){
    
    pca2plot <- data.frame(mat_pca$x[,c(axis1,axis2)],
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- ggplot(pca2plot,aes(x=get(axis1),
                             y=get(axis2),
                             col=Color,
                             shape=Shape)) +
      geom_point(size=5) +
      # geom_label_repel(aes(label=Sample)) +
      # geom_path(aes(color=Type,
      #               group=Type)) +
      xlab(sprintf("%s: %1.2f%% of variance", axis1, mat_pca_percentVar[axis1]*100)) + 
      ylab(sprintf("%s: %1.2f%% of variance", axis2, mat_pca_percentVar[axis2]*100)) + 
      labs(color=color_factor, shape=shape_factor) +
      coord_fixed()
    
  }else{
    
    pca2plot <- data.frame(mat_pca$x[,c(axis1,axis2)],
                           Color=metadata[,color_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    p <- ggplot(pca2plot,aes(x=get(axis1),
                             y=get(axis2),
                             col=Color)) +
      geom_point(size=5) +
      # geom_label_repel(aes(label=Sample)) +
      # geom_path(aes(color=Type,
      #               group=Type)) +
      xlab(sprintf("%s: %1.2f%% of variance", axis1, mat_pca_percentVar[axis1]*100)) + 
      ylab(sprintf("%s: %1.2f%% of variance", axis2, mat_pca_percentVar[axis2]*100)) + 
      labs(color=color_factor, shape=shape_factor) +
      coord_fixed()
    
  }
  
  return(p)
  
}


# Plot 3D PCA
plot3D_pca <- function(mat, metadata,  axis1 = "PC1", axis2 = "PC2", axis3 = "PC3", color_factor,
                     shape_factor = NULL, size_factor = NULL, sample_column='Sample.ID',
                     syms2use=NULL, colors2use=rainbow(4), sizes2use=NULL){
  # mat - matrix with samples as rows and column as features
  # Example
  # my_symbols <- c(NVNB='square',CAR='circle', NV='triangle-up')
  # 
  # # The size of the points mxust be specied as a ranged and not as explicit numbers
  # sizes <- c(20,60)
  # plot3D_pca(mat = t(clr_table), metadata = no_control_metadata,
  #            axis1 = "PC1", axis2 = "PC2", axis3 = "PC3",
  #            color_factor=color_factor, shape_factor = shape_factor,
  #            size_factor = size_factor, sample_column='Sample.ID',
  #            syms2use=my_symbols, 
  #            colors2use=annotation_colors[[color_factor]], 
  #            sizes2use=sizes)
  # 
  common.ids <- intersect(rownames(metadata), rownames(mat))
  metadata <- metadata[common.ids,]
  samples <- rownames(metadata)
  mat <- mat[common.ids,]
  
  mat_pca <- prcomp(mat)
  mat_pca_percentVar <- mat_pca$sdev^2/sum(mat_pca$sdev^2)
  names(mat_pca_percentVar) <- colnames(mat_pca$x)
  
  # Color, shape and size factors are provided
  if(!is.null(shape_factor) && !is.null(size_factor)){
    
    pca3plot <- data.frame(mat_pca$x[,c(axis1, axis2, axis3)],
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Size=as.numeric(metadata[,size_factor]),
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                  color = ~Color,
                  colors = colors2use,
                  text= ~Sample, 
                  symbol= ~Shape,
                  symbols = syms2use,
                 marker = list(sizemode = 'diameter'),
                  size= ~Size,
                  sizes= sizes2use) 
    
    # Color and shape factors are provided
  } else if(!is.null(shape_factor)){
    
    pca3plot <- data.frame(mat_pca$x[,c(axis1, axis2, axis3)],
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                 color = ~Color,
                 colors = colors2use,
                 text= ~Sample, 
                 symbol= ~Shape,
                 symbols = syms2use)
    
    # Only a Color factor is provided
  }else{
    
    pca3plot <- data.frame(mat_pca$x[,c(axis1, axis2, axis3)],
                           Color=metadata[,color_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))

    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z=~PC3,
                 mode = "markers",
                  colors = colors2use,
                  text=~Sample, 
                  color = ~Color)
    
  }
  
  p <-  p  %>%
    layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                     mat_pca_percentVar[axis1]*100)),
                        yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                     mat_pca_percentVar[axis2]*100)),
                        zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                     mat_pca_percentVar[axis3]*100))
                        
                        )
           
           )
  
  return(p)
  
}

plot_3D_pcoa <- function(distance_res, metadata, 
                         axis1 = "PC1", axis2 = "PC2", 
                         axis3 = "PC3", color_factor,
                         shape_factor = NULL,
                         size_factor = NULL, 
                         sample_column='Sample.ID',
                         syms2use=NULL, 
                         colors2use=rainbow(4),
                         sizes2use=NULL){
  # distance_res - A list generated after running read_qza(pcoa_results)
  # Example
  # my_symbols <- c(NVNB='square',CAR='circle', NV='triangle-up')
  # 
  # # The size of the points mxust be specied as a ranged and not as explicit numbers
  # sizes <- c(20,60)
  # plot3D_pca(mat = t(clr_table), metadata = no_control_metadata,
  #            axis1 = "PC1", axis2 = "PC2", axis3 = "PC3",
  #            color_factor=color_factor, shape_factor = shape_factor,
  #            size_factor = size_factor, sample_column='Sample.ID',
  #            syms2use=my_symbols, 
  #            colors2use=annotation_colors[[color_factor]], 
  #            sizes2use=sizes)
  
  
  
  
  vectors_df <- distance_res$Vectors
  rownames(vectors_df) <- vectors_df[,"SampleID"]
  vectors_df <- vectors_df[rownames(metadata), c(axis1,  axis2,  axis3)] %>%
    bind_cols(metadata)
  
  
  # Color, shape and size factors are provided
  if(!is.null(shape_factor) && !is.null(size_factor)){
    
    pca3plot <- data.frame(vectors_df,
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Size=as.numeric(metadata[,size_factor]),
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                 color = ~Color,
                 colors = colors2use,
                 text= ~Sample, 
                 symbol= ~Shape,
                 symbols = syms2use,
                 marker = list(sizemode = 'diameter'),
                 size= ~Size,
                 sizes= sizes2use) 
    
    # Color and shape factors are provided
  } else if(!is.null(shape_factor)){
    
    pca3plot <- data.frame(vectors_df,
                           Color=metadata[,color_factor],
                           Shape=metadata[,shape_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z = ~PC3, 
                 mode = "markers",
                 color = ~Color,
                 colors = colors2use,
                 text= ~Sample, 
                 symbol= ~Shape,
                 symbols = syms2use)
    
    # Only a Color factor is provided
  }else{
    
    pca3plot <- data.frame(vectors_df,
                           Color=metadata[,color_factor],
                           Sample=str_replace(string = metadata[,sample_column],
                                              pattern = "^X",
                                              replacement = ""))
    
    p <- plot_ly(pca3plot, x = ~PC1, y = ~PC2, z=~PC3,
                 mode = "markers",
                 colors = colors2use,
                 text=~Sample, 
                 color = ~Color)
    
  }
  p <-  p  %>%
    layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                     distance_res$ProportionExplained[axis1]*100)),
                        yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                     distance_res$ProportionExplained[axis2]*100)),
                        zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                     distance_res$ProportionExplained[axis3]*100))
                        
    )
    
    )
  
  return(p)
  
}



pca_combination <- function(mat, metadata,combinations,
                            annotation_colors,
                            sample_column = "Sample.ID", 
                            custom_theme=publication_format,
                            plot_dir=NULL, save_plots=TRUE, ...){ 
  # # mat - matrix with samples as rows and column as features
  # metadata - samples metadata
  # combinations - a matrix of factor combinations to be plotted
  # Example:
  # Plot without saving with a column in the metadata called "Sample.ID" containing the sample names
  # pca_combination(mat=t(clr_table), metadata = no_control_metadata,
  # combinations= combinations, annotation_colors=annotation_colors, save_plots=FALSE)
  # PLot and save plot
  # pca_combination(mat=t(clr_table), metadata = no_control_metadata,
  #                combinations= combinations, annotation_colors=annotation_colors,
  # plot_dir=plot_dir, save_plots=TRUE)
  
  All_plots = list()
  for(column in 1:ncol(combinations)){
    
    color_factor <- combinations[1,column]
    shape_factor <- combinations[2,column]
    plot_title <- str_c(color_factor,' Vs ',shape_factor) %>% 
      str_replace_all('_',' ') %>% 
      toTitleCase()
    
    All_p <- plot_pca(mat = mat, 
                      metadata = metadata, 
                      color_factor = color_factor, 
                      shape_factor = shape_factor,
                      sample_column = sample_column, ...) + 
      scale_color_manual(values = annotation_colors[[color_factor]]) + 
      custom_theme + 
      labs(color=toTitleCase(color_factor) %>% str_replace_all('_',' '),
           shape=toTitleCase(shape_factor) %>% str_replace_all('_',' '))
    
    
    All_plots[[plot_title]] <- All_p
    
    if(save_plots && (!is.null(plot_dir))){
      
      ggsave(filename = paste0(plot_dir,'/', color_factor,'-',shape_factor,'-_pca.png'), 
             plot = All_p, device = 'png',
             width = 14, height = 8, dpi = 600 )
      
      
    }
    
    
    
  }
  
  return(All_plots)
}




run_ANCOM_II <- function(feature_df, metadata, ...) {
  
  # feature_df - a dataframe or matrix with samples as columns and features as rows
  # metadata -  a dataframe with samples as rownames and variables to be tested
  
  feature_df <- feature_df[,rownames(no_control_metadata)] 
  
  res <- ANCOM(feature_df,metadata,...)
  
  return(res)
  
}

run_adlex_on_deseq <- function(ps, group, test, ...){
  
  # ps - a phyloseq object containing the feature table and a metadat containing the group to test 
  # group - a string specifying the name of the group to test that is within the phloseq metadata
  # test - the test to run - options are eithier 'kw' or  't' for multilevel or or two level group
  # comparisons respectively
  
  # EAMPLE:
  # result <- run_adlex_on_deseq(ps, "BioMedia_Event", test='kw') # Kruskal wallis test
  # resultt = 'kw <- run_adlex_on_deseq(ps, "Biofilter_media", test = 't') # t-test
  
  # Differnential abundance testing with ALDEx2
  # ALDEx2 (Anova-Like Differential Expression, 2) is a
  # bioconductor package that approaches the problem of differential
  # abundance of sequence feature counts as inherently compositional
  # data. In this approach, data are transformed to a Centered
  # Log-Ratio (CLR), while also making use of Monte Carlo resampling
  # from a Dirichlet Multinomial to address sparsity.
  
  #Load ALDEx2 package.
  library(ALDEx2)
  otu_table.m <- as.matrix(otu_table(ps))
  if(! taxa_are_rows(ps)) {otu_table.m <- t(otu_table.m)} 
  otu_table.df <- as.data.frame(otu_table.m)
  
  treatments <- as.character(get_variable(ps,group)) 
  stopifnot(length(treatments) == ncol(ps))
  
  # Generate instance of CLR
  # The mc.samples parameter determines the number of Monte Carlo
  # Dirichlet instances to generate. This is the sort of parameter
  # that should not affect the results unless it is too small (the
  # package authors claim 128 is “usually sufficient”), but larger
  # values incur a larger computational burden. When in doubt:
  # increase the value of mc.samples and check if the results
  # changed appreciably.
  
  # x <- aldex.clr (
  #   reads = otu_table.df,
  #   conds = treatments,
  #   mc.samples = 128,
  #   verbose = TRUE,
  #   useMC = TRUE
  # )
  
  
  # # Perform test
  # x.tt <- aldex.ttest(x, paired=TRUE)
  
  # Adlex test
  res_df <- aldex(reads = otu_table.df,conditions = treatments, test= test, ...)
  res_df <- apply(X = res_df, MARGIN = 2, FUN = round, digits=3) %>% 
    as.data.frame() %>%  
    rownames_to_column(var = "Taxa") 
  # sort the result table based the most significant test
  if(test == 'kw' ){
    res_df <- res_df %>% 
      arrange(glm.eBH)
  }else if(test == 't'){
    res_df <- res_df %>% 
      arrange(wi.eBH)
  }
  
  return(res_df)
  
}


#function for comparing between 2 groups - it performs both t-test and wilcox test
tTest_analysis <- function(x,y,Data){
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  
  #if all the values are the same
  #if(var(Data[,y]) == 0 || abs(max(Data[,y]) - min(Data[,y])) < 0.02){
  if(var(Data[,y]) == 0){
    indicator.var <- NA; shapiro_pval<- NA; tTest_pval <- 1; tTest_statistic <- NA;
    wilcox_pval <- 1; wilcox_statistic <- NA;
  }else if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    indicator.var <-var.test(Data[,y] ~ Data[,x])$p.value
    
    #check that you have unique values for y
    ifelse(test = length(unique(Data[,y])) > (length(Data[,y])/2),
           yes = shapiro_pval<- (shapiro.test(Data[,y]))$p.value, no = shapiro_pval<-0)
    
    tTest_pval <- t.test(Data[,y] ~ Data[,x])$p.value
    
    tTest_statistic <- t.test(Data[,y] ~ Data[,x])$statistic
    
    wilcox_pval<- wilcox.test(Data[,y] ~ Data[,x])$p.value
    
    wilcox_statistic<- wilcox.test(Data[,y] ~ Data[,x])$statistic
  }else{
    indicator.var <- NA; shapiro_pval<- NA; tTest_pval <- NA; tTest_statistic <- NA;
    wilcox_pval <- NA; wilcox_statistic <- NA;
  }
  results <- data.frame(shapiro_pval = shapiro_pval ,
                        var_pval = indicator.var, 
                        tTest_statistic = tTest_statistic, 
                        tTest_pval = tTest_pval,
                        wilcox_statistic = wilcox_statistic,
                        wilcox_pval = wilcox_pval)
  
  return(results)
}

#function for comparing more than 2 groups, it performs a one way anova and kruskal wallis test
anova_analysis <- function(x,y,Data){
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  if ((length(unique(Data[,y])) > (length(Data[,y])/20))){
    indicator.aov <-aov(Data[,y] ~ Data[,x])
    anova_pval <- anova(indicator.aov)$'Pr(>F)'[1]
    anova_fval <- anova(indicator.aov)$'F value'[1]
    group_DF <- anova(indicator.aov)$Df[1]
    residual_DF <- anova(indicator.aov)$Df[2]
    ifelse(test = length(unique(residuals(indicator.aov))) > (length(residuals(indicator.aov))/2),
           yes = shapiro_pval<- (shapiro.test(residuals(indicator.aov)))$p.value, no = shapiro_pval<-0)
    Levene_pval<- LeveneTest(lm(Data[,y] ~ Data[,x], data = Data))$'Pr(>F)'[1]
    
    kruskal_statistic <- kruskal.test(Data[,y] ~ Data[,x], data = Data)$statistic
    kruskal_pval<- kruskal.test(Data[,y] ~ Data[,x], data = Data)$p.value
  }else{
    anova_pval <- NA; anova_fval <- NA; group_DF <- NA; residual_DF <- NA; shapiro_pval <- NA;
    Levene_pval <- NA; kruskal_statistic <- NA; kruskal_pval<- NA;
  }
  results <- data.frame(group_DF = group_DF ,
                        residual_DF = residual_DF, 
                        anova_fval = anova_fval ,
                        anova_pval = anova_pval,
                        shapiro_pval = shapiro_pval,
                        Levene_pval = Levene_pval , 
                        kruskal_statistic = kruskal_statistic, 
                        kruskal_pval = kruskal_pval)
  return(results)
}

# A function for pairwise comparison between groups using 
# wilcoxcon, t.test, anova and kruslkal wallis
compare_groups <- function(x,y,Data, ...){
  library(ggpubr)
  #x is a string specifying the name of the independent variable
  #y is a string specifying the name of the dependent variable 
  #Data is a dataframe containing the variables to be analysed
  # .. any arguemnt that can be passed to ggpubr's compare_means function
  if((length(unique(Data[,y])) > (length(Data[,y])/10))){
    Formula <- paste(sprintf("%s ~ ", y), x, sep = "")
    Forma <- as.formula(Formula)
    comparison <- compare_means(Forma,data = Data, ...) 
    
  }else{
    comparison <- data.frame(.y. = c(NA,NA), group1=c(NA,NA), 
                             group2=c(NA,NA), p = c(NA,NA), p.adj=c(NA,NA),
                             p.format=c(NA,NA), p.signif=c(NA,NA), method=c(NA,NA))
    
  }
  
  return(comparison)
}

multiple_posHoc_test <- function(data, group, colname="Taxon", post_hoc_method ="dunn.test"){
  #funtion to perform posthoc tests on multiple variable in a dataframe
  
  # data - is a dataframe that contains the group and the numeric variables to test
  # group - a a string that specifies the name of the grouping variable
  # colname - is a string that specifies the general name for the numeric variables
  # post_hoc_method - is a string that specifies the posthoc method either "dunn.test" or "tukeyHSD"
  library(FSA)
  numeric_variables <- colnames(data)[-(which(colnames(data) == group))]
  all_comparisons <- data.frame()
  
  
  if(post_hoc_method == "dunn.test"){
    for (variable in numeric_variables) {
      res<- dunnTest(data[,variable],data[,group])
      current_comparison <-cbind(rep(x = variable, times=nrow(res$res)),res$res)
      all_comparisons <- rbind(all_comparisons,current_comparison)
    }
  }else if(post_hoc_method == "tukeyHSD"){
    
    for (variable in numeric_variables) {
      aov.model <- aov(formula = reformulate(termlabels = group,response = variable), data = data)
      res <- TukeyHSD(aov.model)
      current_comparison <- cbind(rep(x = variable, times=nrow(res[[group]])),res[[group]])
      colnames(current_comparison)[1] <- colname
      current_comparison <- as.data.frame(current_comparison)
      current_comparison$comparison <- rownames(current_comparison)
      current_comparison <-  current_comparison[,c(colname,"comparison","diff","lwr","upr","p adj")]
      
      all_comparisons <- rbind(all_comparisons,current_comparison)
    }
    
  }else{
    stop("post_hoc_method must be either dunn.test or tukeyHSD")
  }
  rownames(all_comparisons) <- 1:nrow(all_comparisons)
  colnames(all_comparisons)[1] <- colname
  return(all_comparisons)
}

replace_if_na <- function(df, replacement){
  
  
  stopifnot(!is.null(df))
  if (nrow(df) > 0) {
    df[is.na(df) | df < 0] <- replacement
  }
  return(df)
}

#  A function to perform pairwise comparison between levels in a factor
# using either 'Dunn', 'T-test', or 'Wilcoxon' tests.
# A little more work is needed
pairwise_test <- function(Data, group, y_col, metric, method="Dunn"){
  #  A function to perform pairwise comparison between levels in a factor
  # using either 'Dunn', 'T-test', or 'Wilcoxon' tests.
  # Impossible comparisons are return as 1 meaning no significant difference
  
  # Example:
  # Dunn Test
  # pairwise_test(Data=data2use, group="Group", y_col="value", metric="Observed", method="Dunn")
  # T-test
  # pairwise_test(Data=data2use, group="Group", y_col="value", metric="Observed", method="T-test")
  # Wilcoxon test
  # pairwise_test(Data=data2use, group="Group", y_col="value", metric="Observed", method="Wilcoxon")
  
  
  methods <- c("Dunn", "T-test", "Wilcoxon")
  if(!(method %in% methods)){
    stop("method can only be one of 'Dunn', 'T-test', or 'Wilcoxon' ")
  }
  
  if(method == "Dunn"){
    # Dunn test
    sig_table <- multiple_posHoc_test(data = Data[,c(group,y_col)], group = group,colname = y_col) %>% 
      separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>% 
      dplyr::mutate(Metric=metric, method=method)  %>% 
      dplyr::select(-!!sym(y_col)) %>% 
      dplyr::rename(p=P.unadj, p.adj=P.adj) %>% 
      dplyr::mutate(p.format=round(p,digits = 2))
    
    sig_table <- rstatix::add_significance(data = sig_table,p.col = 'p',output.col = 'p.signif')
    
    sig_table <- sig_table %>% 
      dplyr::select(Metric,group1, group2, Z, p, p.adj, p.format, p.signif, method)
    
  }else{
   
    # Initializations for t-test or wilcoxon test
    group_combination <- combn(droplevels(Data)[,group] %>% as.character, 2)
    
    sig_table <- t(group_combination) %>% as.data.frame(stringAsFactors=FALSE, check.names=FALSE)
    colnames(sig_table) <- c("group1", "group2")
    
    # Delete comparisons between the same factor level
    sig_table <- sig_table %>% dplyr::filter(group1 != group2)
    
    sig_table$p <- 8888
    
    if( method == "T-test"){  
      # T-test
      res <- pairwise.t.test(x =  Data[,y_col], g =  Data[,group], p.adjust.method = "none")$p.value %>%
        as.data.frame(stringAsFactors=FALSE, check.names=FALSE) %>% 
        replace_if_na(replacement = 1)
      
    }else{
      # Wilcoxon test
      res <- pairwise.wilcox.test(x =  Data[,y_col], g =  Data[,group], p.adjust.method = "none")$p.value %>%
        as.data.frame(stringAsFactors=FALSE, check.names=FALSE) %>% 
        replace_if_na(replacement = 1)
    }
    
    # Extract the p-vales for every comparison 
    for( row in (1:nrow(sig_table)) ){
      
      group1 <- sig_table[row,"group1"]
      group2 <- sig_table[row,"group2"]
      
      
      if( (!is.null(res[group1,group2]))){
        sig_table[row,'p'] <- res[group1,group2]
        
      }else{
        
        # Assign a very large number to be filtered out for non existent comparions
        
        sig_table[row,'p'] <- 9999
        
      }
      
    }
      # Remove non-existing comparisons
      sig_table <- sig_table %>% 
                    dplyr::filter(!is.na(p)) %>% 
                    dplyr::distinct()
    
      # Remove the null rows that were replaced with 9999 and comparison that did not
      # Exist having the vale 8888
      sig_table <- sig_table %>% dplyr::filter(!(p %in% c(9999,8888)))
      
      sig_table <- sig_table %>% 
        dplyr::mutate(Metric=metric, method=method)  %>% 
        dplyr::mutate(p.format=round(p,digits = 2),
                      p.adj=p.adjust(p,method = "holm"))
      
      sig_table <- rstatix::add_significance(data = sig_table,p.col = 'p',output.col = 'p.signif')
      sig_table <- sig_table %>% 
        dplyr::select(Metric,group1, group2, p, p.adj, p.format, p.signif, method)
      
    
    
  }
  
  return(sig_table)
}




#rarefy the otu table at the lowest sampling depth or at a specified depth
# In this example, the rarefaction depth chosen is the 90% of the minimum
#sample depth in the dataset without replacement

rarefaction <- function(ps,depth=0.9*min(sample_sums(ps)), with_replace=FALSE, ...) {
  #1. ~ ps is a phyloseq object
  #2. ~ depth is the refaction depth e.g 0.9*min(sample_sums(ps)), min(sample_sums(ps))
  #    a number of your choice e.g 19449
  #3 .. - any arguement that can be passed to vegan::rarecure
  #rarefy the otu table
  ps.rarefied <- rarefy_even_depth(physeq = ps, 
                                   sample.size = depth,
                                   rngseed = 1, 
                                   replace = with_replace, 
                                   verbose = FALSE)
  
  
  
  #plot rarefaction curve
  # calclate the a rough estimate of the step sample step size for plotting
  # This meant to allow the time for plotting to be relatively constant
  # regardless of the sample depth
  step = (50*depth)/1000
  rare_gg <- rarecurve(t(otu_table(ps.rarefied)), step=step, cex=0.5, ...)
  #print(rare_gg)
  result <- list(ps.rarefied, rare_gg)
  names(result) <- c("ps.rarefied","rare_curve")
  return(result)
}

# Alpha diversity
alpha_diversity <- function(ps.rarefied, metadata, group='treatments', 
                            metrics=c("Observed", "Shannon"),
                            group_order=NULL){
  
  #ps.rarefied - rarefied phyloseq object
  # metatdata - metadata table describing each sample
  # group - independent variable within you metadata for which the results should be grouped
  # group_order - a vector specifying the order within the group for plotting
  # metrics - a vector specifying the alpha diversity metrices to plot
  
  
  original_sample_names <- rownames(sample_data(ps.rarefied))
  #get all  sorts of alpha diversity measurements 
  diversity.df <- estimate_richness(ps.rarefied,
                                    measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  # Do this to make sure that the sample names of diversity.df and that of the metadata are in sync 
  rownames(metadata) <- make.names(rownames(metadata))
  #get the samples present in both the metadata and diversity table
  common.ids <- intersect(rownames(diversity.df), rownames(metadata))
  diversity.df <- diversity.df[common.ids,]
  sample_names <- rownames(diversity.df)
  metadata <- as.data.frame(metadata[common.ids,])
  rownames(metadata) <- sample_names
  colnames(metadata) <- group
  
  rownames(metadata) <- original_sample_names
  rownames(diversity.df) <- original_sample_names
  # REmove the prefix X from the sample names that may have been added when running the commands above
  # rownames(diversity.df) <- gsub(pattern = "^X",replacement = '', rownames(diversity.df))
  # rownames(metadata) <- gsub(pattern = "^X",replacement = '', rownames(metadata))
  diversity.df <- data.frame(metadata,diversity.df, check.names = FALSE)
  
  #mean and standard deviation of the diversity matrices
  diversity_meanSd <- diversity.df %>%  
    dplyr::select(-c(se.chao1,se.ACE))%>% 
    dplyr::group_by(!!sym(group))%>% 
    dplyr::summarise_all(list(mean=mean,std=sd))
  
  #number of samples per treatment group 
  number_of_samples <- diversity.df %>%  
    dplyr::select(-c(se.chao1,se.ACE)) %>% 
    dplyr::group_by(!!sym(group))  %>% 
    dplyr::summarise(N= n())
  
  
  #function to calculate standard error
  SE<- function(x){
    sd(x)/sqrt(length(x))
  }
  
  diversity_se <- diversity.df %>%  
    dplyr::select(-c(se.chao1,se.ACE)) %>% 
    dplyr::group_by(!!sym(group))  %>% 
    dplyr::summarise_all(list(se=SE))
  
  
  
  
  #re-order the levels of the grouping factor if need be
  if(is.null(group_order)){
    sample_data(ps.rarefied) <- metadata
  }else{
    metadata[,group] <- factor(x =metadata[,group], levels = group_order)
    # Add metadata to the phyloseq object
    sample_data(ps.rarefied) <- metadata
  }
  
  
  #plot alpha diversity
  richness_gg <- plot_richness(ps.rarefied,
                               x=group,
                               measures=metrics) + geom_boxplot()
  
  result <- list(diversity.df,metadata, number_of_samples, diversity_meanSd,diversity_se,richness_gg)
  names(result) <- c("diversity.df","metadata", "N", "mean_sd","SE","plot")
  return(result)
}

set_ylim_faceted_box_plot <- function(data=bacteria_alpha_hydro$plot$data,
                                      x='MembraneType',y='value',
                                      facet_variable='variable',
                                      sub_level_variables= c('Observed','Shannon'),
                                      ylimits=list(c(2800,3600),c(5.2,6)),ncol=2,nrow=1,
                                      xLABs=c('',''), YLABs=c('','')){
  #Parameters
  # data is a dataframe that contains the variables to be plotted
  # x is a string representing the x-axis variable
  # y is a string representing the y-axis variable
  # facet_variable is a string representing the variable in data that contains the 
  #  faceting variables
  # sub_level_variables is a sting vector of the faceting variables
  # ylimits is a list containing vectors of the ylimits for as many sub_level_variables
  # ncol is the number of columns for arranging the plots
  # nrow is the number of rows for arranging the plots
  # XLABs is a character vector that contains the labels for your x-axis
  # YLABs is a character vector that contains the labels for your y-axis
  
  library(ggplot2)
  
  publication_format <- theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.ticks.length=unit(-0.15, "cm"),
          axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
          axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
          axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
          axis.text = element_text(size = 16,face ='bold', color = 'black'),
          legend.position = 'right', legend.title = element_text(size = 15,face ='bold', color = 'black'),
          legend.text = element_text(size = 14,face ='bold', color = 'black'),
          strip.text =  element_text(size = 14,face ='bold', color = 'black'))
  
  # Create a list that will contain all the plots  
  graphs <- vector(mode = 'list',length = length(sub_level_variables))
  names(graphs) <- sub_level_variables
  
  
  # Make and set the limit for each box plot
  for(index in seq_along(sub_level_variables)){
    
    current_data <- data[data[,facet_variable] == sub_level_variables[index],]
    
    graphs[[sub_level_variables[index]]] <- ggplot(data = current_data,
                                                   aes_string(x=x, y=y)) + 
      geom_point() +
      geom_boxplot() + 
      facet_wrap(as.formula(paste0('~',facet_variable))) +
      coord_cartesian(ylim=ylimits[[index]]) + 
      labs(x=xLABs[index],y=YLABs[index]) + 
      publication_format
  }
  
  # Arrange the plots in grids
  gg <- gridExtra::grid.arrange(grobs=graphs,ncol=ncol, nrow=nrow)
  
  return(gg)
}

# Beta diversity
betadiversity <- function(ps.rarefied,
                          method = "PCoA", 
                          distance= "bray",
                          groups = "treatments",
                          metadata, ellipse=FALSE, 
                          color_var='MembraneType',
                          shape_var=NULL,
                          palette = c("red","green","blue","yellow","pink")){
  library(ggpubr)
  library(tools)
  library(phyloseq)
  
  # Add metadata to the phyloseq object
  #note a column in your metadata must be called samples else you'll get an error
  #here I add sample column which is the rownames of the passed metadata
  metadata <- data.frame(sample=rownames(metadata),metadata)
  sample_data(ps.rarefied) <- metadata
  common.ids <- intersect(rownames(metadata), colnames(otu_table(ps.rarefied)))
  metadata <- metadata[common.ids,]
  
  # make pcoa plot
  pcoa <- ordinate(ps.rarefied,method = method, distance = distance)
  
  
  pcoa_gg <- plot_ordination(physeq = ps.rarefied,
                             ordination = pcoa,
                             type = "samples",
                             axes = c(1,2), 
                             color = color_var,
                             shape = shape_var)
  
  if(length(groups)== 2){
    
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var, shape = shape_var, size = 6, 
                         ellipse = ellipse, palette = palette)
  }else{
    pcoa_gg <- ggscatter(data = pcoa_gg$data, x = 'Axis.1', y = 'Axis.2',
                         color = color_var,size = 6, ellipse = ellipse, palette = palette)  
  }
  
  
  #plot the axes ticks inwards 
  pcoa_gg <- pcoa_gg + theme(axis.ticks.length=unit(-0.15, "cm"),
                             axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                             axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
  
  
  if(method == "PCoA"){
    Cumul_eig <- ncol(pcoa$values) - 1
    pcoa_gg <- pcoa_gg +
      labs(x= sprintf('PC1 [%s%%]', round(pcoa$values[1, Cumul_eig] * 100, digits = 2)),
           y= sprintf('PC2 [%s%%]', round((pcoa$values[2, Cumul_eig] - pcoa$values[1, Cumul_eig]) * 100,
                                          digits = 2)))
  }
  
  pcoa_gg <- ggpar(pcoa_gg,font.x = c(18,'bold.italic','black'), font.tickslab =  c(16,'bold','black'),
                   font.y = c(18,'bold.italic','black'),
                   legend.title = toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups)), 
                   legend = 'right', font.legend = c(14, "bold", "black"))
  legend_titles <- toTitleCase(gsub(pattern = '\\.', replacement = ' ', x = groups))
  
  if(length(groups) == 2){
    pcoa_gg <- pcoa_gg + labs(color=legend_titles[1], shape=legend_titles[2])
  }
  
  #find if there is a significant clustering between treatments using adonis test aka 
  #PERMANOVA
  #calculate and get your distance matrix - here we use bray curtis but check
  # the help page for more options
  bray_dist <- phyloseq::distance(ps.rarefied, method="bray")
  permanova<- vector(mode = 'list', length = length(groups))
  names(permanova) <- groups
  
  for (group in groups){
    permanova[[group]] <- adonis(bray_dist  ~ metadata[,group])
  }
  #pcoa_gg
  result <- list(pcoa,pcoa_gg,permanova,bray_dist)
  names(result) <- c("ordination","plot","permanova","bray_dist")
  
  return(result)
  
}


draw_distance_boxplot <- function(distance_matrix=as.matrix(bacteria_beta_mem$bray_dist),
                                  metadata= bacteria_alpha_mem$metadata,
                                  group= "MembraneType",
                                  YLAB= 'Bray Curtis distance') {
  library(ggpubr)
  library(tidyverse)
  #Parameters
  # distance_matrix - is a distance matrix comparing samples with zeros diagonally
  # metatdata - is a dataframe with ONLY ONE factor column (categorical variable)
  # mapping samples to treatment names
  # group - is a string representing the name of the categorical variable as it appears
  # in the supplied metadata
  # YLAB - is a string to used as the label for the y-axis
  #publication ready theme
  publication_format <- theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.ticks.length=unit(-0.15, "cm"),
          axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
          axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
          axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
          axis.text = element_text(size = 16,face ='bold', color = 'black'),
          legend.position = 'right', 
          legend.text = element_text(size = 14,face ='bold', color = 'black'),
          strip.text =  element_text(size = 14,face ='bold', color = 'black'))
  
  # combine the metatdata with the distace matrix
  distance_df <- cbind(metadata,distance_matrix)
  
  ####### Get each within group distance matrix
  #within group distance is the distance between samples from the same group
  within_groups_df <- data.frame()
  groups <- levels(metadata[,group])
  
  for(sub_level in groups){
    group_samples <- rownames(distance_df)[distance_df[,group] == sub_level]
    group_dist_mat <- distance_df[group_samples,group_samples]
    # Drop the last column because it doesn't give any new information
    group_dist_mat <- group_dist_mat[,-(ncol(group_dist_mat))]
    
    # Get the distance vector for the current sub_level 
    group_dist_vector <- unlist(sapply(X = group_dist_mat, 
                                       FUN = function(x) x[(which(x == 0)+1):length(x)]))
    #row bind the previous group distance to the current group distance
    within_groups_df <- rbind(within_groups_df,
                              data.frame(groups=rep(x = sub_level,times=length(group_dist_vector)),
                                         distance=group_dist_vector))
  }
  
  
  ###### Get the distance between groups
  between_group_dist_df <- data.frame() #initialize an empty dataframe 
  
  #get all possible 2 way comparisons between each group
  group_combinations <- combn(groups,2)
  
  for(index in seq_along(group_combinations)){
    group1 <- group_combinations[1,index]
    group2 <- group_combinations[2,index]
    group1_samples <- rownames(distance_df)[distance_df[,group] == group1]
    group2_samples <- rownames(distance_df)[distance_df[,group] == group2]
    
    group_dist_mat <- distance_df[group1_samples,group2_samples]
    current_between_groups <- gather(group_dist_mat) %>% 
      mutate(groups=rep(x = paste(group1,group2, sep = " Vs "),times=nrow(.))) %>% 
      select(-key,distance=value) %>% 
      select(groups,distance)
    between_group_dist_df <- rbind(between_group_dist_df,current_between_groups)
  }
  
  
  
  #make within and between groups box plots  
  gg_within <- ggboxplot(data = within_groups_df, x = 'groups',
                         y = 'distance',
                         xlab =sprintf('Within %s',gsub(pattern ='\\.',replacement = " ",x = group)),
                         ylab = YLAB,
                         add = "mean_se") +  ylim(c(0.0,1.0)) + publication_format
  
  gg_between <- ggboxplot(data = between_group_dist_df, x = 'groups',
                          y = 'distance',
                          xlab = sprintf('Between %s',gsub(pattern ='\\.',replacement = " ",x = group)),
                          ylab = YLAB,
                          add = "mean_se") +  ylim(c(0.0,1.0)) + publication_format +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #check if there is only one possible comparison between the groups in order to adjust
  # the x-axis text 
  if(ncol(group_combinations) == 1){
    gg_between <- gg_between  + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  }
  
  plots <- list(gg_within,gg_between)
  names(plots) <- c("within","between")
  
  results <- list(within_groups_df,between_group_dist_df,plots)
  names(results) <- c("within_groups_dataframe","between_groups_dataframe", "plots")
  
  return(results)
}

# Complete this function
group_rare_longF <- function(data, grouping_col, value_col, threshold){
  
  abundant_features <- data %>% 
    group_by(!!grouping_col) %>% 
    summarise(max=max(!!value_col)) %>%
    arrange(desc(max)) %>% filter(max >= !!threshold)
  
  rare_taxa <- data %>% 
    group_by(!!grouping_col) %>% 
    summarise(max=max(!!value_col)) %>%
    arrange(desc(max)) %>% 
    filter(max < !!threshold) %>% 
    pull(!!grouping_col)
  
}

#function to group rare taxa or return a table with the rare taxa
group_low_abund_taxa <- function(abund_table,threshold=0.05, rare_taxa=FALSE) {
  # abund_table is a relative abundance matrix with taxa as columns and  samples as rows
  #rare_taxa is a boolean specifying if only rare taxa should be returned
  #If set to TRU then a table with only the rare taxa will be returned 
  #intialize an empty vector that will contain the indices for the
  #low abundance columns/ taxa to group
  taxa_to_group <- c()
  #intialize the index variable of species with low abundance (taxa/columns)
  index <- 1
  
  #loop over every column or taxa check to see if the max abundance is less than the set threshold
  #if true save the index in the taxa_to_group vector variable
  while(TRUE){
    
    for (column in seq_along(abund_table)){
      if(max(abund_table[,column]) < threshold ){
        taxa_to_group[index] <- column
        index = index + 1
      }
    }
    if(is.null(taxa_to_group)){
      threshold <- readline("please provide a higher threshold for grouping rare taxa, only numbers are allowed   ")
      threshold <- as.numeric(threshold)
    }else{
      break
    }
    
  }
  
  
  if(rare_taxa){
    abund_table <- abund_table[,taxa_to_group,drop=FALSE]
  }else{
    #remove the low abundant taxa or columns
    abundant_taxa <-abund_table[,-(taxa_to_group), drop=FALSE]
    #get the rare taxa
    # rare_taxa <-abund_table[,taxa_to_group]
    rare_taxa <- subset(x = abund_table, select = taxa_to_group)
    #get the proportion of each sample that makes up the rare taxa
    rare <- rowSums(rare_taxa)
    #bind the abundant taxa to the rae taxa
    abund_table <- cbind(abundant_taxa,rare)
    #rename the columns i.e the taxa
    colnames(abund_table) <- c(colnames(abundant_taxa),"Rare")
  }
  
  return(abund_table)
}

#function to generate an ordered list of microbes from an oTU table based on relative abundance
sort_microbes_on_relative_abundance <- function(otu_table, ORDER= 'decreasing') {
  
  #this function generates an ordered list of microbes from an oTU table based on relative abundance
  #inputs
  #1. otu table - with samples as columns and microbes as row names  
  #2. ORDER is a string that species the order of the microbes based on relative abundance
  #the options is 'decreasing' or 'increasing' for decreasing and increase list based on relative abundance respectively
  
  #ensure the OTU table is a dataframe
  otu_table<-as.data.frame(otu_table)
  #re-arranging the OTU table
  colNames<- colnames(otu_table) #get the colnames of the otu table
  otu_table$microbes <- rownames(otu_table) #create a new column that will contaisn the names of the microbes
  
  #re-order the columns of the otutable
  col.order <- c("microbes", colNames)
  otu_table <- otu_table[,col.order,drop=FALSE]
  
  #stack the otu table by microbes
  stacked_table <- Stack_data(Data = otu_table, 
                              lab_columns = 'microbes', 
                              num_of_numericVariables = length(2:ncol(otu_table)), 
                              col_num_of_numericVariables = 2:ncol(otu_table))
  
  
  colnames(stacked_table) <- c("microbes","samples","relative_abundance")
  
  #calculate the mean abundance for each microbe
  microbes_mean<-aggregate(.~microbes,data = stacked_table, FUN = mean)
  if(ORDER == 'decreasing'){
    #order the dataframe by microbe relative abundance from the highest to the lowest
    ordered_microbes <-microbes_mean[order(microbes_mean$relative_abundance, decreasing = TRUE),]
  }else if(ORDER == 'increasing'){
    ordered_microbes <-microbes_mean[order(microbes_mean$relative_abundance, decreasing = FALSE),]
  }else{
    #cat("You have not selected a valid order either ORDER = 'decreasing' or 'increasing' ")
    stop("You have not selected a valid order either ORDER = 'decreasing' or 'increasing' ")
  }
  #get the list of microbes with decreasing abundance
  ordered_microbes_list <- ordered_microbes[,'microbes']
  
  return(ordered_microbes_list)
}


make_abundance_plots <- function(otu_table, grouping_info,group, 
                                 taxon_level,Xlab, position = 'right',
                                 split_by = ';',plot_type='box',
                                 group_rare_taxa =TRUE, 
                                 threshold = 0.05,
                                 convert_percent=TRUE,
                                 sort_box_by_abundance=TRUE, 
                                 box_order = 'decreasing',
                                 box_legend_title=NULL,
                                 palette= c("blue","purple","red","orange","green",
                                            "brown","skyblue","gold","grey","black"),
                                 fun=sum, write_plot=FALSE) {
  #this function makes abundance plots from qiime absolute rarefied output taxon tables 
  #and can optionally group rare taxa at a defined threshhold for plotting
  
  ####inputs
  #1. otu_table - the absolute rarefied output taxon table file produced by running summarize_taxa.py with
  #the -a option on the rarefied OTU generated by QIIME after core diversity analysis
  #2. group-  is a string specifying the independent variable in the mapping file that is currently under consideration
  #3. taxon_level - a string specifying the label for the taxon level under consideration
  #4. Xlab - a string specifying the X-axis label for plotting
  #5. position - a string either c("top","bottom","left","right", "none") specifying the position of the legend
  #6. split_by - a string specifying the delimiter for spliting the taxa names either "_" or ";"
  #7. group_rare_taxa - a boolean either TRUE or FALSE if rare taxa at the specified threshold should be grouped
  #8. threshold - a float or integer between 0 and 1 specifying the threshold for grouping rare taxa
  #9.plot_type - the type of plot to be made the only option is 'bar' if not passed a box plot will be made instead 
  #10. sort by abundance = should your box plot be sorted by abundance?
  #11. box_order - how should the boxes in the boxplot be ordered ? c('decreasing or increasing')
  #12. box_legend_title - title to use for boxplot legend
  #13. palette <- set a differnt color palette for the box plot
  #14.fun - a function, without brackets, to apply in order to collapse the samples when making stacked bar plots
  #set the default Y-axis label
  Ylab <- 'Relative abundance'
  
  
  #split the names of the taxonomy assignment
  if (split_by == '_'){
    split_names <- strsplit(colnames(otu_table), split = '_')
  }else if(split_by == ';'){
    split_names <- strsplit(colnames(otu_table), split = ';')
  }else{
    split_names <- colnames(otu_table)
  }
  #intialize an empty vector that will contain the formated microbe name
  microbe <- c()
  #loop over the split_file_name list of lists getting the first element of each list
  #and setting the respective treatment value to the rightful index in the vector
  #Treatments
  for (index in 1:length(colnames(otu_table))) {
    k <- split_names[[index]]
    microbe[index] = split_names[[index]][length(k)]
  }
  
  #convert the microbe vector of names to a dataframe
  
  microbe <- as.data.frame(microbe)
  colnames(microbe) <- taxon_level #rename the column
  
  #rename the columns of the otu table as the formatted microbes name
  colnames(otu_table)<- microbe[,taxon_level]
  
  
  #convert the oTU table to a dataframe - this trick helps with getting unique microbe name
  #because it forces the column names to be unique
  #this trick allows for unique names of the columns when
  #we have duplicated names like Ambigous taxa, other and uncultured bacterium
  otu_table2 <- data.frame(otu_table)
  taxon_names <-colnames(otu_table2)
  ##################### debugging code########################
  # #print the duplicated values                            #
  # print(taxon_names[duplicated(taxon_names)])             #
  # #print the length of the duplicated values              #
  # print(length(taxon_names[duplicated(taxon_names)]))     #
  # # if(length(taxon_names[duplicated(taxon_names)] > 0)){ #
  # #   next                                                #
  ########################################################## # # }
  
  
  #get the dataframe of the independent variable under consideration
  independent_variable.df <- subset(x = grouping_info, select=group)
  
  common.ids <- intersect(rownames(independent_variable.df), rownames(otu_table2))
  independent_variable.df <- independent_variable.df[common.ids,,drop=FALSE]
  otu_table2 <- otu_table2[common.ids,,drop=FALSE]
  #column bind the independent variable column to the otu_table
  otu_table.df <- cbind(independent_variable.df,otu_table2)
  
  if(plot_type == 'bar'){
    #sum the values for each treatment group
    otu_table.df<- aggregate(reformulate(termlabels = group, response = '.'), 
                             data = otu_table.df, FUN = fun)
    #store the independent variable you working on in a dataframeme to be used later
    new_indar <- subset(otu_table.df, select = group)
    #rename the rownames as the variable for the new grouped independent variable
    rownames(otu_table.df) <- otu_table.df[,1]
    
    #drop that independent variable so that it can be converted to a matrix
    otu_table.df <- otu_table.df[,-1]
    
    #calculate the relative abundance of the taxa for each group
    otu_table_abund.df <- as.data.frame(t(apply(X = otu_table.df, MARGIN = 1, FUN = function(x) x/sum(x))))
    
    if(group_rare_taxa){
      #group low abundance taxa
      otu_table_abund.df <- group_low_abund_taxa(abund_table = otu_table_abund.df, threshold = threshold)
    }
    #add the independent variable to the abundance table
    otu_table_abund.df <- cbind(new_indar,otu_table_abund.df)
    #stack the data
    stacked_table <- Stack_data(Data =otu_table_abund.df, 
                                lab_columns = group, 
                                num_of_numericVariables = length(2:ncol(otu_table_abund.df)), 
                                col_num_of_numericVariables = 2:ncol(otu_table_abund.df))
    
    i <- which(colnames(stacked_table) == 'Pathogens')
    colnames(stacked_table)[i] <- taxon_level
    
    if(convert_percent){
      stacked_table$Relative.abundance <- stacked_table$Relative.abundance *100
      Ylab <- 'Relative abundance (%)'
    }
    #get my custum made palette
    custom_palette <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                        "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                        "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                        "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033",
                        "#999933","#FF9933","#FFFAFAFF",colors())
    #plotting bar plots
    gg <- ggbarplot(data = stacked_table, 
                    x = group, 
                    y = 'Relative.abundance', 
                    fill = taxon_level,
                    xlab = Xlab,
                    ylab = Ylab, 
                    palette = custom_palette)
    #plot the axes ticks inwards 
    gg <- gg +theme(axis.ticks.length=unit(-0.15, "cm"),
                    axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                    axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
    
    
    gg <- ggpar(gg,font.x = c(18,'bold.italic','black'), font.tickslab =  c(16,'bold','black'),
                font.y = c(18,'bold.italic','black'),legend.title = toTitleCase(taxon_level), 
                legend = position, font.legend = c(14, "bold", "black"))
    if(write_plot){
      ggsave(filename = (paste('./',paste(group,taxon_level, 'abundance_bar_plot.png', sep = '_'), sep = '')),
             plot = gg,
             device = 'png',
             units = 'in',
             width = 14,
             height = 8, dpi = 600)
    }
  }else{
    #box plotting
    #calculate the relative abundance of the taxa for each sample row
    otu_table_abund.df <- as.data.frame(t(apply(X = otu_table, MARGIN = 1, 
                                                FUN = function(x) x/sum(x))))
    
    if(group_rare_taxa){
      #grouup rare otus / taxa
      otu_table_abund.df <-  group_low_abund_taxa(abund_table = otu_table_abund.df, threshold = threshold)
    }
    common.ids <- intersect(rownames(independent_variable.df), rownames(otu_table_abund.df))
    independent_variable.df <- independent_variable.df[common.ids,,drop=FALSE]
    otu_table_abund.df <- otu_table_abund.df[common.ids,,drop=FALSE]
    #add the independent variable to the abundance table
    otu_table_abund.df <- cbind(independent_variable.df,otu_table_abund.df)
    #stack the data
    stacked_table<- Stack_data(Data =otu_table_abund.df, 
                               lab_columns = group, 
                               num_of_numericVariables = length(2:ncol(otu_table_abund.df )), 
                               col_num_of_numericVariables = 2:ncol(otu_table_abund.df ))
    
    #change the name of column from 'Pathogens' to the taxon level
    i <- which(colnames(stacked_table) == 'Pathogens')
    colnames(stacked_table)[i] <- taxon_level
    #convert the relative abundance to percentage if required
    if(convert_percent){
      stacked_table$Relative.abundance <- stacked_table$Relative.abundance *100
      Ylab <- 'Relative abundance (%)'
    }
    
    if(sort_box_by_abundance){
      #remove the independent variable column
      OTU_TABLE<- otu_table_abund.df[,-1] 
      #get the list of microbes ordered by relative abundance
      ordered_microbes_list <- sort_microbes_on_relative_abundance(otu_table = t(OTU_TABLE),
                                                                   ORDER = box_order)
      #check if we group rare taxa
      if ((grep(pattern = "Rare", x = ordered_microbes_list)) > 0){
        #temporarily remove the Rare group
        ordered_microbes_list  <- ordered_microbes_list [-(grep(pattern = "Rare", x = ordered_microbes_list))]
        #append the Rare group to the end of the list of microbes
        ordered_microbes_list <- c(ordered_microbes_list,"Rare")
      }
      #order the taxon level variable based on descending abundance
      stacked_table[,taxon_level] <- factor(x = stacked_table[,taxon_level], levels = ordered_microbes_list)
    }
    
    #box plot
    gg <- ggboxplot(data = stacked_table, 
                    x = taxon_level, 
                    y = 'Relative.abundance', 
                    color = group,
                    xlab = Xlab,
                    ylab = Ylab, 
                    add = 'mean_se', 
                    palette = palette)
    
    #plot the axes ticks inwards 
    gg <- gg +theme(axis.ticks.length=unit(-0.15, "cm"),
                    axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                    axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
    
    #format the legend and ticks labs
    gg <- ggpar(gg,font.x = c(18,'bold.italic','black'), font.tickslab =  c(16,'bold','black'),
                font.y = c(18,'bold.italic','black'),legend.title = box_legend_title, 
                legend = position, font.legend = c(14, "bold", "black"), x.text.angle = 45)
    if(write_plot){
      ggsave(filename = (paste('./',paste(group,taxon_level, 'abundance_box_plot.png', sep = '_'), sep = '')), 
             plot = gg, 
             device = 'png', 
             units = 'in', 
             width = 14, 
             height = 8, dpi = 600)
    }
  }
  
  
  final_result <- list(otu_table2,stacked_table, gg)
  names(final_result) <- c('unstacked_table', 'stacked_table','abundance_plot')
  return(final_result)
}


long2wide <- function(long_abund_table, row_names='samples', col_names="Phylum", values="relativeAbundance"){
  
  df <- pivot_wider(data = long_abund_table, id_cols = row_names,
                    names_from = col_names, values_from =values) %>% 
    as.data.frame(check.names=FALSE, stringAsFactor=FALSE)
  
  rownames(df) <- df[,row_names]
  
  col2drop <- which(colnames(df) == row_names)
  
  df <- df[,-col2drop]
  
  return(as.matrix(df))
}


make_area_plot <- function(stacked_data.df, x,  XLAB, 
                           group="Phylum", y="Relative.abundance", YLAB='Relative abundance (%)',
                           palette=custom_palette){
  
  # stacked_data.df - A long formatted dataframe with three columns - 
  # 1) independent-variable (x) 2) group for filling the area e.g microbe names 
  # 3) Relativae abundan vales for plotting
  # x - A string with the name of the variable (independent variable) on the x-axis
  # y - A string with the name of the variable (dependent variable) on the y-axis
  # group - A string specifying the name of the variable to be used for filling the area plot
  
  # x axis data must be numeric in order to make an area plot
  stacked_data.df[,x] <- as.numeric(stacked_data.df[,x])
  
  
  # Change area plot fill colors by groups
  # Note to use aes_string the strings must be double quoted
  ga<-ggplot(data = stacked_data.df, 
             aes_string(x=x,
                        y=y,
                        fill=group)) +
    # Use the data unchanged for plotting
    geom_area()+
    #add palette manually
    scale_fill_manual(values=custom_palette)+
    labs(x= XLAB, y=YLAB)+
    theme_classic()
  # Plot the axes ticks inwards 
  ga <- ga +theme(axis.ticks.length=unit(-0.15, "cm"),
                  axis.text.x=element_text(margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
                  axis.text.y=element_text(margin=margin(t=0,r=0.5,b=0,l=0,unit ="cm")))
  
  return(ga)
}


################################################################################
# A custom geometric mean function, with zero/NA tolerance.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
################################################################################
# phyloseq_CLR code
# Requires zCompositions package
################################################################################
zero_comp = function(x){
  if(taxa_are_rows(x)){x <- t(x)}
  matx = otu_table(x)
  # `zCompositions::cmultRepl` expects the samples to be in rows and OTUs to be in columns 
  matxzc = zCompositions::cmultRepl(matx, method="CZM", output="p-counts")
  otu_table(x) <- otu_table(matxzc, taxa_are_rows = FALSE)
  return(x)
}
# CLR definition
geometric_mean = function(x){
  exp(mean(log(x)))
}
clr = function(x, base=2){
  x <- log((x / geometric_mean(x)), base)
}
phyloseq_CLR = function(physeq){
  suppressMessages({physeq <- zero_comp(physeq)})
  return(transform_sample_counts(physeq, fun = clr))
}

# Function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
collapse_samples <- function(taxon_table,metadata,group,fun=sum, convertToRelativeAbundance=FALSE){
  # function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
  # taxon_table - a matrix count table with samples as rows and features/OTUs as columns
  # metadata - a dataframe to containing the group to collapse samples by. Sample names must be the rownames of the metadata
  # group - an independent factor variable within the metadata to collapse the samples by
  # fun - a function without brackets to apply in order to collapse the samples
  # convertToRelativeAbundance - a boolean set to TRUE OR FALSE if the taxon_table shout be converted to relative abundance
  # default is FALSE
  common.ids <- intersect(rownames(taxon_table),rownames(metadata))
  metadata <- droplevels(metadata[common.ids,,drop=FALSE])
  taxon_table <- taxon_table[common.ids,,drop=FALSE]
  taxon_table <- cbind(subset(x = metadata, select=group),taxon_table)
  
  taxon_table <- aggregate(reformulate(termlabels = group, response = '.'),
                           data = taxon_table, FUN = fun)
  rownames(taxon_table) <- taxon_table[,1]
  taxon_table <- taxon_table[,-1]
  if(convertToRelativeAbundance){
    taxon_table <- t(apply(X = taxon_table, MARGIN = 1, FUN = function(x) x/sum(x)))
  }
  
  final <- list(taxon_table,metadata)
  names(final) <- c("taxon_table","metadata")
  return(final)
}

summarise_group <- function(x_vect, y_vect, FUNC=mean, FUNC_string="mean", add_func_name=FALSE){
  # A function to return a one row summary dataframe of
  # y-vect grouped by x_vect using a specified function FUNC
  
  df <- data.frame(x=x_vect, y=y_vect)
  
  res <- aggregate(y~x, data=df, FUN=FUNC)
  rownames(res) <- res$x
  res <- res[,-1, drop=F] %>% t
  if(add_func_name){
    colnames(res) <- paste0(colnames(res),sep=sprintf("_%s",FUNC_string))
  }
  return(res)
}

summarise_groups <- function(long_df, metrics, metrics_column,  group, y_col){
  
  # Example :- summarise_groups(diversity_long.df, diversity_metrics,"measures", "BioFilter", "value")
  # long_df - A long formatted dataframe
  # metrics - A vector of subgroups in metrics_column e.g c('Observed', 'Shannon')
  # metrics_column - A string specifying the column with separate groups describe in metrics
  # group - A string specifying 
  
  df <- map_dfr(.x = metrics, .f = function(metric=.x){
    data2use <- long_df %>% dplyr::filter(!!sym(metrics_column) ==  metric) %>% as.data.frame(check.names=FALSE)
    data.frame(Metric=metric, summarise_group(data2use[,group], data2use[, y_col]), check.names = FALSE)
  }
  )
  
  rownames(df) <- df[,"Metric"]
  
  return(df)
  
}

foldchange <- function(group1, group2, metric, mean_table){
  # the function must be run only after genetrating the mean table
  # using the summarise_groups function
  
  # The ratio of the changes between final and intial over initial
  # i.e the ratio of the changes between final and intial over initial
  group1_val <- mean_table[metric, group1]
  group2_val <- mean_table[metric, group2]
  #res <- (group1_val - group2_val)/group2_val
  res <- group1_val/group2_val
  
  return(res)
  
}

add_foldchange_to_pairwise_table <- function(df, metric_column, group1_column, group2_column, mean_table){
  
  df <- as.data.frame(df, check.names=FALSE)
  
  df %>% 
    dplyr::mutate(foldChange=pmap_dbl(.l = list(metric=df[,metric_column], 
                                                group1=df[,group1_column], 
                                                group2=df[,group2_column]), 
                                      .f = function(metric, group1, group2){
                                        foldchange(group1 = group1, group2 = group2, 
                                                   metric = metric, 
                                                   mean_table = mean_table)
                                        
                                      } )) %>%  
    dplyr::select(Metric, group1, group2, foldChange, everything())
  
}

# Get a list of named vectors of letters representing
# significant differences between treatments
get_letters_vec <- function(metrics, comp_df) {
  
  #Example:
  # metrics - a vector of metrics within the column "Metric" e.g "Shannon", "pH" etc
  # comp_df - a dataframe generated after running the function "pairwise_test"   
  map(.x = metrics, function(metric=.x){
    
    sub_comp <- comp_df %>% filter(Metric == metric)
    
    # Retrieve the p-values for each pairwise comparison
    p_values <- sub_comp$p
    names(p_values) <- paste(sub_comp$group1,sub_comp$group2, sep = "-")
    
    # get the letters
    multcompView::multcompLetters(p_values)$Letters
    
    
  })
  
}



group_bar_plot <- function(df, x, y, groups, palette, x_lab, y_lab, 
                           bar_color, facet_plot= TRUE, summarise_by="mean",
                           text_size=6,
                           facet_by=NULL, letters_vec=NULL){
  # letters_vec - a vector with Treatments as it names and letters as its values
  # groups must be arranges in the order of the treatments i.e group1:group2:groupi
  # Example :
  # group_bar_plot(df=no_control_diversity, x="Event", y="Observed",
  # groups=c("Biofilter_media", "Biofilter_configuration", "Event"),
  # palette=annotation_colors$Biofilter_media,
  # x_lab = "Event", y_lab="Observed", bar_color="Biofilter_media",
  # facet_by=reformulate("Biofilter_configuration"))
  
  
  
  
  summary_table <- desc_statby(df, measure.var = y, grps = groups) %>% 
    pavian::zero_if_na() %>%
    droplevels()
  
  e <- try(summary_table %>% unite(col = "Treatment", 
                                   all_of(groups),
                                   sep = ":", remove = FALSE), 
           silent = TRUE)
  
  if(!inherits(e, "try-error")){
    
    summary_table <-  summary_table  %>%  unite(col = "Treatment", 
                                                all_of(groups),
                                                sep = ":", remove = FALSE)
    
  }
  
  
  if(!is.null(letters_vec)){
    
    labels <- letters_vec %>%
      as.data.frame() %>%
      set_names(., "labels") %>%
      rownames_to_column("Treatment")
    
    summary_table <-  summary_table %>% left_join(labels)
  }
  
  
  if(summarise_by != "median"){
    p <- ggplot(data = summary_table, 
                mapping = aes_string(x=x, y="mean", fill=bar_color))+ 
      geom_bar(position="dodge", stat="identity") + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                    width=0.2, position=position_dodge(0.9)) 
  }else{
    
    p <- ggplot(data = summary_table,
                mapping = aes_string(x=x, y="median", fill=bar_color))+ 
      geom_bar(position="dodge", stat="identity") + 
      geom_errorbar(aes(ymin=median-mad, ymax=median+mad),
                    width=0.2, position=position_dodge(0.9)) 
    
  }
  p <- p + 
    scale_fill_manual(values = palette)
  
  # add text to plot
  if(!is.null(letters_vec)){
    
    toAdd <- median(summary_table$range) - min(summary_table$range)
    max_label_len <- str_length(summary_table$labels) %>% max()
    if(max_label_len  < 3){
      
      p <- p + geom_text(mapping = aes(y=mean+se+toAdd, 
                                       label=labels,
                                       fontface = "bold"),
                         size = text_size, position=position_dodge(0.9))
      
    }else{
      text_size <- 4
      p <- p + ggrepel::geom_text_repel(mapping = aes(y=mean+se+toAdd, 
                                                      label=labels,
                                                      fontface = "bold"),
                                        size = text_size, position=position_dodge(0.9))
      
    }
    
    
  }
  
  if(facet_plot && (!is.null(facet_by)) ){ p <- p + 
    facet_grid(facet_by, drop = TRUE, scale="free") }
  
  legend_title <-  bar_color %>% str_replace_all('_',' ') %>% toTitleCase 
  p <- p +
    labs(x=x_lab, y=y_lab, fill=legend_title) +
    publication_format
  
  return(p)
  
}

##### Perform constrained analysis either dbRDA or CCA ######################

# Function to generate an ordered list of microbes from an oTU table based 
# on relative abundance
sort_microbes_on_relative_abundance <- function(otu_table, 
                                                ORDER= 'decreasing') {
  
  #this function generates an ordered list of microbes from an oTU table based on relative abundance
  #inputs
  #1. otu table - with samples as columns and microbes as row names  
  #2. ORDER is a string that species the order of the microbes based on relative abundance
  #the options is 'decreasing' or 'increasing' for decreasing and increase list based on relative abundance respectively
  
  #ensure the OTU table is a dataframe
  otu_table<-as.data.frame(otu_table)
  #re-arranging the OTU table
  colNames<- colnames(otu_table) #get the colnames of the otu table
  otu_table$microbes <- rownames(otu_table) #create a new column that will contaisn the names of the microbes
  
  #re-order the columns of the otutable
  col.order <- c("microbes", colNames)
  otu_table <- otu_table[,col.order]
  
  #stack the otu table by microbes
  stacked_table <- Stack_data(Data = otu_table, 
                              lab_columns = 'microbes', 
                              num_of_numericVariables = length(2:ncol(otu_table)), 
                              col_num_of_numericVariables = 2:ncol(otu_table))
  
  
  colnames(stacked_table) <- c("microbes","samples","relative_abundance")
  
  #calculate the mean abundance for each microbe
  microbes_mean<-aggregate(.~microbes,data = stacked_table, FUN = mean)
  if(ORDER == 'decreasing'){
    #order the dataframe by microbe relative abundance from the highest to the lowest
    ordered_microbes <-microbes_mean[order(microbes_mean$relative_abundance, decreasing = TRUE),]
  }else if(ORDER == 'increasing'){
    ordered_microbes <-microbes_mean[order(microbes_mean$relative_abundance, decreasing = FALSE),]
  }else{
    #cat("You have not selected a valid order either ORDER = 'decreasing' or 'increasing' ")
    stop("You have not selected a valid order either ORDER = 'decreasing' or 'increasing' ")
  }
  #get the list of microbes with decreasing abundance
  ordered_microbes_list <- ordered_microbes[,'microbes']
  
  return(ordered_microbes_list)
}
make_otu_table_taxa_name_short <- function(otu_table,split_by='_',
                                           taxon_level='genus'){
  #1. otu_table - the absolute rarefied output taxon table file produced by running summarize_taxa.py with
  #the -a option on the rarefied OTU generated by QIIME after core diversity analysis
  #2. split_by - a string specifying the delimiter for spliting the taxa names either "_" or ";"
  #3. taxon_level - a string specifying the label for the taxon level under consideration
  
  #note the microbe namse will be adjust in the final output for uniqueness if duplicate names exists
  #split the names of the taxonomy assignment
  if (split_by == '_'){
    split_names <- strsplit(colnames(otu_table), split = '_')
  }else{
    split_names <- strsplit(colnames(otu_table), split = ';')
  }
  #intialize an empty vector that will contain the formated microbe name
  microbe <- c()
  #loop over the split_names list of lists getting the last element of each list
  #and setting the respective microbe_name value to the rightful index in the vector
  for (index in 1:length(colnames(otu_table))) {
    k <- split_names[[index]]
    microbe[index] = split_names[[index]][length(k)]
  }
  
  #convert the microbe vector of names to a dataframe
  
  microbe <- as.data.frame(microbe)
  colnames(microbe) <- taxon_level #rename the column
  
  #rename the columns of the otu table as the formatted microbes name
  colnames(otu_table)<- microbe[,taxon_level]
  
  
  #convert the oTU table to a dataframe - this trick helps with getting unique microbe name
  #because it forces the column names to be unique
  #this trick allows for unique names of the columns when
  #we have duplicated names like Ambigous taxa, other and uncultured bacterium
  otu_table2 <- data.frame(otu_table)
  return(otu_table2)
}
Stack_data <- function(Data,lab_columns,num_of_numericVariables,
                       col_num_of_numericVariables) {
  #data is the dataframe to stack
  #lab_colums is a vector of numbers or strings specifying the columns 
  #containing the labels or independent variables
  
  #num_of_numericVariables is a number specifying how many 
  #numeric / denpendent variables you want to stack 
  
  #col_num_of_numericVariables is a vector of numbers or strings specifying the columns 
  #containing the numeric or dependent variables 
  
  #subsetting the non-numeric variable columns
  lab <- subset(x = Data, select = lab_columns)
  
  #initialize an empty dataframe
  tab <- data.frame()
  #bind the non-numeric columns row-wise 18 times
  for (i in 1:num_of_numericVariables){
    tab <- rbind(tab,lab)
  }
  #stack the numeric variable columns together
  j <- stack(Data[,col_num_of_numericVariables])
  #rename the column names
  colnames(j) <- c("Relative.abundance", "Pathogens")
  #re-arrange the columns
  k <- data.frame(Pathogens = j$Pathogens, Relative.abundance = j$Relative.abundance)
  #combine the stacked variable to form a new stacked dataframe
  stacked_data2 <- cbind(tab,k)
  colnames(stacked_data2) <- c(lab_columns,colnames(k))
  rownames(stacked_data2) <- 1:length(stacked_data2[,1])
  return(stacked_data2)
}

#' @title Get percent of total inertia explained by principal or constrained axes
#' @param mod cca.object
#' @param axes A numeric vector indicating the axes for which percent explained inertia should be returned for
#' @return Returns a named vector containing explained inertia as percentages. Returns constrained inertia fo CCA and unconstrained for CA
axis.expl <- function(mod, axes = 1:2) {
  
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
  
}

constrained_analysis <- function(otu_table_file,table_type= 'OTU', 
                                 seq_chem_data.df,env_columns=6:13,
                                 fill, category_columns=2:5,shape,
                                 lab_shape,color,lab_color,
                                 legend_title=NULL,sort_by_abundance = FALSE, 
                                 palette = c('chocolate4','chocolate3',
                                             'chocolate1'),
                                 constrain_type='dbrda', add_species = FALSE,
                                 file_name, result_dir="./", save_figure=FALSE,
                                 size=NULL) {
  
  # This function performs constrained ordinations c("dbRDA" or "CCA") and 
  # returns the plot along with the significance statistics
  
  # EXAMPLE
  # result <- constrained_analysis(otu_table_file = feature_table, 
  #                                seq_chem_data.df = metadata,
  #                                fill = "Diet", shape = "TreatmentControl", 
  #                               lab_shape = "TreatmentControl",file_name = 'soiltype',
  #                                color = "Diet", lab_color = "Diet", result_dir = plot_dir,
  #                                env_columns = 7:22,category_columns = 3:6,
  #                                palette= annotation_colors[['Diet']])
  # 
  ###### Parameters #########
  # ~ otu_table_file - a string specifying the path to your otu table biom file
  # or a matrix with features as rows and samples as columns
  # ~ table_type - 
  # ~ seq_chem_data.df - is a dataframe with categorical variables to analysed and the coresponding
  #   environmental variables columns/variables
  # ~ env_columns - is a string or numeric vector indicating the columns in seq_chem_data.df
  #   containing the environmental variables
  # ~ category_columns - is a string or numeric vector indicating the columns in seq_chem_data.df
  #   containing the categorical / independent variables
  # ~ fill - what color or categorical variable should be used fill the points in your constrained
  #   ordination plot
  # ~ color - what color or categorical variable should be used around the points in your constrained
  #   ordination plot
  # ~ lab_color _ a string specifying the name to be used in the legend for the variable specifying the 
  #   color of your points
  # ~ shape - a string specifying the categorical variable / column in to be represented as the sahpe of 
  #   of your points in the ordination plot
  # ~ lab_shape _ a string specifying the name to be used in the legend for the variable specifying the 
  #   shape of your points
  # ~ constrain_type - the type of contrained analysis to perform either 'dbrda' for dbRDA or 'cca' for CCA
  # ~ palette - a vector of color palettes to use for plotting
  # ~ sort_by_abundance  - Either true or false to get only the most abundant genera (20 by default)
  # ~ add_species  - a logical specipying if a triplot with species included should be made
  # ~ file_name - a string with the .png extension specifying the name of the png file to save 
  #   the ordination as
  # ~ result_dir -  a directory to output the results of the analysis 
  # ~ save_figure - should the ordination be saved? TRUE OR  FALSE
  
  
  #load in the required libraries
  library(vegan) 
  library(biomformat)
  library(ggpubr)
  
  if(is.character(otu_table_file)){
    #read-in the bion table file    
    otu_table.biom<- read_biom(biom_file = otu_table_file)
    #convert the table to a matrix then transpose
    otu_table.m <- t(as.matrix(biom_data(otu_table.biom)))
  }else{
    
    otu_table.m <- t(otu_table_file)
  }
  
  
  
  #if the table supplied is an OTU table not a taxon level table do this
  if (table_type == 'OTU'){
    #convert to relative abundance
    otu_table.m <- t(apply(X = t(otu_table.m),MARGIN = 2, FUN = function(x) x/sum(x)))
    otu_table.df <- as.data.frame(otu_table.m)
  }
  
  #if the table supplied is a taxon table do this
  if (table_type != 'OTU'){
    otu_table.df <- make_otu_table_taxa_name_short(otu_table = otu_table.m,taxon_level = 'genus')
    
    if(sort_by_abundance == TRUE){
      #get only the top twenty most abundant genera
      abund_genera <- sort_microbes_on_relative_abundance(otu_table = t(otu_table.df))[1:20]
      otu_table.df <- otu_table.df[,abund_genera]
    }
    
  }
  
  
  
  #get the common sample names
  common.ids <- intersect(rownames(seq_chem_data.df),rownames(otu_table.df))
  #subset the chemical dataframe to contain only the data for the samples present in the otu_table 
  seq_chem_data.df <- seq_chem_data.df[common.ids,]
  #subset the otu matrix to contain only the data for the samples present in the chemical dataframe
  otu_table.df <- otu_table.df[common.ids,] 
  #get only the environmental variables columns
  meta_table <- seq_chem_data.df[,env_columns]
  
  
  #Use adonis to find significant environmental variables
  otu_table.adonis <- adonis(otu_table.df ~ ., data=meta_table)
  
  #Check if there are significant environmental variables
  CHECK <- sum(na.omit(otu_table.adonis$aov.tab$"Pr(>F)"<=0.05))
  if( CHECK == 0) { print("None of the environmetal variables were significant") }
  
  
  #Check if there are significant environmental variables and if you want to perform CCA analysis
  if((CHECK > 0) && (constrain_type =='cca')){
    
    #Extract the best variables
    bestEnvVariables<-rownames(otu_table.adonis$aov.tab)[otu_table.adonis$aov.tab$"Pr(>F)"<=0.05]
    
    #Last two are NA entries, so we have to remove them
    bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
    #set the x and y variable strings for plotting
    x_axis <- 'CCA1'
    y_axis <- 'CCA2'
    #We are now going to use only those environmental variables in cca that were found significant
    eval(parse(text=paste("sol <- cca(otu_table.df ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=meta_table)",sep="")))
    
  }else{
    #do this if dbRDA contraint ordination is the selected option 
    if(constrain_type =='dbrda'){
      #set the x and y variable strings for plotting
      x_axis <- 'dbRDA1'
      y_axis <- 'dbRDA2'
      
      sol <-dbrda(otu_table.df ~ .,data=meta_table, distance = "bray")
      
    }else{
      #set the x and y variable strings for plotting
      x_axis <- 'CCA1'
      y_axis <- 'CCA2'
      #for all the chemical data
      sol<-cca(otu_table.df ~ ., data=meta_table)
      
    }
    
    
  }
  
  #get the percentage of variance explained by each contrained axis to be used for plotting
  labs <- axis.expl(sol)
  
  #get the scores for species, biplots e.t.c
  scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
  
  #write out the significance of the CCA to a file
  sink(file = glue::glue('{result_dir}/CCA_significance.txt'))
  
  print('Model Significance')
  #calculate the significance of your CCA model
  print(model_significance <- anova.cca(sol))
  print('axes Significance')
  #calculate the significance of your axes
  print(axes_significance <- anova.cca(sol,by = "axis",permutations = 500))
  #calculate the significance of your environmental variables
  print("Environmental Variables Significance")
  print(env_significance <- anova.cca(sol,by = "terms",permutations = 200))
  ## calculating loadings/environmental correlations with the axes
  site_scores <- scrs$sites
  site_scores_environment<-cbind(site_scores,meta_table) ## merge
  print("Correlation of the environmal variables with the constrained axes")
  print(correlations <- cor(site_scores_environment)) ## calculate correlations
  
  sink()
  
  
  
  
  #Extract site(samples) data
  df_sites<-data.frame(scrs$sites)
  
  #add treatment information to the samples data
  treatment_info <- subset(x=seq_chem_data.df, select=category_columns)
  df_sites <- cbind(treatment_info,df_sites)
  
  
  ####### Plotting #########
  
  #Draw the  sites / samples
  p <- ggscatter(data = df_sites,
                 x =  x_axis,
                 y =  y_axis,
                 fill = fill,
                 color = color ,
                 size = ifelse(test = is.null(size),yes = 5,no =  size),
                 palette = palette,
                 xlab = paste0(names(labs[1]), " (", sprintf("%.1f", labs[1]), "%)"),
                 ylab = paste0(names(labs[2]), " (", sprintf("%.1f", labs[2]), "%)"),
                 shape = shape)
  
  
  # Draw biplots with arrows indicating the environmental variables
  multiplier <- vegan:::ordiArrowMul(scrs$biplot)
  df_arrows<- scrs$biplot*multiplier
  # Rename the columns as x and y
  colnames(df_arrows)<-c("x","y")
  # Create a data frame for the arrows
  df_arrows<-as.data.frame(df_arrows)
  
  # Get a copy of the significance of the environmental variables
  env<- env_significance 
  
  # Get the environmental variables that are common to both dataframes
  common_envs <- intersect(rownames(df_arrows),rownames(env))
  # Subset the dataframes to contain the environmental variables that are common to both
  df_arrows <- df_arrows[common_envs,]
  env <- env[common_envs,]
  # Create a new column in the df_arrows dataframe that will hold the p-values for each
  # environmental variable and assign their p-values
  df_arrows$p_value <- env$`Pr(>F)`
  # Create a new column in the df_arrows dataframe that will hold the colors for each
  # environmental variable and assign their colors red for significant variables and 
  # black for insignificant variables
  df_arrows$arrow_color <- ''
  df_arrows[which(df_arrows$p_value < 0.05),'arrow_color'] <- 'red' #significant envs
  df_arrows[-(which(df_arrows$p_value < 0.05)),'arrow_color'] <- 'black' #insignificant envs
  arrow_color<- df_arrows$arrow_color
  
  
  if( CHECK == 0) { 
    p <- p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                        arrow = arrow(length = unit(0.2, "cm")),color='black',alpha=1)
    
    # Add the environmental variable name to the arrow
    p <-p+geom_text(data=as.data.frame(df_arrows[,c("x","y")]*1.1),
                    aes(x, y, label = rownames(df_arrows)),
                    color='black',
                    alpha=1)
  }else{
    
    # Draw the arrows on the existing plot
    p <- p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                        arrow = arrow(length = unit(0.2, "cm")),color=arrow_color,alpha=1)
    # Add the environmental variable name to the arrow
    p<- p+geom_text(data=as.data.frame(df_arrows[,c("x","y")]*1.1),
                    aes(x, y, label = rownames(df_arrows)),
                    color=arrow_color,
                    alpha=1)
  }
  
  
  # # Draw species
  if(add_species){
    df_species<- as.data.frame(scrs$species)
    colnames(df_species)<-c("x","y")
    
    # Either choose text or points
    p<- ggtext(data =df_species, x = 'x', y = 'y',label = rownames(df_species), ggp = p,size = 14)
    #  The shape arguement here is just for name purpose
    # the scale_shape_manual defines the type shape to plotted in this case 2 represents triangles
    p <-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2)
    p <-p+theme_bw()
  }
  
  # Format the figure font sizes and legend
  p <- ggpar(p = p,
             font.x = c(20, "bold", "black"),
             font.y = c(20, "bold", "black"),
             font.legend = c(18, "bold", "black"), 
             font.tickslab = c(18, "bold", "black"), 
             legend = 'right')
  
  # Label the legends
  p <- p + labs(shape=lab_shape, colour=lab_color,fill = lab_color)
  
  if(save_figure){
    # save the plot to a png file
    png(filename = glue::glue("{result_dir}/{file_name}_{constrain_type}.png"),
        width = 14,
        height = 8,
        units = "in",
        res = 600)
    
    print(p)
    dev.off()
  }
  
  #save the significance results to be ouputed and name them in the list
  significance <- list(model_significance,axes_significance, env_significance,correlations)
  names(significance) <- c('model','axes','envs','correlations')
  #add the significance and the graph in the final list
  final <- list(p,sol,significance)
  names(final) <- c('plot','cca_object','significance_tests')
  #return the final list
  return(final)
}


# Step wise VIF function for detecting and removing multicollinearity

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
   # Example:
   # https://www.r-bloggers.com/2013/02/collinearity-and-stepwise-vif-selection/
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}


# Modified correlation plot function from the corrplot package


draw_grid = function(coords, fg) {
  symbols(coords, add = TRUE, inches = FALSE, fg = fg, bg = NA,
          rectangles = matrix(1, nrow = nrow(coords), ncol = 2))
}

corrplot2  <- function (corr, method = "circle", 
                        type = "full", add = FALSE, 
                        col = NULL, col.lim = NULL, bg = "white", title = "", 
                        is.corr = TRUE, diag = TRUE, outline = FALSE, mar = c(0, 
                                                                              0, 0, 0),
                        addgrid.col = NULL, addCoef.col = NULL, addCoefasPercent = FALSE, 
                        order = "original", hclust.method ="complete", 
                        addrect = NULL, rect.col = "black", rect.lwd = 2, tl.pos = NULL, 
                        tl.cex = 1, tl.col = "red", tl.offset = 0.4, tl.srt = 90, 
                        cl.pos = NULL, cl.length = NULL, cl.cex = 0.8, cl.ratio = 0.15, 
                        cl.align.text = "c", cl.offset = 0.5, number.cex = 1, 
                        number.font = 2, number.digits = NULL, addshade = "all", shade.lwd = 1, shade.col = "white", 
                        p.mat = NULL, sig.level = 0.05, insig = "pch", pch = 4, 
                        pch.col = "black", pch.cex = 3, plotCI = "n", lowCI.mat = NULL, 
                        uppCI.mat = NULL, na.label = "?", na.label.col = "black", 
                        win.asp = 1, ...) 
{
  # method = match.arg(method)
  # type = match.arg(type)
  # order = match.arg(order)
  # hclust.method = match.arg(hclust.method)
  # addshade = match.arg(addshade)
  # insig = match.arg(insig)
  # plotCI = match.arg(plotCI)
  if (win.asp != 1 && !(method %in% c("circle", "square"))) {
    stop("Parameter 'win.asp' is supported only for circle and square methods.")
  }
  asp_rescale_factor = min(1, win.asp)/max(1, win.asp)
  stopifnot(asp_rescale_factor >= 0 && asp_rescale_factor <= 
              1)
  if (!is.matrix(corr) && !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  if (is.null(addgrid.col)) {
    addgrid.col = switch(method, color = NA, shade = NA, 
                         "grey")
  }
  if (any(corr[!is.na(corr)] < col.lim[1]) || any(corr[!is.na(corr)] > 
                                                  col.lim[2])) {
    stop("color limits should cover matrix")
  }
  if (is.null(col.lim)) {
    if (is.corr) {
      col.lim = c(-1, 1)
    }
    else {
      if (!diag) {
        diag(corr) = NA
      }
      col.lim = c(min(corr, na.rm = TRUE), max(corr, na.rm = TRUE))
    }
  }
  SpecialCorr = 0
  if (is.corr) {
    if (min(corr, na.rm = TRUE) < -1 - .Machine$double.eps^0.75 || 
        max(corr, na.rm = TRUE) > 1 + .Machine$double.eps^0.75) {
      stop("The matrix is not in [-1, 1]!")
    }
    SpecialCorr = 1
    if (col.lim[1] < -1 | col.lim[2] > 1) {
      stop("col.lim should be within the interval [-1,1]")
    }
  }
  intercept = 0
  zoom = 1
  if (!is.corr) {
    c_max = max(corr, na.rm = TRUE)
    c_min = min(corr, na.rm = TRUE)
    if (diff(col.lim)/(c_max - c_min) > 2) {
      warning("col.lim interval too wide, please set a suitable value")
    }
    if (c_max <= 0 | c_min >= 0) {
      intercept = -c_min
      zoom = 1/(diff(col.lim))
      if (col.lim[1] * col.lim[2] < 0) {
        warning("col.lim interval not suitable to the matrix")
      }
    }
    else {
      stopifnot(c_max * c_min < 0)
      stopifnot(c_min < 0 && c_max > 0)
      intercept = 0
      zoom = 1/max(abs(col.lim))
      SpecialCorr = 1
    }
    if (zoom == Inf) {
      stopifnot(c_max == 0 && c_min == 0)
      zoom = 0
    }
    corr = (intercept + corr) * zoom
  }
  col.lim2 = (intercept + col.lim) * zoom
  int = intercept * zoom
  if (is.null(col)) {
    col = colorRampPalette(c("#67001F", "#B2182B", 
                             "#D6604D", "#F4A582", "#FDDBC7", 
                             "#FFFFFF", "#D1E5F0", "#92C5DE", 
                             "#4393C3", "#2166AC", "#053061"))(200)
  }
  n = nrow(corr)
  m = ncol(corr)
  min.nm = min(n, m)
  ord = seq_len(min.nm)
  if (order != "original") {
    ord = corrMatOrder(corr, order = order, hclust.method = hclust.method)
    corr = corr[ord, ord]
  }
  if (is.null(rownames(corr))) {
    rownames(corr) = seq_len(n)
  }
  if (is.null(colnames(corr))) {
    colnames(corr) = seq_len(m)
  }
  apply_mat_filter = function(mat) {
    x = matrix(1:n * m, nrow = n, ncol = m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) < 
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) = Inf
    }
    return(mat)
  }
  getPos.Dat = function(mat) {
    tmp = apply_mat_filter(mat)
    Dat = tmp[is.finite(tmp)]
    ind = which(is.finite(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    PosName = ind
    PosName[, 1] = colnames(mat)[ind[, 2]]
    PosName[, 2] = rownames(mat)[ind[, 1]]
    return(list(Pos, Dat, PosName))
  }
  getPos.NAs = function(mat) {
    tmp = apply_mat_filter(mat)
    ind = which(is.na(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    return(Pos)
  }
  testTemp = getPos.Dat(corr)
  Pos = getPos.Dat(corr)[[1]]
  PosName = getPos.Dat(corr)[[3]]
  if (any(is.na(corr)) && is.character(na.label)) {
    PosNA = getPos.NAs(corr)
  }
  else {
    PosNA = NULL
  }
  AllCoords = rbind(Pos, PosNA)
  n2 = max(AllCoords[, 2])
  n1 = min(AllCoords[, 2])
  nn = n2 - n1
  m2 = max(AllCoords[, 1])
  m1 = min(AllCoords[, 1])
  mm = max(1, m2 - m1)
  expand_expression = function(s) {
    ifelse(grepl("^[:=$]", s), parse(text = substring(s, 
                                                      2)), s)
  }
  newrownames = sapply(rownames(corr)[(n + 1 - n2):(n + 1 - 
                                                      n1)], expand_expression)
  newcolnames = sapply(colnames(corr)[m1:m2], expand_expression)
  DAT = getPos.Dat(corr)[[2]]
  len.DAT = length(DAT)
  rm(expand_expression)
  assign.color = function(dat = DAT, color = col, isSpecialCorr = SpecialCorr) {
    if (isSpecialCorr) {
      newcorr = (dat + 1)/2
    }
    else {
      newcorr = dat
    }
    newcorr[newcorr <= 0] = 0
    newcorr[newcorr >= 1] = 1 - 1e-16
    color[floor(newcorr * length(color)) + 1]
  }
  col.fill = assign.color()
  isFALSE = function(x) identical(x, FALSE)
  isTRUE = function(x) identical(x, TRUE)
  if (isFALSE(tl.pos)) {
    tl.pos = "n"
  }
  if (is.null(tl.pos) || isTRUE(tl.pos)) {
    tl.pos = switch(type, full = "lt", lower = "ld", 
                    upper = "td")
  }
  if (isFALSE(cl.pos)) {
    cl.pos = "n"
  }
  if (is.null(cl.pos) || isTRUE(cl.pos)) {
    cl.pos = switch(type, full = "r", lower = "b", 
                    upper = "r")
  }
  if (isFALSE(outline)) {
    col.border = col.fill
  }
  else if (isTRUE(outline)) {
    col.border = "black"
  }
  else if (is.character(outline)) {
    col.border = outline
  }
  else {
    stop("Unsupported value type for parameter outline")
  }
  oldpar = par(mar = mar, bg = par()$bg)
  on.exit(par(oldpar), add = TRUE)
  if (!add) {
    plot.new()
    xlabwidth = max(strwidth(newrownames, cex = tl.cex))
    ylabwidth = max(strwidth(newcolnames, cex = tl.cex))
    laboffset = strwidth("W", cex = tl.cex) * tl.offset
    for (i in 1:50) {
      xlim = c(m1 - 0.5 - laboffset - xlabwidth * (grepl("l", 
                                                         tl.pos) | grepl("d", tl.pos)), m2 + 0.5 + 
                 mm * cl.ratio * (cl.pos == "r") + xlabwidth * 
                 abs(cos(tl.srt * pi/180)) * grepl("d", 
                                                   tl.pos))
      ylim = c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b") - 
                 laboffset, n2 + 0.5 + laboffset + ylabwidth * 
                 abs(sin(tl.srt * pi/180)) * grepl("t", 
                                                   tl.pos) + ylabwidth * abs(sin(tl.srt * pi/180)) * 
                 (type == "lower") * grepl("d", tl.pos))
      plot.window(xlim, ylim, asp = 1, xaxs = "i", 
                  yaxs = "i")
      x.tmp = max(strwidth(newrownames, cex = tl.cex))
      y.tmp = max(strwidth(newcolnames, cex = tl.cex))
      laboffset.tmp = strwidth("W", cex = tl.cex) * 
        tl.offset
      if (max(x.tmp - xlabwidth, y.tmp - ylabwidth, laboffset.tmp - 
              laboffset) < 0.001) {
        break
      }
      xlabwidth = x.tmp
      ylabwidth = y.tmp
      laboffset = laboffset.tmp
      if (i == 50) {
        warning(c("Not been able to calculate text margin, ", 
                  "please try again with a clean new empty window using ", 
                  "{plot.new(); dev.off()} or reduce tl.cex"))
      }
    }
    if (.Platform$OS.type == "windows") {
      grDevices::windows.options(width = 7, height = 7 * 
                                   diff(ylim)/diff(xlim))
    }
    xlim = xlim + diff(xlim) * 0.01 * c(-1, 1)
    ylim = ylim + diff(ylim) * 0.01 * c(-1, 1)
    plot.window(xlim = xlim, ylim = ylim, asp = win.asp, 
                xlab = "", ylab = "", xaxs = "i", 
                yaxs = "i")
  }
  laboffset = strwidth("W", cex = tl.cex) * tl.offset
  symbols(Pos, add = TRUE, inches = FALSE, rectangles = matrix(1, 
                                                               len.DAT, 2), bg = bg, fg = bg)
  if (method == "circle" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = asp_rescale_factor * 
              0.9 * abs(DAT)^0.5/2, fg = col.border, bg = col.fill)
  }
  if (method == "ellipse" && plotCI == "n") {
    ell.dat = function(rho, length = 99) {
      k = seq(0, 2 * pi, length = length)
      x = cos(k + acos(rho)/2)/2
      y = cos(k - acos(rho)/2)/2
      cbind(rbind(x, y), c(NA, NA))
    }
    ELL.dat = lapply(DAT, ell.dat)
    ELL.dat2 = 0.85 * matrix(unlist(ELL.dat), ncol = 2, byrow = TRUE)
    ELL.dat2 = ELL.dat2 + Pos[rep(1:length(DAT), each = 100), 
    ]
    polygon(ELL.dat2, border = col.border, col = col.fill)
  }
  if (is.null(number.digits)) {
    number.digits = switch(addCoefasPercent + 1, 2, 0)
  }
  stopifnot(number.digits%%1 == 0)
  stopifnot(number.digits >= 0)
  if (method == "number" && plotCI == "n") {
    x = (DAT - int) * ifelse(addCoefasPercent, 100, 1)/zoom
    text(Pos[, 1], Pos[, 2], font = number.font, col = col.fill, 
         labels = format(round(x, number.digits), nsmall = number.digits), 
         cex = number.cex)
  }
  NA_LABEL_MAX_CHARS = 2
  if (is.matrix(PosNA) && nrow(PosNA) > 0) {
    stopifnot(is.matrix(PosNA))
    if (na.label == "square") {
      symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                                                               nrow(PosNA)), bg = na.label.col, fg = na.label.col)
    }
    else if (nchar(na.label) %in% 1:NA_LABEL_MAX_CHARS) {
      symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                                                               nrow(PosNA)), fg = bg, bg = bg)
      text(PosNA[, 1], PosNA[, 2], font = number.font, 
           col = na.label.col, labels = na.label, cex = number.cex, 
           ...)
    }
    else {
      stop(paste("Maximum number of characters for NA label is:", 
                 NA_LABEL_MAX_CHARS))
    }
  }
  if (method == "pie" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, 
                                                           len.DAT) * 0.85, fg = col.border)
    pie.dat = function(theta, length = 100) {
      k = seq(pi/2, pi/2 - theta, length = 0.5 * length * 
                abs(theta)/pi)
      x = c(0, cos(k)/2, 0)
      y = c(0, sin(k)/2, 0)
      cbind(rbind(x, y), c(NA, NA))
    }
    PIE.dat = lapply(DAT * 2 * pi, pie.dat)
    len.pie = unlist(lapply(PIE.dat, length))/2
    PIE.dat2 = 0.85 * matrix(unlist(PIE.dat), ncol = 2, byrow = TRUE)
    PIE.dat2 = PIE.dat2 + Pos[rep(1:length(DAT), len.pie), 
    ]
    polygon(PIE.dat2, border = "black", col = col.fill)
  }
  if (method == "shade" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = addgrid.col)
    shade.dat = function(w) {
      x = w[1]
      y = w[2]
      rho = w[3]
      x1 = x - 0.5
      x2 = x + 0.5
      y1 = y - 0.5
      y2 = y + 0.5
      dat = NA
      if ((addshade == "positive" || addshade == 
           "all") && rho > 0) {
        dat = cbind(c(x1, x1, x), c(y, y1, y1), c(x, 
                                                  x2, x2), c(y2, y2, y))
      }
      if ((addshade == "negative" || addshade == 
           "all") && rho < 0) {
        dat = cbind(c(x1, x1, x), c(y, y2, y2), c(x, 
                                                  x2, x2), c(y1, y1, y))
      }
      return(t(dat))
    }
    pos_corr = rbind(cbind(Pos, DAT))
    pos_corr2 = split(pos_corr, 1:nrow(pos_corr))
    SHADE.dat = matrix(na.omit(unlist(lapply(pos_corr2, shade.dat))), 
                       byrow = TRUE, ncol = 4)
    segments(SHADE.dat[, 1], SHADE.dat[, 2], SHADE.dat[, 
                                                       3], SHADE.dat[, 4], col = shade.col, lwd = shade.lwd)
  }
  if (method == "square" && plotCI == "n") {
    draw_method_square(Pos, DAT, asp_rescale_factor, col.border, 
                       col.fill)
  }
  if (method == "color" && plotCI == "n") {
    draw_method_color(Pos, col.border, col.fill)
  }
  draw_grid(AllCoords, addgrid.col)
  if (plotCI != "n") {
    if (is.null(lowCI.mat) || is.null(uppCI.mat)) {
      stop("Need lowCI.mat and uppCI.mat!")
    }
    if (order != "original") {
      lowCI.mat = lowCI.mat[ord, ord]
      uppCI.mat = uppCI.mat[ord, ord]
    }
    pos.lowNew = getPos.Dat(lowCI.mat)[[1]]
    lowNew = getPos.Dat(lowCI.mat)[[2]]
    pos.uppNew = getPos.Dat(uppCI.mat)[[1]]
    uppNew = getPos.Dat(uppCI.mat)[[2]]
    if (!method %in% c("circle", "square")) {
      stop("Method shoud be circle or square if drawing confidence intervals.")
    }
    k1 = (abs(uppNew) > abs(lowNew))
    bigabs = uppNew
    bigabs[which(!k1)] = lowNew[!k1]
    smallabs = lowNew
    smallabs[which(!k1)] = uppNew[!k1]
    sig = sign(uppNew * lowNew)
    color_bigabs = col[ceiling((bigabs + 1) * length(col)/2)]
    color_smallabs = col[ceiling((smallabs + 1) * length(col)/2)]
    if (plotCI == "circle") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(bigabs)^0.5/2, 
              bg = ifelse(sig > 0, col.fill, color_bigabs), 
              fg = ifelse(sig > 0, col.fill, color_bigabs))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(smallabs)^0.5/2, 
              bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                                                                      0, col.fill, color_smallabs))
    }
    if (plotCI == "square") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(bigabs)^0.5, bg = ifelse(sig > 
                                                                       0, col.fill, color_bigabs), fg = ifelse(sig > 
                                                                                                                 0, col.fill, color_bigabs))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(smallabs)^0.5, 
              bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                                                                      0, col.fill, color_smallabs))
    }
    if (plotCI == "rect") {
      rect.width = 0.25
      rect(pos.uppNew[, 1] - rect.width, pos.uppNew[, 2] + 
             smallabs/2, pos.uppNew[, 1] + rect.width, pos.uppNew[, 
                                                                  2] + bigabs/2, col = col.fill, border = col.fill)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + DAT/2, pos.lowNew[, 1] + rect.width, pos.lowNew[, 
                                                                                                             2] + DAT/2, col = "black", lwd = 1)
      segments(pos.uppNew[, 1] - rect.width, pos.uppNew[, 
                                                        2] + uppNew/2, pos.uppNew[, 1] + rect.width, 
               pos.uppNew[, 2] + uppNew/2, col = "black", 
               lwd = 1)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + lowNew/2, pos.lowNew[, 1] + rect.width, 
               pos.lowNew[, 2] + lowNew/2, col = "black", 
               lwd = 1)
      segments(pos.lowNew[, 1] - 0.5, pos.lowNew[, 2], 
               pos.lowNew[, 1] + 0.5, pos.lowNew[, 2], col = "grey70", 
               lty = 3)
    }
  }
  if (!is.null(addCoef.col) && method != "number") {
    text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = round((DAT - 
                                                                  int) * ifelse(addCoefasPercent, 100, 1)/zoom, number.digits), 
         cex = number.cex, font = number.font)
  }
  if (!is.null(p.mat) && insig != "n") {
    if (order != "original") {
      p.mat = p.mat[ord, ord]
    }
    if (!is.null(rownames(p.mat)) | !is.null(rownames(p.mat))) {
      if (!all(colnames(p.mat) == colnames(corr)) | !all(rownames(p.mat) == 
                                                         rownames(corr))) {
        warning("p.mat and corr may be not paired, their rownames and colnames are not totally same!")
      }
    }
    pos.pNew = getPos.Dat(p.mat)[[1]]
    pNew = getPos.Dat(p.mat)[[2]]
    if (insig == "label_sig") {
      if (!is.character(pch)) 
        pch = "*"
      place_points = function(sig.locs, point) {
        text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
             labels = point, col = pch.col, cex = pch.cex, 
             lwd = 2)
      }
      if (length(sig.level) == 1) {
        place_points(sig.locs = which(pNew < sig.level), 
                     point = pch)
      }
      else {
        l = length(sig.level)
        for (i in seq_along(sig.level)) {
          iter = l + 1 - i
          pchTmp = paste(rep(pch, i), collapse = "")
          if (i == length(sig.level)) {
            locs = which(pNew < sig.level[iter])
            if (length(locs)) {
              place_points(sig.locs = locs, point = pchTmp)
            }
          }
          else {
            locs = which(pNew < sig.level[iter] & pNew > 
                           sig.level[iter - 1])
            if (length(locs)) {
              place_points(sig.locs = locs, point = pchTmp)
            }
          }
        }
      }
    }
    else {
      ind.p = which(pNew > sig.level)
      p_inSig = length(ind.p) > 0
      if (insig == "pch" && p_inSig) {
        points(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
               pch = pch, col = pch.col, cex = pch.cex, lwd = 2)
      }
      if (insig == "p-value" && p_inSig) {
        text(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
             round(pNew[ind.p], number.digits), col = pch.col)
      }
      if (insig == "blank" && p_inSig) {
        symbols(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
                inches = FALSE, squares = rep(1, length(pos.pNew[, 
                                                                 1][ind.p])), fg = addgrid.col, bg = bg, add = TRUE)
      }
    }
  }
  if (cl.pos != "n") {
    colRange = assign.color(dat = col.lim2)
    ind1 = which(col == colRange[1])
    ind2 = which(col == colRange[2])
    colbar = col[ind1:ind2]
    if (is.null(cl.length)) {
      cl.length = ifelse(length(colbar) > 20, 11, length(colbar) + 
                           1)
    }
    labels = seq(col.lim[1], col.lim[2], length = cl.length)
    if (cl.pos == "r") {
      vertical = TRUE
      xlim = c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
      ylim = c(n1 - 0.5, n2 + 0.5)
    }
    if (cl.pos == "b") {
      vertical = FALSE
      xlim = c(m1 - 0.5, m2 + 0.5)
      ylim = c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn * 
                 0.02)
    }
    colorlegend(colbar = colbar, labels = round(labels, 2), 
                offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex, 
                xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }
  if (tl.pos != "n") {
    pos.xlabel = cbind(m1:m2, n2 + 0.5 + laboffset)
    pos.ylabel = cbind(m1 - 0.5, n2:n1)
    if (tl.pos == "td") {
      if (type != "upper") {
        stop("type should be 'upper' if tl.pos is 'dt'.")
      }
      pos.ylabel = cbind(m1:(m1 + nn) - 0.5, n2:n1)
    }
    if (tl.pos == "ld") {
      if (type != "lower") {
        stop("type should be 'lower' if tl.pos is 'ld'.")
      }
      pos.xlabel = cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
    }
    if (tl.pos == "d") {
      pos.ylabel = cbind(m1:(m1 + nn) - 0.5, n2:n1)
      pos.ylabel = pos.ylabel[1:min(n, m), ]
      symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE, 
              bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1, 
                                                                       length(pos.ylabel[, 1])))
      text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n, 
                                                                     m)], col = tl.col, cex = tl.cex, ...)
    }
    else {
      text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames, 
           srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5, 
                                                     0), c(0, 0)), col = tl.col, cex = tl.cex, offset = tl.offset, 
           ...)
      text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames, 
           col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset, 
           ...)
    }
  }
  title(title, ...)
  if (type == "full" && plotCI == "n" && !is.null(addgrid.col)) {
    rect(m1 - 0.5, n1 - 0.5, m2 + 0.5, n2 + 0.5, border = addgrid.col)
  }
  if (!is.null(addrect) && order == "hclust" && type == 
      "full") {
    corrRect.hclust(corr, k = addrect, method = hclust.method, 
                    col = rect.col, lwd = rect.lwd)
  }
  corrPos = data.frame(PosName, Pos, DAT)
  colnames(corrPos) = c("xName", "yName", "x", 
                        "y", "corr")
  # if (!is.null(p.mat)) {
  #   corrPos = cbind(corrPos, pNew)
  #   colnames(corrPos)[6] = c("p.value")
  # }
  corrPos = corrPos[order(corrPos[, 3], -corrPos[, 4]), ]
  rownames(corrPos) = NULL
  res = list(corr = corr, corrPos = corrPos, arg = list(type = type))
  invisible(res)
}



mantel_test <- function(chem_data,otu_table,nutrient_columns, method="spearman", ...){
  
  library(vegan) 
  
######## parameters ################
  #1. ~ chem_data - is a dataframe of envonmental variables / nutrients along with the description of the samples
  #with sample names as row names
  #2. ~ otu_table - an otu table matrix with samples as rows and OTUs or features as columns
  #3. ~ nutrient_columns - is a string or numeric vector of the numeric colums for the nutrients /
  #    environmental variables
  #4. ~ method the correlation method to use
  #5. ~ ... any arguments that can be passed to vegan::mantel
  
  common.ids <- intersect(rownames(chem_data),rownames(otu_table))
  #get the chemical data and otu tables for the sequenced samples
  chem_data <- chem_data[common.ids,]
  otu_table <- otu_table[common.ids,] 
  
  #calculate the bray-curtis distance between samples' OTUs 
  otu.dist <- vegdist(otu_table) # Bray-Curtis
  
  meta_table <- chem_data[,nutrient_columns]
  #get the environmenta variable names
  nutrients <- colnames(meta_table) 
  
  names(nutrients) <- nutrients
  
  
  stats_table <- map_dfr(.x = nutrients,.f = function(nutrient){
    
    
      nutrient_dist_meta_table <- subset(meta_table,select = nutrient)
      nutrient_dist <- vegdist(x = scale(nutrient_dist_meta_table),
                               method = "euclid")
      mantel_result <- mantel(otu.dist, nutrient_dist,  method=method, ...)
      
      data.frame(env_variable=nutrient,
                 r=mantel_result$statistic,
                 p_value=mantel_result$signif)
    
  }) %>% 
    mutate_if(.predicate = is.numeric,.funs = round, digits=3) 
  
  
  
  return(stats_table)
}


# Clean ancom microbe names
clean_name <- function(microbe_name){ microbe_name %>%
    str_replace_all(pattern = "^X|^X\\.", replacement = "") %>%
    str_replace_all(pattern = "\\._", replacement = ";_") %>% 
    str_replace_all(pattern = "\\.{1,}", replacement = " ") %>%  
    str_replace_all(pattern = "anthom", replacement = "Xanthom")}



group_box_plot <- function(df, x, y, groups, palette, x_lab, y_lab, box_color,
                           facet_plot= TRUE, text_size=6,
                           facet_by=NULL,  letters_vec=NULL){
  
  # Example :
  
  
  
  summary_table <- desc_statby(df, measure.var = y, grps = groups) %>% 
    pavian::zero_if_na() %>%
    droplevels()
  
  e <- try(summary_table %>% unite(col = "Treatment", 
                                   all_of(groups),
                                   sep = ":", remove = FALSE), 
           silent = TRUE)
  
  if(!inherits(e, "try-error")){
    
    summary_table <-  summary_table  %>%  unite(col = "Treatment", 
                                                all_of(groups),
                                                sep = ":", remove = FALSE)
    
  }
  
  if(!is.null(letters_vec)){
    
    labels <- letters_vec %>%
      as.data.frame() %>%
      set_names(., "labels") %>%
      rownames_to_column("Treatment")
    
    summary_table <-  summary_table %>% left_join(labels)
  }
  
  
  p <- ggplot(data = df, mapping = aes_string(x=x, y=y, fill=box_color))+ 
    geom_boxplot() + 
    scale_fill_manual(values = palette)
  
  
  # add text to plot
  if(!is.null(letters_vec)){
    
    
    toAdd <- if_else(condition = max(summary_table$range) <= 5,
                     true = min(summary_table$range),
                     false = median(summary_table$range) - min(summary_table$range))
    
    max_label_len <- str_length(summary_table$labels) %>% max()
    if(max_label_len  < 3){
      
      p <- p + geom_text(data=summary_table,
                         mapping = aes(y=max+toAdd, 
                                       label=labels,
                                       fontface = "bold"),
                         size = text_size, position=position_dodge(0.9))
      
    }else{
      
      text_size <- 4
      
      p <- p + ggrepel::geom_text_repel(data=summary_table,
                                        mapping = aes(y=max+toAdd, 
                                                      label=labels,
                                                      fontface = "bold"),
                                        size = text_size, position=position_dodge(0.9))
      
    }
  }
  
  
  
  if(facet_plot && (!is.null(facet_by)) ){
    p <- p +  facet_grid(facet_by, drop = TRUE, scale="free") }
  
  legend_title <-  box_color %>% str_replace_all('_',' ') %>% toTitleCase 
  p <- p +
    labs(x=x_lab, y=y_lab, fill=legend_title) +
    publication_format
  
  return(p)
  
}





# replace missing values using predictive modelling

replace_missing_values <- function(data, method=NULL, action=1, ...){
  
  # data  DATAFRAME - A dataframe with missing values to be imputed
  # action VECTOR -please see ?mice::complete
  # ... any arguement that can be passed to mice::mice
  
  # method can be one of the following - please see ?mice::mice for details
  
  #Built-in univariate imputation methods are:
  
  # pmm	  any	  Predictive mean matching
  # midastouch	any	  Weighted predictive mean matching
  # sample	  any	  Random sample from observed values
  # cart	  any	  Classification and regression trees
  # rf	  any	  Random forest imputations
  # mean  	numeric	  Unconditional mean imputation
  # norm	  numeric	  Bayesian linear regression
  # norm.nob	  numeric	Linear regression ignoring model error
  # norm.boot	 numeric	Linear regression using bootstrap
  # norm.predict	numeric	Linear regression, predicted values
  # lasso.norm	 numeric	Lasso linear regression
  # lasso.select.norm	 numeric	Lasso select + linear regression
  # quadratic	numeric	 Imputation of quadratic terms
  # ri	numeric	Random  indicator for nonignorable data
  # logreg	 binary	Logistic regression
  # logreg.boot	binary	Logistic regression with bootstrap
  # lasso.logreg	binary	Lasso logistic regression
  # lasso.select.logreg	binary	Lasso select + logistic regression
  # polr	ordered	Proportional odds model
  # polyreg	unordered	Polytomous logistic regression
  # lda	unordered	Linear discriminant analysis
  # 2l.norm	numeric	Level-1 normal heteroscedastic
  # 2l.lmer	numeric	Level-1 normal homoscedastic, lmer
  # 2l.pan	numeric	Level-1 normal homoscedastic, pan
  # 2l.bin	binary	Level-1 logistic, glmer
  # 2lonly.mean	numeric	Level-2 class mean
  # 2lonly.norm	numeric	Level-2 class normal
  # 2lonly.pmm	any	Level-2 class predictive mean matching
  
  
  # cart means classification and regression tres 
  mice_mod <- mice(data, method = method, ...)
  
  # action 1 means to return the first imputed datate
  data_mod <- complete(data = mice_mod, action = action)
  
  return(data_mod)
} 


# Get a list of a unique combination of elements in a vector

unique_vector_combination <- function(vec){
  
  cmp_list <- utils::combn(vec %>% as.character, 2) %>%
    as.data.frame() %>%
    map(function(col) if(col[1] != col[2]) col)
  
  cmp_list <- cmp_list[lengths(cmp_list) > 0] %>% S4Vectors::unique.Vector()
  return(cmp_list)
}


# A function to replot lefser plot for biomarker identification
# with names ofd the indentified biomarkers shortened
replot_lfserPlot<- function(res, theme2use, vec, legend_title="Treatment", isASV=TRUE){
  
  # - Vec is a vector of two factor groups where 
  # - the first group = 0 an dthe second group = 1
  
  if(isASV){
    
  taxon_levels <- c("Kingdom","Phylum", "Class", 
                    "Order", "Family", "Genus", "Species" )
  
  df <- res %>% separate(Names, into = taxon_levels, sep = "\\|")
  
  row_index <- rownames(df)
  names(row_index) <- row_index
  
  short_names <- future_map_chr(row_index, function(row){
    
    row_df <- df[row,taxon_levels] 
    
    last_unknown_column <- ( row_df %>%
                               grep(x=., pattern="__([0-9]+)?$"))[1]
    
    if(is_empty(last_unknown_column) || is.na(last_unknown_column) ){
      
      new_name <-  row_df[,"Species"]
      
    }else{
      new_name <- row_df[,last_unknown_column-1]
      
    }
    return(new_name)
  }) %>% as.data.frame() %>%
    set_names("short_name") %>%
    rownames_to_column("id")
  
  data2plot <- res %>% 
    rownames_to_column("id") %>%
    left_join(short_names)
  
  }else{
    
    data2plot <- res %>% dplyr::rename(short_name=Names)
  }
  

  vec <- sort(vec)
  colors <-  c("red", "forestgreen")
  names(colors)<- vec
  

  data2plot <- data2plot %>% 
    mutate(group=if_else(scores < 0, vec[1], vec[2]))
  
  plt <- ggplot(data2plot, aes(x=reorder(short_name, scores), scores)) + 
    geom_bar(stat = "identity", aes(fill = group)) +
    scale_fill_manual(values = colors) + 
    labs(x=NULL, y="LDA SCORE (log 10)", fill=legend_title) + 
    coord_flip() + 
    theme2use
  
  return(plt)
}


# # Format asv table to a format required to run lefse on the commadline
# 
# tabs2lefse <- function(asv_table, taxonomy_table, metadata){
#   
#   features <- taxonomy_table$Feature.ID
#   
#   
#   samples <- intersect(rownames(metadata), colnames(asv_table))
#   
#   asv_table <- asv_table[features, samples]
#   metadata <- metadata[samples,]
#   
#   tax_table <- taxonomy_table %>%  mutate(Taxon=gsub(sep, "|", x=Taxon))
#   
#   
# } 



run_lefse <- function(asv_table, metadata, taxonomy_table=NULL, 
                      group="Treatment", sep=";", 
                      theme2use=publication_format,isASV=TRUE,...){
  
  library(lefser)
  # Set see for reproducibility of results
  set.seed(151) 
  
  if(isASV){

    features <- taxonomy_table$Feature.ID
     
  }else{
    
    features <- rownames(asv_table)
}

  samples <- intersect(rownames(metadata), colnames(asv_table))
  metadata <- metadata[samples,]
  asv_table <- asv_table[features, samples]
  non_constant <- apply(X = asv_table,
                        MARGIN = 1, 
                        FUN = function(row) length(unique(row)) > 1 )
  asv_table <- asv_table[non_constant,,drop=FALSE] # Drop features with constant values across samples
  
  if(isASV){
  features <- intersect(rownames(taxonomy_table), rownames(asv_table))
  taxonomy_table <- taxonomy_table[features,]
  asv_table <- asv_table[features,]
  tax_table <- taxonomy_table %>%  mutate(Taxon=gsub(sep, "|", x=Taxon))
  # deduplicate taxonomy names
  rownames(asv_table) <- make.unique(tax_table$Taxon, sep = "_")
  }
  
  
  # Get a unique combination of a factor level and give them newly
  # created list names
  cmp_list <- unique_vector_combination(metadata[,group])
  
  names(cmp_list) <- cmp_list[lengths(cmp_list) > 0] %>%
    unique.Vector() %>% map_chr(function(x) paste(x[1], x[2], sep = ".vs."))
  
  # -------- Debug  -----------#
 # cmp_list  <- list("SA.vs.CMR"=c("SA", "CMR")) 

 #########################
  result <- future_map(cmp_list, .f = function(vec){
    
    sub_metadata <- metadata[which(metadata[,group] %in% vec),] %>%
      droplevels()
    vec <- sort(vec)
    sub_metadata[,group] <- factor(sub_metadata[,group], levels = vec)
    
    sub_asv_table <- asv_table[,rownames(sub_metadata)]
    
    non_constant <- apply(X = sub_asv_table,
                          MARGIN = 1, 
                          FUN = function(row) length(unique(row)) > 1 )
    
    sub_asv_table <- sub_asv_table[non_constant,,drop=FALSE] # Drop features with constant values across samples
    
    se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(sub_asv_table)),
                                colData=sub_metadata)
    
    res <- try(lefser(expr = se0, groupCol = group, ...))
    
    if(class(res) == "data.frame" && nrow(res) == 0){
      
      gg <- "No Significant diffrence detected"
      
    }else if(class(res) ==  "try-error"){
      
      gg <- "LDA error occurred: One variable appears to be constant within groups"
    
      }else{
      
    gg <- try(replot_lfserPlot(res, theme2use, vec, group, isASV))
    
    }
    
    return(list(result=res, plot=gg))
  
  })
  

  return(result)
}
