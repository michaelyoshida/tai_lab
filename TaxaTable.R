### SCRIPT TO CREATE TABLE-STYLE RELATIVE ABUNDANCE BY TAXA
### CREATED BY MICHAEL YOSHIDA FOR THE TAI LAB, 2025.

library(phyloseq)
library(dplyr)
library(tidyr)

taxa_abundance_table <- function(
    physeq,
    fill          = "Phylum",    # taxonomic rank used to aggregate ASVs
    setting       = "samples",   # "samples" or "grouped"
    min_abundance = 0.01,        # mean rel. abundance threshold for "Other"
    kingdom       = "Fungi",     # kept for interface parity; not used here
    save_csv      = NULL,        # e.g. "relabund_table.csv" (NULL = don't save)
    wide          = FALSE        # TRUE = samples/group x taxa wide matrix
) {
  
  # aggregate ASVs at the chosen taxonomic rank
  physeq_glommed <- tax_glom(physeq, taxrank = fill)
  
  # convert sample counts to relative abundance 
  physeq_relA <- transform_sample_counts(physeq_glommed, function(x) x / sum(x))
  
  # Melt (re-format) to long data frame
  df <- psmelt(physeq_relA)
  
  # Keep treatment labels if present (if present was debug)
  df$treatmentnum <- factor(
    as.character(df$treatmentnum),
    levels = c("1", "2"),
    labels = c("Control", "Garden")
  )
  
  sample_col <- "Sample"
  
  # use individual samples (not grouped) if specified
  x_var <- sample_col
  
  if (setting == "samples") {
    
    # Order individual samples chronologically by year, season, and sample ID
    ord <- order(
      df$year,
      df$season,
      df[[sample_col]]
    )
    
    df[[sample_col]] <- factor(
      df[[sample_col]],
      levels = unique(df[[sample_col]][ord])
    )
    
    #
  } else if (setting == "grouped") {
    
    x_var <- "GroupID"
    
    # Order grouped observations chronologically by year, season, and GroupID
    ord <- order(
      df$year,
      df$season,
      df$GroupID
    )
    
    df$GroupID <- factor(
      df$GroupID,
      levels = unique(df$GroupID[ord])
    )
    
  } else {
    stop('Invalid setting. Choose either "samples" or "grouped".')
  }
  
  # calculate mean relative abundance for each taxa across all sampling times and filter based on threshold argument, filtered out go in "other"
  taxa_means <- df %>%
    dplyr::group_by(.data[[fill]]) %>%
    dplyr::summarize(mean_abundance = mean(Abundance), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abundance))
  
  taxa_to_keep <- taxa_means[[fill]][taxa_means$mean_abundance >= min_abundance]
  
  df <- df %>%
    dplyr::mutate(Taxon = ifelse(.data[[fill]] %in% taxa_to_keep,
                                 .data[[fill]], "Other"))
  
  # order taxon levels by total abundance (keep "Other" last)
  legend_order <- df %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarize(total_abundance = sum(Abundance), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_abundance)) %>%
    dplyr::pull(Taxon)
  
  # place other catgoery at the end
  if ("Other" %in% legend_order) {
    legend_order <- c(setdiff(legend_order, "Other"), "Other")
  }
  df$Taxon <- factor(df$Taxon, levels = legend_order)
  
  # build taxa table to export 
  keep_cols <- intersect(
    c("GroupID", "SampleID", "season", "year",
      "site", "treatmentnum", "SampleType", "sample.type"),
    names(df)
  )
  
  # create one row per sample/group and taxon combination
  out_long <- df %>%
    dplyr::transmute(
      !!x_var       := .data[[x_var]],   # this is the only "ID" column
      Taxon         = Taxon,
      Abundance     = Abundance,
      Abundance_pct = Abundance * 100
    ) %>%
    dplyr::bind_cols(df[, setdiff(keep_cols, x_var), drop = FALSE]) %>%
    dplyr::relocate(any_of(c(x_var, "Taxon", "Abundance", "Abundance_pct")))
  
  # export to CSV 
  if (!is.null(save_csv)) {
    write.csv(out, file = save_csv, row.names = FALSE)
    message("Saved relative abundance table to: ", normalizePath(save_csv))
  }
  
  #return completed abundance table
  return(out)
}


group <- average_merge_phyloseq(ginseng_rare)

test <- taxa_abundance_table(
  group,
  fill = "Order",
  setting = "grouped",
  min_abundance = 0.01,
  kingdom = "Fungi",
  save_csv = "ginseng_order_grouped_relabund.csv",
  wide = FALSE
)

test

bacteria <- average_merge_phyloseq(bacteria_rare)

test <- taxa_abundance_table(
  bacteria,
  fill = "Class",
  setting = "grouped",
  min_abundance = 0.01,
  kingdom = "Fungi",
  save_csv = "bacteria_class_grouped_relabund.csv",
  wide = FALSE
)

test



