### VERSION 1 OF TAXA BAR PLOTS --> Specifically for facetting by treatment column in ginseng phyloseq format
### CREATED BY MICHAEL YOSHIDA, FOR THE TAI LAB, 2025

library(ggplot2)
library(phyloseq)

taxa_bar_plot <- function(
    physeq,
    fill = "Phylum",
    title = "",
    show_legend = TRUE,
    setting = "samples",
    min_abundance = 0.01,  # taxa with mean relative abundance below this threshold will be lumped as "Other"
    kingdom = "Fungi"      # specify whether the data are "Fungi" or "Bacteria"
) {
  # 1) Aggregate (glom) ASVs by the chosen taxonomic level
  physeq_glommed <- tax_glom(physeq, taxrank = fill)
  
  # 2) Convert to relative abundance (per sample)
  physeq_relA <- transform_sample_counts(physeq_glommed, function(x) x / sum(x))
  
  # 3) Count how many unique taxa we have at this level
  num_taxa <- length(get_taxa_unique(physeq_relA, taxonomic.rank = fill))
  message("Number of unique ", fill, " taxa: ", num_taxa)
  
  # 4) Convert the phyloseq object into a long-format data frame
  physeq_df <- psmelt(physeq_relA)
  
  # Make site character
  physeq_df$treatmentnum <- as.character(physeq_df$treatmentnum)
  
  physeq_df$treatmentnum <- factor(
    physeq_df$treatmentnum,
    levels = c(1, 2),
    labels = c("Control", "Garden")
  )
  
  
  # 5) Decide what goes on the x-axis (Samples or GroupID) and order accordingly
  x_var <- "Samples"
  if (setting == "samples") {
    # Order samples by seasonYear (assuming that column exists in physeq_df)
    physeq_df$Samples <- factor(
      physeq_df$Samples,
      levels = unique(physeq_df$Samples[order(physeq_df$seasonYear)])
    )
  } else if (setting == "grouped") {
    # Order group IDs by seasonYear
    physeq_df$GroupID <- factor(
      physeq_df$GroupID,
      levels = unique(physeq_df$GroupID[order(physeq_df$seasonYear)])
    )
    x_var <- "GroupID"
  } else {
    stop('Invalid setting. Choose either "samples" or "grouped".')
  }
  
  # 6) Lump taxa with mean relative abundance below the threshold into "Other"
  taxa_means <- physeq_df %>%
    dplyr::group_by(.data[[fill]]) %>%
    dplyr::summarize(mean_abundance = mean(Abundance)) %>%
    dplyr::arrange(desc(mean_abundance))
  
  taxa_to_keep <- taxa_means[[fill]][taxa_means$mean_abundance >= min_abundance]
  
  physeq_df <- physeq_df %>%
    dplyr::mutate(Legend = ifelse(.data[[fill]] %in% taxa_to_keep, .data[[fill]], "Other"))
  
  legend_order <- physeq_df %>%
    dplyr::group_by(Legend) %>%
    dplyr::summarize(total_abundance = sum(Abundance)) %>%
    dplyr::arrange(desc(total_abundance)) %>%
    dplyr::pull(Legend)
  
  if ("Other" %in% legend_order) {
    legend_order <- c(setdiff(legend_order, "Other"), "Other")
  }
  
  physeq_df$Legend <- factor(physeq_df$Legend, levels = legend_order)
  
  # after your mutate() and before you build p:
  #physeq_df$sample_type <- physeq_df$sample.type
  
  # 1) What are the names in your data?
  #print(names(physeq_df))
  
  # 2) Do you really have any non-NA sample types?
  #print(table(physeq_df$sample_type, useNA = "always"))
  
  # 3) Quick peek at the first few rows:
  #print(head(physeq_df[, c("Sample","site","sample.type","sample_type","GroupID")]))
  
  # 7) Create the bar plot with normal stacking (no reverse)
  p <- ggplot(physeq_df, aes_string(x = x_var, y = "Abundance", fill = "Legend")) +
    geom_bar(stat = "identity") +
    labs(
      x = x_var,
      y = "Relative Abundance",
      title = paste(title, fill, "Level")
    ) +
    theme_minimal() +
    facet_wrap(~ treatmentnum, scales = "free_x") +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # 8) Define hardcoded color palettes for fungi and bacteria
  
  # Fungal palettes
  class_colors <- c(
    "Agaricomycetes"             = "#4E79A7",
    "Dothideomycetes"            = "#F28E2B",
    "Eurotiomycetes"             = "#E15759",
    "Lecanoromycetes"            = "#76B7B2",
    "Leotiomycetes"              = "#59A14F",
    "Other"                      = "#EDC948",
    "Pezizomycetes"              = "#B07AA1",
    "Rhizophlyctidomycetes"      = "#FF9DA7",
    "Saccharomycetes"            = "#9C755F",
    "Sordariomycetes"            = "#BAB0AC",
    "Unidentified_Ascomycota"    = "#C7E9B4",
    "Unidentified_Basidiomycota" = "#FCAE91",
    "Unidentified_Fungi"         = "#FEE0D2"
  )
  
  order_colors <- c(
    "Agaricales"                   = "#4E79A7",
    "Capnodiales"                  = "#F28E2B",
    "Chaetothyriales"              = "#FF9D9A",
    "Eurotiales"                   = "#E15759",
    "Helotiales"                   = "#59A14F",
    "Hysteriales"                  = "#76B7B2",
    "Hypocreales"                  = "#BAB0AC",
    "Other"                        = "#EDC948",
    "Pezizales"                    = "#B07AA1",
    "Pleosporales"                 = "#9C755F",
    "Rhizophlyctidales"            = "#FFBE7D",
    "Saccharomycetales"            = "#D37295",
    "Spizellomycetales"            = "#8CD17D",
    "Unidentified_Agaricomycetes"  = "#C7E9B4",
    "Unidentified_Ascomycota"      = "#FCAE91",
    "Unidentified_Basidiomycota"   = "#FEE0D2",
    "Unidentified_Dothideomycetes" = "#AF7AA1",
    "Unidentified_Fungi"           = "#D4A6C8",
    "Unidentified_Sordariomycetes" = "#F1CE63"
  )
  
  family_colors <- c(
    "Nectriaceae"                           = "#4E79A7",
    "Aspergillaceae"                        = "#F28E2B",
    "Pleosporaceae"                         = "#E15759",
    "Hysteriaceae"                          = "#76B7B2",
    "Clavicipitaceae"                       = "#59A14F",
    "Dermateaceae"                          = "#EDC948",
    "Mucoraceae"                            = "#B07AA1",
    "Cadophoraceae"                         = "#FF9DA7",
    "Lophiostomataceae"                     = "#9C755F",
    "Cladosporiaceae"                       = "#BAB0AC",
    "Herpotrichiellaceae"                   = "#FFBE7D",
    "Pezizaceae"                           = "#D37295",
    "Unidentified_Agaricales"               = "#8CD17D",
    "Rhizophlyctidaceae"                    = "#CAB2D6",
    "Unidentified_Pleosporales"             = "#A6CEE3",
    "Saccharomycetales_fam_Incertae_sedis"  = "#1F78B4",
    "Spizellomycetaceae"                    = "#33A02C",
    "Rhizophlyctidales_fam_Incertae_sedis"    = "#FB9A99",
    "Other"                                 = "#D37295",
    "Unidentified_Fungi"                    = "#8CD17D",
    "Unidentified_Agaricomycetes"           = "#C7E9B4",
    "Unidentified_Ascomycota"               = "#FCAE91",
    "Unidentified_Basidiomycota"            = "#FEE0D2",
    "Unidentified_Helotiales"               = "#AF7AA1",
    "Unidentified_Hypocreales"              = "#D4A6C8",
    "Unidentified_Sordariomycetes"          = "#F1CE63",
    "Unidentified_Dothideomycetes"          = "#CAB2D6"
  )
  
  # Bacterial palettes 
  bacteria_class_colors <- c(
    "Alphaproteobacteria" = "#66c2a5",
    "Betaproteobacteria"  = "#fc8d62",
    "Gammaproteobacteria" = "#8da0cb",
    "Bacilli"             = "#e78ac3",
    "Clostridia"          = "#a6d854",
    "Actinobacteria"      = "#ffd92f",
    "Thermoleophilia"     = "#1f78b4",
    "Bacteroidia"         = "#33a02c",
    "Chloroflexia"        = "#6a3d9a",
    "Gemmatimonadetes"    = "#ff7f00",
    "Nitrososphaeria"     = "#b15928",
    "Anaerolineae"        = "#cab2d6",
    "Ktedonobacteria"     = "#fb9a99",
    "Planctomycetes"      = "#fdbf6f",
    "Polyangia"           = "#b2df8a",
    "Verrucomicrobiae"    = "#a6cee3",
    "Vicinamibacteria"    = "#bc80bd",
    "Blastocatellia"      = "#ffff99",
    "Acidobacteriae"      = "#8c8c8c",  
    "KD4-96"              = "#bc8f8f",  
    "Other"               = "#e5c494"
  )
  
  
  
  
  bacteria_order_colors <- c(
    "Enterobacterales"       = "#66c2a5",
    "Pseudomonadales"        = "#fc8d62",
    "Lactobacillales"        = "#8da0cb",
    "Bacillales"             = "#e78ac3",
    "Bacteroidales"          = "#a6d854",
    "Micrococcales"          = "#ffd92f",
    "Alicyclobacillales"     = "#1f78b4",
    "Paenibacillales"        = "#33a02c",
    "Rhizobiales"            = "#6a3d9a",
    "Sphingomonadales"       = "#ff7f00",
    "Burkholderiales"        = "#b15928",
    "Solirubrobacterales"    = "#cab2d6",
    "Frankiales"             = "#fb9a99",
    "Gemmatimonadales"       = "#fdbf6f",
    "Gaiellales"             = "#b2df8a",
    "Pseudonocardiales"      = "#a6cee3",
    "Propionibacteriales"    = "#bc80bd",
    "Nitrososphaerales"      = "#ffff99",
    "Clostridiales"          = "#e31a1c",
    "Chloroflexales"         = "#fb8072",
    "Chitinophagales"        = "#d9d9d9",
    "Streptomycetales"       = "#fdae61",
    "Ktedonobacterales"      = "#abd9e9",
    "Micromonosporales"      = "#b3b3b3",
    "Flavobacteriales"       = "#8dd3c7",
    "Pirellulales"           = "#e6ab02",
    "Vicinamibacterales"     = "#bcbd22",
    "Gammaproteobacteria Incertae Sedis" = "#999999",
    "Anaerolineales"         = "#e7298a",
    "Streptosporangiales"    = "#66c2a5",
    "Chthoniobacterales"     = "#f781bf",
    "Corynebacteriales"      = "#a65628",
    "Thermomicrobiales"      = "#1b9e77",
    "Other"                  = "#e5c494"
  )
  
  
  bacteria_family_colors <- c(
    "Enterobacteriaceae"        = "#393b79",
    "Pseudomonadaceae"          = "#5254a3",
    "Lactobacillaceae"          = "#6b6ecf",
    "Bacillaceae"               = "#9c9ede",
    "Bacteroidaceae"            = "#637939",
    "Other"                     = "#8ca252",
    "Chitinophagaceae"          = "#b5cf6b",
    "Sphingomonadaceae"         = "#cedb9c",
    "Gemmatimonadaceae"         = "#8c6d31",
    "Xanthobacteraceae"         = "#bd9e39",
    "Pirellulaceae"             = "#e7ba52",
    "Pyrinomonadaceae"          = "#e7cb94",
    "Nitrosomonadaceae"         = "#843c39",
    "Alicyclobacillaceae"       = "#ad494a",
    "Nitrososphaeraceae"        = "#d6616b",
    "Vicinamibacteraceae"       = "#e7969c",
    "Micrococcaceae"            = "#7b4173",
    "Paenibacillaceae"          = "#a55194",
    "Solirubrobacteraceae"      = "#ce6dbd",
    "Nocardioidaceae"           = "#de9ed6",
    "KD4-96"                    = "#3182bd",
    "67-14"                     = "#6baed6",
    "Micromonosporaceae"        = "#9ecae1",
    "Beijerinckiaceae"          = "#c6dbef",
    "Blastocatellaceae"         = "#e6550d",
    "Roseiflexaceae"            = "#fd8d3c",
    "Pseudonocardiaceae"        = "#fdae6b",
    "Planococcaceae"            = "#fdd0a2",
    "Gaiellaceae"               = "#31a354",
    "Anaerolineaceae"           = "#74c476",
    "Comamonadaceae"            = "#a1d99b",
    "C0119"                     = "#c7e9c0",
    "Bryobacteraceae"           = "#756bb1",
    "Rhizobiales_Incertae_Sedis" = "#9e9ac8"
  )
  
  
  # Bacterial phylum palette for manual scaling when fill is "Phylum"
  bacteria_phylum_colors <- c(
    "Actinobacteriota"   = "#8dd3c7",
    "Firmicutes"         = "#fc8d62",
    "Proteobacteria"     = "#bebada",
    "Chloroflexi"        = "#fb8072",
    "Bacteroidota"       = "#80b1d3",
    "Acidobacteriota"    = "#fdb462",
    "Myxococcota"        = "#bc80bd",
    "Planctomycetota"    = "#ccebc5",
    "Crenarchaeota"      = "#fccde5",
    "Verrucomicrobiota"  = "#d9d9d9",
    "Gemmatimonadota"    = "#ffed6f",
    "Other"              = "#b3b3b3"
  )
  
  # 9) Apply the conditional color scale based on kingdom and fill level
  if (kingdom == "Fungi") {
    if (fill == "Class") {
      p <- p + scale_fill_manual(values = class_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(class_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else if (fill == "Order") {
      p <- p + scale_fill_manual(values = order_colors)
    } else if (fill == "Family") {
      p <- p + scale_fill_manual(values = family_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(family_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else {
      p <- p + scale_fill_brewer(palette = "Set3")
    }
  } else if (kingdom == "Bacteria") {
    if (fill == "Phylum") {
      p <- p + scale_fill_manual(values = bacteria_phylum_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(bacteria_phylum_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else if (fill == "Class") {
      p <- p + scale_fill_manual(values = bacteria_class_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(bacteria_class_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else if (fill == "Order") {
      p <- p + scale_fill_manual(values = bacteria_order_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(bacteria_order_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else if (fill == "Family") {
      p <- p + scale_fill_manual(values = bacteria_family_colors)
      unmapped_levels <- setdiff(levels(physeq_df$Legend), names(bacteria_family_colors))
      if (length(unmapped_levels) > 0) {
        message("Unmapped Legend levels: ", paste(unmapped_levels, collapse = ", "))
      }
    } else {
      p <- p + scale_fill_brewer(palette = "Set3")
    }
  } else {
    p <- p + scale_fill_brewer(palette = "Set3")
  }
  
  return(p)
}


