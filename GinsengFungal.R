### THIS SCRIPT IS THE MAIN ONE FOR THE GINSENG FUNGAL PROJECT.
### IT CREATES ALL PHYLOSEQ OBJECTS NEEDED FOR DOWNSTREAM ANALYSIS

###LIBRARIES SECTION

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(vegan)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(ggrepel)

### IMPORT AND FORMAT PHYLOSEQ OBJECT

setwd("/Users/michaelyoshida/Downloads/TAI LAB WORK")

fungi <- qza_to_phyloseq(
  features="264196shortginsengfungal-table.qza",
  tree="264196shortginsengfungal-rooted-tree.qza",
  taxonomy="264196shortginsengfungal-taxonomyv2.qza",
  metadata = "ginsengfungalmetadata.tsv"
)

sample_data(fungi)

#getting rid of the positive controls
samplesToRemove <- c("positive1","positive2","positive3")
sampleNames <- sample_names(fungi)
keep <- !(sampleNames %in% samplesToRemove)
fungi <- subset_samples(fungi, keep)

#subset out non-fungi from taxonomy
fungi_only <- subset_taxa(fungi, Kingdom!= "Eukaryota_kgd_Incertae_sedis" & Kingdom!="Eukaryota_unidentified" 
                          & Kingdom!="Unassigned" & Kingdom!="Alveolata" & Kingdom!="Alveolata_unidentified" 
                          & Kingdom!="Viridiplantae" & Kingdom!= "Viridiplantae_unidentified" 
                          & Kingdom!= "Rhizaria" & Kingdom!= "Rhizaria_unidentified" 
                          & Kingdom!= "Stramenopila" & Kingdom!= "Stramenopila_unidentified" 
                          & Kingdom!= "Ichthyosporia" & Kingdom!= "Metazoa_unidentified") 

fungi_only_samples <- as.data.frame(sample_data(fungi_only))
fungi_only_samples$seasonYear <- paste0(fungi_only_samples$year, "-", fungi_only_samples$season)
fungi_only_samples$seasonYear <- factor(fungi_only_samples$seasonYear, 
                                        levels = unique(fungi_only_samples$seasonYear[order(fungi_only_samples$year, 
                                                                                            factor(fungi_only_samples$season,                                                                       levels = c("spring", "summera", "summerb","fall")))]))
#updating phyloseq object to have column with year and season 
sample_data(fungi_only) <- fungi_only_samples

nsamples(fungi_only) # --> should be 244 before read count filtering

fungi_only <- prune_samples(sample_sums(fungi_only) > 1000, fungi_only)

nsamples(fungi_only) # --> 233 remain, 11 removed


# 1) Remember the original taxa order
orig_taxa <- taxa_names(fungi_only)

# 2) Extract taxonomy table as a data.frame
tax.clean <- as.data.frame(
  tax_table(fungi_only),
  stringsAsFactors = FALSE
)

# 3) Ensure every column is character
tax.clean[] <- lapply(tax.clean, as.character)

# 4) Replace NA and "__" placeholders with blanks
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean == "__"] <- ""

# 5) Strip out unwanted “uncultured…”, “unidentified…”, “metagenome…” bits
bad_patterns <- c(
  "(?i)\\buncultured[_A-Za-z0-9]*",
  "(?i)\\bunidentified[_A-Za-z0-9]*",
  "(?i)\\bmetagenome[_A-Za-z0-9]*"
)
tax.clean[] <- lapply(tax.clean, function(col) {
  for (pat in bad_patterns) {
    col <- gsub(pat, "", col, perl = TRUE)
  }
  # Trim any leading/trailing underscores or spaces
  gsub("^[_ ]+|[_ ]+$", "", col)
})

# 6) Fill down missing ranks
for (i in seq_len(nrow(tax.clean))) {
  if (tax.clean[i, "Phylum"] == "") {
    new <- paste0("Unidentified_", tax.clean[i, "Kingdom"])
    tax.clean[i, c("Phylum","Class","Order","Family","Genus","Species")] <- new
    
  } else if (tax.clean[i, "Class"] == "") {
    new <- paste0("Unidentified_", tax.clean[i, "Phylum"])
    tax.clean[i, c("Class","Order","Family","Genus","Species")] <- new
    
  } else if (tax.clean[i, "Order"] == "") {
    new <- paste0("Unidentified_", tax.clean[i, "Class"])
    tax.clean[i, c("Order","Family","Genus","Species")] <- new
    
  } else if (tax.clean[i, "Family"] == "") {
    new <- paste0("Unidentified_", tax.clean[i, "Order"])
    tax.clean[i, c("Family","Genus","Species")] <- new
    
  } else if (tax.clean[i, "Genus"] == "") {
    new <- paste0("Unidentified_", tax.clean[i, "Family"])
    tax.clean[i, c("Genus","Species")] <- new
    
  } else if (tax.clean[i, "Species"] == "") {
    tax.clean[i, "Species"] <- paste0("Unidentified_", tax.clean[i, "Genus"])
  }
}

# 7) Restore and reorder row names to match the original phyloseq object
rownames(tax.clean) <- orig_taxa
tax.clean <- tax.clean[orig_taxa, ]

# 8) Assign the cleaned taxonomy back into your phyloseq object
tax_table(fungi_only) <- as.matrix(tax.clean)


### MAKE SUBSET OBJECTS + RAREFY OBJECTS

ginseng_be_garden_s1 <- subset_samples(fungi_only, site == 1 & year %in% c("2018", "2022") & season == "fall" & treatment == "garden")
ginseng_be_garden_s2 <- subset_samples(fungi_only, site == 2 & year %in% c("2018", "2022") & season == "fall" & treatment == "garden")
ginseng_be_garden_s3 <- subset_samples(fungi_only, site == 3 & year %in% c("2018", "2022") & season == "fall" & treatment == "garden")

ginseng_be_control_s1 <- subset_samples(fungi_only, site == 1 & year %in% c("2018", "2022") & season == "fall" & treatment == "control")
ginseng_be_control_s2 <- subset_samples(fungi_only, site == 2 & year %in% c("2018", "2022") & season == "fall" & treatment == "control")
ginseng_be_control_s3 <- subset_samples(fungi_only, site == 3 & year %in% c("2018", "2022") & season == "fall" & treatment == "control")

ginseng_site1_garden <- subset_samples(fungi_only, site == 1 & treatment == "garden")
ginseng_site2_garden <- subset_samples(fungi_only, site == 2 & treatment == "garden")
ginseng_site3_garden <- subset_samples(fungi_only, site == 3 & treatment == "garden")

ginseng_site1_control <- subset_samples(fungi_only, site == 1 & treatment == "control")
ginseng_site2_control <- subset_samples(fungi_only, site == 2 & treatment == "control")
ginseng_site3_control <- subset_samples(fungi_only, site == 3 & treatment == "control")


ginseng_rare = rarefy_even_depth(fungi_only)

average_merge_phyloseq <- function(phyloseq_obj) {
  # firstly create Group ID column for which merge_samples will group on
  temp <- as.data.frame(sample_data(phyloseq_obj))
  temp$GroupID <- sub("([0-9])$", "", rownames(temp))
  sample_data(phyloseq_obj) <- sample_data(temp)
  
  merged_phyloseq <- merge_samples(phyloseq_obj, group = "GroupID", fun = mean)
  
  # transform to relative abundance, unfortunately destroys Group ID
  transformed_phyloseq <- transform_sample_counts(merged_phyloseq, function(x) x / sum(x))
  
  # Reassign the new sample names and remake 'GroupID' column in sample_data
  new_sample_data <- as.data.frame(sample_data(transformed_phyloseq))
  new_sample_data$GroupID <- rownames(new_sample_data)
  sample_data(transformed_phyloseq) <- sample_data(new_sample_data)
  
  return(transformed_phyloseq)
}

ginseng_rare_site1 = subset_samples(ginseng_rare, site == 1)
ginseng_rare_site2 = subset_samples(ginseng_rare, site == 2)
ginseng_rare_site3 = subset_samples(ginseng_rare, site == 3)

ginseng_rare_garden = subset_samples(ginseng_rare, treatment == "garden")
ginseng_rare_site1_garden = subset_samples(ginseng_rare_garden, site == 1)
ginseng_rare_site2_garden = subset_samples(ginseng_rare_garden, site == 2)
ginseng_rare_site3_garden = subset_samples(ginseng_rare_garden, site == 3)

ginseng_rare_control = subset_samples(ginseng_rare, treatment == "control")
ginseng_rare_site1_control = subset_samples(ginseng_rare_control, site == 1)
ginseng_rare_site2_control = subset_samples(ginseng_rare_control, site == 2)
ginseng_rare_site3_control = subset_samples(ginseng_rare_control, site == 3)

ginseng_rare_beginning_end <- subset_samples(ginseng_rare, 
                                             year %in% c("2018", "2022") & season == "fall" )

ginseng_rare_beginning_end_garden <- subset_samples(ginseng_rare_beginning_end, treatment == "garden")
ginseng_rare_beginning_end_control <- subset_samples(ginseng_rare_beginning_end, treatment == "control")

ginseng_rare_be_garden_site1 = subset_samples(ginseng_rare_beginning_end_garden, site == 1)
ginseng_rare_be_garden_site2 = subset_samples(ginseng_rare_beginning_end_garden, site == 2)
ginseng_rare_be_garden_site3 = subset_samples(ginseng_rare_beginning_end_garden, site == 3)

ginseng_rare_be_control_site1 = subset_samples(ginseng_rare_beginning_end_control, site == 1)
ginseng_rare_be_control_site2 = subset_samples(ginseng_rare_beginning_end_control, site == 2)
ginseng_rare_be_control_site3 = subset_samples(ginseng_rare_beginning_end_control, site == 3)


### ALPHA DIVERSITY SECTION


plot_shannon_alphadiv <- function(phyloseq_object, column, plot_title) {
  alpha_div <- estimate_richness(phyloseq_object, measures = "Shannon")
  
  alpha_div[[column]] <- sample_data(phyloseq_object)[[column]]
  alpha_div$treatment <- sample_data(phyloseq_object)$treatment
  
  custom_colors <- c("control" = "sienna1", "garden" = "cornflowerblue")
  
  #boxplot
  p <- ggplot(alpha_div, aes_string(x = column, y = "Shannon", color = "treatment")) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_color_manual(values = custom_colors) +
    ggtitle(plot_title) +
    theme_classic(base_size = 16) +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted")
    ) +
    labs(
      x = column,
      y = "Shannon Diversity Index",
      color = "Treatment"
    )
  
  return(p)
}
plot_shannon_alphadiv(ginseng_rare, "seasonYear", "Shannon Alpha Diversity by Year")


#### BETA DIVERSITY SECTION


# one beta diversity plot per sample site and all together (for 2018 and 2022 end)
# helper function I made to try to fix repeitive legend issues
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

beta_pcoa_plot <- function(physeq_obj, title, show_legend = TRUE, shape = "treatment") {
  # Calculate beta diversity and perform PCoA
  betaDiversity <- distance(physeq_obj, method = "bray")
  ordination <- ordinate(physeq_obj, method = "PCoA", distance = betaDiversity)
  eigenvalues <- ordination$values$Relative_eig * 100
  x_variation <- round(eigenvalues[1], 2)  # First axis
  y_variation <- round(eigenvalues[2], 2)  # Second axis
  
  sample_data(physeq_obj)$year <- factor(sample_data(physeq_obj)$year)
  sample_data(physeq_obj)[[shape]]        <- factor(sample_data(physeq_obj)[[shape]])
  
  # Automatically generate shape mapping based on the unique levels of the chosen column
  shape_levels <- unique(as.character(sample_data(physeq_obj)[[shape]]))
  auto_shape_mapping <- setNames(seq(0, length(shape_levels) - 1), shape_levels)
  
  p <- plot_ordination(physeq_obj, ordination, color = "year", shape = shape) + 
    geom_point(size = 3, alpha = 0.7) + 
    ggtitle(title) + 
    scale_color_manual(
      values = c(
        "2018" = "#E41A1C",  # Red
        "2019" = "#FF7F00",  # Orange
        "2020" = "#4DAF4A",  # Green
        "2021" = "#984EA3",  # Purple
        "2022" = "#377EB8"   # Blue
      ),
      breaks = c("2018", "2019", "2020", "2021", "2022"),
      drop = FALSE
    ) + 
    #scale_shape_manual(values = auto_shape_mapping) + 
    theme_classic(base_size = 16) + 
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12)
    ) + 
    labs(
      x = paste0("PCoA Axis 1 (", x_variation, "%)"),
      y = paste0("PCoA Axis 2 (", y_variation, "%)"),
      color = "Year",
      shape = shape
    )
  
  #uncomment for sample tags
  #return(p + geom_text_repel(aes(label = sample_names(physeq_obj)), size = 3))
  return(p)
}


# Generate PCoA plots for each site (without legend)
p1 <- beta_pcoa_plot(ginseng_rare_site1_garden, "Site 1", show_legend = FALSE)
p2 <- beta_pcoa_plot(ginseng_rare_site2_garden, "Site 2", show_legend = FALSE)
p3 <- beta_pcoa_plot(ginseng_rare_site3_garden, "Site 3", show_legend = FALSE)
p4 <- beta_pcoa_plot(ginseng_rare_garden, "All Sites", show_legend = TRUE)  

# Extract the legend separately before removing it from p4
legend <- get_legend(p4)

# delete P4's legend now that "legend" contains it, to prevent prior duplication issues
p4 <- beta_pcoa_plot(ginseng_rare_garden, "All Sites", show_legend = FALSE)

# ggplot version of par
combined_plot <- ((p1 + p2) / (p3 + p4))  

# Combine the figure with the extracted legend
final_plot <- wrap_plots(combined_plot, legend, ncol = 1, heights = c(4, 0.3))

print(final_plot)

p5 <- beta_pcoa_plot(ginseng_rare_beginning_end_garden, "PCoA - Beta Diversity of Garden Samples (Bray-Curtis), 2018 to 2022", show_legend = TRUE, shape = "site") 
print(p5)
#save as PDF 
#ggsave("PCoA_Beta_Diversity_Sites_With_Legend.pdf", plot = final_plot, width = 12, height = 14, dpi = 300)

### TAXA BAR PLOT SECTION
taxa_bar_plot <- function(
    physeq,
    fill = "Phylum",
    title = "",
    show_legend = TRUE,
    setting = "samples",
    min_abundance = 0.01   # taxa with mean relative abundance below this threshold will be lumped as "Other"
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
  
  # 5) Decide what goes on the x-axis (Samples or GroupID) and order accordingly
  x_var <- "Sample"
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
  # Calculate the mean abundance (across samples) for each taxon at the chosen level
  taxa_means <- physeq_df %>%
    dplyr::group_by(.data[[fill]]) %>%
    dplyr::summarize(mean_abundance = mean(Abundance)) %>%
    dplyr::arrange(desc(mean_abundance))
  
  # Identify taxa that meet or exceed the threshold
  taxa_to_keep <- taxa_means[[fill]][taxa_means$mean_abundance >= min_abundance]
  
  # Create a new column that assigns taxa not in the list to "Other"
  physeq_df <- physeq_df %>%
    dplyr::mutate(Legend = ifelse(.data[[fill]] %in% taxa_to_keep, .data[[fill]], "Other"))
  
  # 7) Make the bar plot with the new "Legend" as fill
  p <- ggplot(physeq_df, aes_string(x = x_var, y = "Abundance", fill = "Legend")) +
    geom_bar(stat = "identity") +
    labs(
      x = x_var,
      y = "Relative Abundance",
      title = paste(title, fill, "Level")
    ) +
    theme_minimal() +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_brewer(palette = "Set3")
  
  return(p)
}



t1 <- taxa_bar_plot(ginseng_rare_be_garden_site1, title="Site 1 Garden", show_legend = FALSE)
t2 <- taxa_bar_plot(ginseng_rare_be_garden_site2, title="Site 2 Garden", show_legend = FALSE)
t3 <- taxa_bar_plot(ginseng_rare_be_garden_site3, title="Site 3 Garden", show_legend = FALSE)
t4 <- taxa_bar_plot(ginseng_rare, title = "Al Sites")

taxa_legend <- get_legend(t4)

t4 <- taxa_bar_plot(ginseng_rare, title = "All Sites", show_legend = FALSE)

combined_taxa_plot <- ((t1 + t2) / (t3 + t4))

final_taxa_plot <- wrap_plots(combined_taxa_plot, taxa_legend, ncol = 1, heights = c(4, 0.3))
print(final_taxa_plot)


t5 <- taxa_bar_plot(ginseng_rare_beginning_end_garden, title="2018 and 2022, Garden", show_legend = FALSE)
t6 <- taxa_bar_plot(ginseng_rare_beginning_end_control, title = "2018 and 2022, Control")
compare_legend <- get_legend(t6)

t6 <- taxa_bar_plot(ginseng_rare_beginning_end_control, title = "2018 and 2022, Control", show_legend = FALSE)
direct_comparison <- (t5 + t6)

compare <- wrap_plots(direct_comparison, compare_legend, ncol = 1, heights = c(4, 0.3))
print(compare)

### ALDEX PLOT SECTION

aldex_volcano_plot <- function(physeq,
                               sample_col = "year",
                               glom_rank = NULL,       # new parameter: taxonomic level to agglomerate by
                               sig_cutoff = 0.05,
                               mc.samples = 128,
                               denom = "iqlr",
                               verbose = TRUE,
                               title = "Volcano Plot of Differential Abundance",
                               max_y = 3) {
  # If a taxonomic level is provided, agglomerate taxa at that level
  if (!is.null(glom_rank)) {
    physeq <- tax_glom(physeq, taxrank = glom_rank)
  }
  
  # 1. Convert the OTU table to a data frame
  samplesALDeX <- as.data.frame(otu_table(physeq))
  
  # 2. Extract the sample grouping vector using the provided sample_data column
  sample_type_vec <- sample_data(physeq)[[sample_col]]
  
  # 3. Perform the CLR transformation using ALDEx2
  x <- aldex.clr(samplesALDeX, sample_type_vec,
                 mc.samples = mc.samples, denom = denom, verbose = verbose)
  
  # 4. Run differential abundance tests: Welch’s t-test and compute effect sizes
  x.tt <- aldex.ttest(x, paired.test = FALSE)
  x.effect <- aldex.effect(x, include.sample.summary = TRUE, verbose = verbose)
  
  # 5. Merge the test results and effect sizes into one data frame
  x.all <- data.frame(x.tt, x.effect)
  x.all$OTU_ID <- rownames(x.all)
  
  # 6. Retrieve and merge taxonomy
  taxonomy_df <- as.data.frame(tax_table(physeq))
  taxonomy_df$OTU_ID <- rownames(taxonomy_df)
  
  # Create a 'TaxonName' column (adjust to your preferred taxonomy ranks)
  taxonomy_df$TaxonName <- paste0(
    ifelse(is.na(taxonomy_df$Family), "", paste0(taxonomy_df$Family, "; ")),
    ifelse(is.na(taxonomy_df$Genus),  "", paste0(taxonomy_df$Genus,  "; ")),
    ifelse(is.na(taxonomy_df$Species), "", taxonomy_df$Species)
  )
  
  x.all_annot <- merge(x.all, taxonomy_df, by = "OTU_ID", all.x = TRUE)
  
  # For any row that doesn't have taxonomy info, fall back to OTU/ASV ID
  no_tax_idx <- which(is.na(x.all_annot$TaxonName) | x.all_annot$TaxonName == "")
  if (length(no_tax_idx) > 0) {
    x.all_annot$TaxonName[no_tax_idx] <- x.all_annot$OTU_ID[no_tax_idx]
  }
  
  # 7. Flag taxa as "Significant" or "Not Significant" based on the cutoff
  x.all_annot$sig <- ifelse(x.all_annot$wi.eBH < sig_cutoff,
                            "Significant", "Not Significant")
  
  # 8. Create a more publication-ready volcano plot
  volcano_plot <- ggplot(
    x.all_annot,
    aes(x = diff.btw, y = -log10(wi.eBH), color = sig)
  ) +
    geom_point(alpha = 0.8, size = 3) +
    theme_minimal(base_size = 16) +
    labs(
      title = title,
      x = expression(paste("Median ", log[2], " Difference")),
      y = "-log10(Median q value)",
      color = "Significance"
    ) +
    scale_color_manual(values = c("Significant" = "#d7191c",
                                  "Not Significant" = "#2b83ba")) +
    geom_hline(yintercept = -log10(sig_cutoff), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = subset(x.all_annot, wi.eBH < sig_cutoff),
      aes(label = TaxonName),
      size = 3,
      max.overlaps = 15,  # reduce clutter
      box.padding = 0.35,
      point.padding = 0.3
    ) +
    ylim(0, max_y) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )
  
  # 9. Print the volcano plot
  print(volcano_plot)
}