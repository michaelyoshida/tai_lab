### BASE SCRIPT FOR BACTERIA WORK --> CREATES ALL PHYLOSEQ OBEJCTS FOR DOWNSTREAM ANALYSIS
### CREATED BY MICHAEL YOSHIDA FOR THE TAI LAB, 2025.

#clustered at 99%

library(phyloseq)
library(qiime2R)

setwd("/Users/michaelyoshida/Downloads")

ps <- qza_to_phyloseq(
  features="fwd-table-openref-99.qza",
  taxonomy="16S_fwd_openref99_taxonomy.qza", 
  metadata = "forward_read_metadata.tsv"
)

#taxonomy

ps_prune_r <- subset_taxa(ps, Kingdom!="Eukaryota", Kingdom!="NA")
ps_prune_f <- subset_taxa(ps_prune_r, Order!="Chloroplast" & Family!="Mitochondria")

ps_prune_fp <- prune_samples(sample_sums(ps_prune_f) >= 1000, ps_prune_f)
prune_reads_sample <- as.data.frame(sample_sums(ps_prune_fp))

psF <- prune_samples(
  ! sample_names(ps_prune_fp) %in% c("B1-Plus","B2-Plus"),
  ps_prune_fp)

orig_taxa <- taxa_names(psF)

# B) Extract the taxonomy table as a data.frame
tax.clean <- as.data.frame(
  tax_table(psF),
  stringsAsFactors = FALSE
)

# C) Define patterns to strip out
bad_patterns <- c(
  "(?i)\\buncultured.*",
  "(?i)\\bunidentified.*",
  "(?i)\\bmetagenome.*"
)

# D) Sweep across all taxonomic ranks and remove unwanted bits
tax.clean <- as.data.frame(
  lapply(tax.clean, function(col) {
    for (pat in bad_patterns) {
      col <- gsub(pat, "", col, perl = TRUE)
    }
    # Trim leftover underscores or spaces
    gsub("^[_ ]+|[_ ]+$", "", col)
  }),
  stringsAsFactors = FALSE
)

#    replace any NA with "" so that =="" tests work
tax.clean[is.na(tax.clean)] <- ""

# E) Fill down missing ranks
for (i in seq_len(nrow(tax.clean))) {
  if (tax.clean[i, "Phylum"] == "") {
    newname <- paste0("Unidentified_", tax.clean[i, "Kingdom"])
    tax.clean[i, c("Phylum","Class","Order","Family","Genus","Species")] <- newname
    
  } else if (tax.clean[i, "Class"] == "") {
    newname <- paste0("Unidentified_", tax.clean[i, "Phylum"])
    tax.clean[i, c("Class","Order","Family","Genus","Species")] <- newname
    
  } else if (tax.clean[i, "Order"] == "") {
    newname <- paste0("Unidentified_", tax.clean[i, "Class"])
    tax.clean[i, c("Order","Family","Genus","Species")] <- newname
    
  } else if (tax.clean[i, "Family"] == "") {
    newname <- paste0("Unidentified_", tax.clean[i, "Order"])
    tax.clean[i, c("Family","Genus","Species")] <- newname
    
  } else if (tax.clean[i, "Genus"] == "") {
    newname <- paste0("Unidentified_", tax.clean[i, "Family"])
    tax.clean[i, c("Genus","Species")] <- newname
    
  } else if (tax.clean[i, "Species"] == "") {
    tax.clean[i, "Species"] <- paste0("Unidentified_", tax.clean[i, "Genus"])
  }
}

# F) Restore original taxa names and ensure correct ordering
rownames(tax.clean) <- orig_taxa
tax.clean <- tax.clean[orig_taxa, ]

# G) Assign cleaned taxonomy back into the phyloseq object
tax_table(psF) <- as.matrix(tax.clean)


psF_names <- sample_names(psF)

# Create a logical mask for sample names that match the 2018 pattern
mask_2018 <- grepl("^F18[1-3][CG]\\d+$", psF_names)

# Replace names matching the pattern using a regex
psF_names_new <- psF_names  # start with the current names

# Remove a trailing "B" 
psF_names_new <- sub("B$", "", psF_names_new)
sample_names(psF) <- psF_names_new

psF_names_new[mask_2018] <- sub("^F18([1-3])([CG])(\\d+)$", "F2018-\\1\\2\\3", psF_names[mask_2018])

# Assign the updated names back to the phyloseq object
sample_names(psF) <- psF_names_new

sample_data(psF)$season <- sub("_.*", "", sample_data(psF)$season)


bacteria_samples <- as.data.frame(sample_data(psF))
bacteria_samples$seasonYear <- paste0(bacteria_samples$year, "-", bacteria_samples$season)
bacteria_samples$seasonYear <- factor(bacteria_samples$seasonYear, 
                                      levels = unique(bacteria_samples$seasonYear[order(bacteria_samples$year, 
                                                                                        factor(bacteria_samples$season, 
                                                                                               levels = c("spring", "summera", "summerb","fall")))]))

bacteria_samples$year <- as.character(bacteria_samples$year)

sample_data(psF) <- bacteria_samples

sample_data(psF)

bacteria_rare = rarefy_even_depth(psF)

site1_rare = subset_samples(bacteria_rare, site!="2" & site!="3")
site2_rare = subset_samples(bacteria_rare, site!="1" & site!="3")
site3_rare = subset_samples(bacteria_rare, site!="2" & site!="1")

site1_control_rare = subset_samples(bacteria_rare, site!="2" & site!="3" & treatment!="experimental")
site1_garden_rare = subset_samples(bacteria_rare, site!="2" & site!="3" & treatment!="control")

site2_control_rare = subset_samples(bacteria_rare, site!="1" & site!="3" & treatment!="experimental")
site2_garden_rare = subset_samples(bacteria_rare, site!="1" & site!="3" & treatment!="control")

site3_control_rare = subset_samples(bacteria_rare, site!="2" & site!="1" & treatment!="experimental")
site3_garden_rare = subset_samples(bacteria_rare, site!="2" & site!="1" & treatment!="control")

garden_rareF18F22 <- subset_samples(bacteria_rare, treatment!="control" & season == "fall" & year %in% c(2018,2022))
control_rareF18F22 <- subset_samples(bacteria_rare, treatment!="garden" & season == "fall" & year %in% c(2018,2022))

site1_experimental_rare_F18F22 <-subset_samples(site1_garden_rare, season == "fall" & year %in% c(2018,2022))
site2_experimental_rare_F18F22 <-subset_samples(site2_garden_rare, season == "fall" & year %in% c(2018,2022))
site3_experimental_rare_F18F22 <-subset_samples(site3_garden_rare, season == "fall" & year %in% c(2018,2022))

site1_control_rare_F18F22 <-subset_samples(site1_control_rare, season == "fall" & year %in% c(2018,2022))
site2_control_rare_F18F22 <-subset_samples(site2_control_rare, season == "fall" & year %in% c(2018,2022))
site3_control_rare_F18F22 <-subset_samples(site3_control_rare, season == "fall" & year %in% c(2018,2022))


experimental_data <- subset_samples(psF, treatment!="control")
control_data <- subset_samples(psF, treatment!="experimental")

experimental_F18F22 <- subset_samples(experimental_data, season == "fall" & year %in% c(2018,2022))
control_F18F22 <- subset_samples(control_data, season == "fall" & year %in% c(2018,2022))

site1_experimental <- subset_samples(psF, site!="2" & site!="3" & treatment!="control")
site2_experimental <- subset_samples(psF, site!="1" & site!="3" & treatment!="control")
site3_experimental <- subset_samples(psF, site!="1" & site!="2" & treatment!="control")

site1_control <- subset_samples(psF, site!="2" & site!="3" & treatment!="garden")
site2_control <- subset_samples(psF, site!="1" & site!="3" & treatment!="garden")
site3_control <- subset_samples(psF, site!="1" & site!="2" & treatment!="garden")

site1 <- subset_samples(psF, site!="2" & site!="3" )
site2 <- subset_samples(psF, site!="1" & site!="3" )
site3 <- subset_samples(psF, site!="1" & site!="2" )

site1_experimental_F18F22 <-subset_samples(site1_experimental, season == "fall" & year %in% c(2018,2022))
site2_experimental_F18F22 <-subset_samples(site2_experimental, season == "fall" & year %in% c(2018,2022))
site3_experimental_F18F22 <-subset_samples(site3_experimental, season == "fall" & year %in% c(2018,2022))

site1_control_F18F22 <-subset_samples(site1_control, season == "fall" & year %in% c(2018,2022))
site2_control_F18F22 <-subset_samples(site2_control, season == "fall" & year %in% c(2018,2022))
site3_control_F18F22 <-subset_samples(site3_control, season == "fall" & year %in% c(2018,2022))

plot_shannon_alphadiv <- function(phyloseq_object, column, plot_title, facet_by = NULL) {
  # 1) compute even-ness
  alpha_div <- estimate_richness(phyloseq_object, measures = "Jaccard")
  
  # 2) bring in x-axis grouping
  alpha_div[[column]] <- sample_data(phyloseq_object)[[column]]
  
  # 3) pull raw treatment values
  raw_treatment <- as.character(sample_data(phyloseq_object)$treatmentnum)
  raw_site <- as.character(sample_data(phyloseq_object)$site)
  
  # 4) if facetting by treatment, recode that column;
  #    otherwise just store a recoded 'treatment' for plotting
  if (!is.null(facet_by) && facet_by == "treatmentnum") {
    alpha_div[[facet_by]] <- factor(
      raw_treatment,
      levels = c("1", "2"),
      labels = c("Control", "Garden")
    )
  } 
  if (!is.null(facet_by) && facet_by == "site") {
    alpha_div[[facet_by]] <- factor(
      raw_site,
      levels = c("1", "2", "3"),
      labels = c("Site 1", "Site 2", "Site 3")
    )
  }
  
  # 6) build the plot
  p <- ggplot(alpha_div, aes_string(x = column, y = "Shannon")) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggtitle(plot_title) +
    theme_classic(base_size = 16) +
    theme(
      plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted")
    ) +
    labs(
      x     = column,
      y     = "Shannon Diversity Index",
      color = "Treatment"
    )
  
  # 7) add faceting if requested
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  return(p)
}



plot_jaccard <- function(physeq, column, plot_title, facet_by = NULL) {
  # 1) compute Jaccard distance + PCoA
  jacc_dist <- distance(physeq, method = "jaccard", binary = TRUE)
  ord       <- ordinate(physeq, method = "PCoA", distance = jacc_dist)
  pct1 <- round(ord$values$Relative_eig[1] * 100, 2)
  # 2) extract PC1 coordinates
  vecs      <- ord$vectors
  scores_df <- data.frame(PC1 = vecs[,1],
                          sample = rownames(vecs),
                          stringsAsFactors = FALSE)
  
  # 3) pull in your grouping column AS A VECTOR
  sd   <- sample_data(physeq)
  grp  <- sd[[column]]                # this is now a vector
  names(grp) <- rownames(sd)          # name it so we can subset by sample name
  scores_df[[column]] <- grp[scores_df$sample]
  
  # 4) optional facet recoding (I put this only to debug)
  if (!is.null(facet_by)) {
    raw <- sd[[facet_by]]
    names(raw) <- rownames(sd)
    raw <- raw[scores_df$sample]
    
    if (facet_by == "treatmentnum") {
      lvls <- c("1","2"); labs <- c("Control","Garden")
    } else if (facet_by == "site") {
      lvls <- c("1","2","3"); labs <- c("Site 1","Site 2","Site 3")
    } else {
      lvls <- sort(unique(raw)); labs <- lvls
    }
    scores_df[[facet_by]] <- factor(as.character(raw),
                                    levels = lvls,
                                    labels = labs)
  }
  
  # 5) plot
  p <- ggplot(scores_df, aes_string(x = column, y = "PC1", color = column)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggtitle(plot_title) +
    theme_classic(base_size = 16) +
    theme(
      plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted")
    ) +
    labs(
      x     = column,
      y     = paste0("PCoA1 (", pct1, "%)"),
      color = column
    )
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_by)))
  }
  
  return(p)
}

plot_evenness <- function(phyloseq_object, column, plot_title, facet_by = NULL) {
  # 1) compute Shannon H and observed richness
  alpha_div <- estimate_richness(phyloseq_object,
                                 measures = c("Shannon", "Observed"))
  
  # 2) calculate Pielou’s evenness
  alpha_div$Evenness <- alpha_div$Shannon / log(alpha_div$Observed)
  
  # 3) bring in x-axis grouping
  alpha_div[[column]] <- sample_data(phyloseq_object)[[column]]
  
  # 4) optional facet recoding
  raw_treatment <- as.character(sample_data(phyloseq_object)$treatmentnum)
  raw_site      <- as.character(sample_data(phyloseq_object)$site)
  
  if (!is.null(facet_by) && facet_by == "treatmentnum") {
    alpha_div[[facet_by]] <- factor(
      raw_treatment,
      levels = c("1", "2"),
      labels = c("Control", "Garden")
    )
  }
  
  if (!is.null(facet_by) && facet_by == "site") {
    alpha_div[[facet_by]] <- factor(
      raw_site,
      levels = c("1", "2", "3"),
      labels = c("Site 1", "Site 2", "Site 3")
    )
  }
  
  # 5) build the plot (y = Evenness)
  p <- ggplot(alpha_div, aes_string(x = column, y = "Evenness")) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggtitle(plot_title) +
    theme_classic(base_size = 16) +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted")
    ) +
    labs(
      x = column,
      y = "Pielou’s Evenness (J)",
      color = column
    )
  
  # 6) add facets if requested
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  return(p)
}


s1 <- plot_jaccard(site1_rare, facet_by = 'treatment', column = "seasonYear", plot_title = "Site 1 Rare Bacteria Jacard Test")
s1

s1 <- plot_evenness(site1_rare, facet_by = 'treatmentnum', column = "seasonYear", plot_title = "Site 1 Rare Bacteria Openref Even-ness Test")
s1

s1 <- plot_evenness(garden_rareF18F22, facet_by = 'site', column = "seasonYear", plot_title = "All Sites F18F22 Bacteria Openref Evenness Test, Garden")
s1

s1 <- plot_evenness(control_rareF18F22, facet_by = 'site', column = "seasonYear", plot_title = "All Sites F18F22 Bacteria Openref Evenness Test, Control")
s1

s1 <- plot_evenness(ginseng_rare_site1, facet_by = 'treatmentnum', column = "seasonYear", plot_title = "Site 1 Rare Fungi Even-ness Test")
s1

s1 <- plot_evenness(ginseng_rare_beginning_end_garden, facet_by = 'site', column = "seasonYear", plot_title = "All Sites F18F22 Fungi Evenness Test, Garden")
s1

s1 <- plot_evenness(ginseng_rare_beginning_end_control, facet_by = 'site', column = "seasonYear", plot_title = "All Sites F18F22 Fungi Evenness Test, Control")
s1

### ALPHA DIVERSITY
s1 <- plot_shannon_alphadiv(site1_rare, "seasonYear", "Shannon Alpha Diversity by Year, Open Ref Fwd Read Bacteria, Site 1", facet_by = "treatmentnum")
s1

f18f22garden <- plot_shannon_alphadiv(garden_rareF18F22, "seasonYear", "Shannon Alpha Diversity Fall 2018 vs Fall 2022, Open Ref Fwd Read Bacteria, Garden", facet_by = "site")
f18f22garden

f18f22control <-plot_shannon_alphadiv(control_rareF18F22, "seasonYear", "Shannon Alpha Diversity Fall 2018 vs Fall 2022, Open Ref Fwd Read Bacteria, Control", facet_by = "site")
f18f22control


s1 <- plot_shannon_alphadiv(site1_rare, "seasonYear", "Shannon Alpha Diversity by Year, Bacteria, Site 1, OpenRef", facet_by = "treatment")
s1

f18f22garden <- plot_shannon_alphadiv(garden_rareF18F22, "seasonYear", "Shannon Alpha Diversity Fall 2018 vs Fall 2022, Bacteria, OpenRef, Experimental", facet_by = "site")
f18f22garden

f18f22control <-plot_shannon_alphadiv(control_rareF18F22, "seasonYear", "Shannon Alpha Diversity Fall 2018 vs Fall 2022, Bacteria, OpenRef, Control", facet_by = "site")
f18f22control




beta_pcoa_plot <- function(physeq_obj, title, show_legend = TRUE, shape = "year") {
  # 1) Compute Bray–Curtis distance and PCoA
  betaDiversity <- distance(physeq_obj, method = "bray")
  ordination   <- ordinate(physeq_obj, method = "PCoA", distance = betaDiversity)
  
  # 2) % variance explained
  eigenvalues <- ordination$values$Relative_eig * 100
  x_variation <- round(eigenvalues[1], 2)
  y_variation <- round(eigenvalues[2], 2)
  
  shape_mapping <- c(
    "2018" = 21,  # solid circle
    "2019" = 15,  # solid square
    "2020" = 3,  # solid diamond CHANGE TO PLUS potentially
    "2021" = 17,  # solid triangle up
    "2022" = 25   # solid triangle down
  )
  
  # 4) Build the plot
  p <- plot_ordination(
    physeq_obj,
    ordination,
    shape = shape
  ) +
    geom_point(
      size   = 3,
      alpha  = 0.8,
      stroke = 0.8
    ) +
    scale_shape_manual(values = shape_mapping) +
    ggtitle(title) +
    theme_classic(base_size = 16) +
    theme(
      legend.position  = ifelse(show_legend, "bottom", "none"),
      legend.title     = element_text(size = 14, face = "bold"),
      legend.text      = element_text(size = 12),
      plot.title       = element_text(size = 18, face = "bold",
                                      hjust = 0.5, margin = margin(b = 10)),
      axis.text.x      = element_text(size = 12),
      axis.text.y      = element_text(size = 12)
    ) +
    labs(
      x     = paste0("PCoA Axis 1 (", x_variation, "%)"),
      y     = paste0("PCoA Axis 2 (", y_variation, "%)"),
      shape = "Year"
    )
  
  return(p)
}

# extra plot: all F1822, shape for time, colour for site
beta_pcoa_plot <- function(physeq_obj,
                           title,
                           show_legend = TRUE,
                           shape       = "year") {
  # 1) PCoA
  bc_dist <- distance(physeq_obj, method = "bray")
  ord     <- ordinate(physeq_obj, method = "PCoA", distance = bc_dist)
  eigs    <- ord$values$Relative_eig * 100
  x_var   <- round(eigs[1], 2)
  y_var   <- round(eigs[2], 2)
  
  # 2) Ensure grouping columns are factors
  sd <- sample_data(physeq_obj)
  # shape column
  if (!shape %in% colnames(sd)) {
    stop(paste("Column", shape, "not found in sample_data"))
  }
  sd[[shape]] <- factor(sd[[shape]])
  # site column
  if ("site" %in% colnames(sd)) {
    sd[["site"]] <- factor(sd[["site"]])
  }
  sample_data(physeq_obj) <- sd
  
  # 3) Define shape mapping
  shape_mapping <- c(
    "2018" = 21,  # filled circle
    "2019" = 15,  # solid square
    "2020" = 3,   # plus
    "2021" = 17,  # solid triangle up
    "2022" = 25   # filled triangle down
  )
  
  # 4) Build the plot
  p <- plot_ordination(physeq_obj, ord,
                       color = "site",    # now discrete
                       shape = shape) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = shape_mapping) +
    ggtitle(title) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title    = element_text(size = 14, face = "bold"),
      legend.text     = element_text(size = 12),
      plot.title      = element_text(size = 18, face = "bold",
                                     hjust = 0.5,
                                     margin = margin(b = 10)),
      axis.text       = element_text(size = 12)
    ) +
    labs(
      x     = paste0("PCoA Axis 1 (", x_var, "%)"),
      y     = paste0("PCoA Axis 2 (", y_var, "%)"),
      color = "Site",
      shape = shape
    )
  
  return(p)
}

# extra plot: control and experimental all in one plot s1
beta_pcoa_plot <- function(physeq_obj,
                           title,
                           show_legend = TRUE) {
  # 1) PCoA
  bc_dist <- distance(physeq_obj, method = "bray")
  ord     <- ordinate(physeq_obj, method = "PCoA", distance = bc_dist)
  eigs    <- ord$values$Relative_eig * 100
  x_var   <- round(eigs[1], 2)
  y_var   <- round(eigs[2], 2)
  
  # 2) Ensure treatment is a factor
  sd <- sample_data(physeq_obj)
  if (!"treatment" %in% colnames(sd)) {
    stop("Column 'treatment' not found in sample_data")
  }
  sd[["treatment"]] <- factor(sd[["treatment"]])
  sample_data(physeq_obj) <- sd
  
  # 3) Auto‐generate a shape mapping for each treatment level
  treat_lvls   <- levels(sd[["treatment"]])
  # use filled symbols 21–25 (will recycle if you have >5 levels)
  shape_values <- setNames(21:(20 + length(treat_lvls)), treat_lvls)
  
  # 4) Plot with only shape aesthetic
  p <- plot_ordination(physeq_obj, ord, shape = "treatment") +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = shape_values) +
    ggtitle(title) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title    = element_text(size = 14, face = "bold"),
      legend.text     = element_text(size = 12),
      plot.title      = element_text(size = 18, face = "bold",
                                     hjust = 0.5, margin = margin(b = 10)),
      axis.text       = element_text(size = 12)
    ) +
    labs(
      x     = paste0("PCoA Axis 1 (", x_var, "%)"),
      y     = paste0("PCoA Axis 2 (", y_var, "%)"),
      shape = "Treatment"
    )
  
  return(p)
}
#bacteria
beta_pcoa_plot <- function(physeq_obj,
                           title,
                           show_legend = TRUE) {
  # 1) PCoA
  bc_dist <- distance(physeq_obj, method = "bray")
  ord     <- ordinate(physeq_obj, method = "PCoA", distance = bc_dist)
  eigs    <- ord$values$Relative_eig * 100
  x_var   <- round(eigs[1], 2)
  y_var   <- round(eigs[2], 2)
  
  # 2) Ensure treatment is a factor
  sd <- sample_data(physeq_obj)
  if (!"treatment" %in% colnames(sd)) {
    stop("Column 'treatment' not found in sample_data")
  }
  sd[["treatment"]] <- factor(sd[["treatment"]])
  sample_data(physeq_obj) <- sd
  
  # 3) Auto‐generate a shape mapping for each treatment level
  st_lvls     <- levels(sd[["treatment"]])
  shape_values <- setNames(21:(20 + length(st_lvls)), st_lvls)
  
  # 4) Build the plot (no color, only shape)
  p <- plot_ordination(physeq_obj, ord, shape = "treatment") +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = shape_values) +
    ggtitle(title) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title    = element_text(size = 14, face = "bold"),
      legend.text     = element_text(size = 12),
      plot.title      = element_text(size = 18, face = "bold",
                                     hjust = 0.5, margin = margin(b = 10)),
      axis.text       = element_text(size = 12)
    ) +
    labs(
      x     = paste0("PCoA Axis 1 (", x_var, "%)"),
      y     = paste0("PCoA Axis 2 (", y_var, "%)"),
      shape = "Sample Type"
    )
  
  return(p)
}



#statistical analysis: repeated measures ANOVA for alpha diversity significance or not
# --> adonis for beta diversity 


### BETA DIVERSITY
s1control <- beta_pcoa_plot(site1_rare, "PCoA Beta Diversity by Year, Site 1 Bacteria, OpenRef")
s1control

s1experimental <- beta_pcoa_plot(garden_rareF18F22, "PCoA Beta Diversity F18 v F22, Bacteria, OpenRef, Garden")
s1experimental

s1controlbrief <- beta_pcoa_plot(control_rareF18F22, "PCoA Beta Diversity F18 vs F22, Bacteria, OpenRef, Control")
s1controlbrief

s1experimentalbrief <- beta_pcoa_plot(site1_experimental_rare_F18F22, "PCoA Beta Diversity Fall 2018 vs Fall 2022, Bacteria, Site 1, OpenRef, Experimental")
s1experimentalbrief

s2controlbrief <- beta_pcoa_plot(site2_control_rare_F18F22, "PCoA Beta Diversity Fall 2018 vs Fall 2022, Bacteria, Site 2, OpenRef, Control")
s2controlbrief

s2experimentalbrief <- beta_pcoa_plot(site2_experimental_rare_F18F22, "PCoA Beta Diversity Fall 2018 vs Fall 2022, Bacteria, Site 2, OpenRef, Experimental")
s2experimentalbrief

s3controlbrief <- beta_pcoa_plot(site3_control_rare_F18F22, "PCoA Beta Diversity Fall 2018 vs Fall 2022, Bacteria, Site 3, OpenRef, Control")
s3controlbrief

s3experimentalbrief <- beta_pcoa_plot(site3_experimental_rare_F18F22, "PCoA Beta Diversity Fall 2018 vs Fall 2022, Bacteria, Site 3, OpenRef, Experimental")
s3experimentalbrief




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
  
  # 8. Plot volcano
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
  
  # 9. Print 
  print(volcano_plot)
}

aldex_taxa_plot(site1_experimental_F18F22, title = "Site 1, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Experimental", sample_col = "year", denom = "all")
aldex_taxa_plot(site2_experimental_F18F22, title = "Site 2, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Experimental", sample_col = "year", denom = "all")
aldex_taxa_plot(site3_experimental_F18F22, title = "Site 3, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Experimental", sample_col = "year", denom = "all")

aldex_taxa_plot(site1_control_F18F22, title = "Site 1, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")
aldex_taxa_plot(site2_control_F18F22, title = "Site 2, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")
aldex_taxa_plot(site3_control_F18F22, title = "Site 3, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")

aldex_taxa_plot(ginseng_rare_be_garden_site1, title = "Site 1, Fungi, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Garden", denom = "all", kingdom = "Fungi")
aldex_taxa_plot(ginseng_rare_be_garden_site2, title = "Site 2, Fungi, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Garden", denom = "all", kingdom = "Fungi")
aldex_taxa_plot(ginseng_rare_be_garden_site3, title = "Site 3, Fungi, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Garden", denom = "all", kingdom = "Fungi")

aldex_taxa_plot(site1_control_F18F22, title = "Site 1, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")
aldex_taxa_plot(site2_control_F18F22, title = "Site 2, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")
aldex_taxa_plot(site3_control_F18F22, title = "Site 3, Bacteria, Taxa Plot of Differential Abundance by Genus, F18 vs F22, Control", sample_col = "year", denom = "all")


aldexsitebacteria2 <- aldex_volcano_plot(site2_experimental_F18F22, title = "Site 2, Bacteria, Volcano Plot of Differential Abundance, F18 vs F22, Experimental", sample_col = "year", denom = "all")
aldexsitebacteria3 <- aldex_volcano_plot(site3_experimental_F18F22, title = "Site 3, Bacteria, Volcano Plot of Differential Abundance, F18 vs F22, Experimental", sample_col = "year", denom = "all")

btaxa1 <- average_merge_phyloseq(site1_rare)
t7 <- taxa_bar_plot(btaxa1, fill='Family',title="Site 1 Open Ref Bacteria, Taxa Bar Plots,", show_legend = TRUE, setting = 'grouped', kingdom = "Bacteria")
t7

s1 <- average_merge_phyloseq(garden_rareF18F22)
t7 <- taxa_bar_plot(s1, fill = 'Family', kingdom = 'Bacteria', title="Open Ref Bacteria Garden Samples, Taxa Bar Plots", show_legend = TRUE, setting = 'grouped')
t7

s1 <- average_merge_phyloseq(control_rareF18F22)
t7 <- taxa_bar_plot(s1, fill = "Family", title="Open Ref Bacteria Control Samples, Taxa Bar Plots", show_legend = TRUE, setting = 'grouped', kingdom = "Bacteria")
t7

btaxa1 <- average_merge_phyloseq(site1_rare)
t7 <- taxa_bar_plot(btaxa1, fill='Family',title="Site 1 Fungi, Taxa Bar Plots,", show_legend = TRUE, setting = 'grouped', kingdom = "Bacteria")
t7

s1 <- average_merge_phyloseq(garden_rareF18F22)
t7 <- taxa_bar_plot(s1, fill = 'Family', kingdom = 'Fungi', title="Fungi Garden Samples, Taxa Bar Plots", show_legend = TRUE, setting = 'grouped')
t7

s1 <- average_merge_phyloseq(control_rareF18F22)
t7 <- taxa_bar_plot(s1, fill = "Family", title="Fungi Control Samples, Taxa Bar Plots", show_legend = TRUE, setting = 'grouped', kingdom = "Fungi")
t7

ginseng_rare_site1_garden = subset_samples(ginseng_rare, site == 1 & treatment == 'garden')
s1 <- average_merge_phyloseq(ginseng_rare_site1_garden)
t7 <- taxa_bar_plot(s1, title="Site 1 Garden Only, Taxa Bar Plots, F18 and F22 ", show_legend = TRUE, setting = 'grouped')
t7

s1control <- beta_pcoa_plot(site1_rare, "PCoA Beta Diversity by Year, Site 1 Bacteria, OpenRef")
s1control
