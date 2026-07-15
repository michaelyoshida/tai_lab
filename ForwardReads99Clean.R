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

# Clean taxonomy labels
#
# - Remove ambiguous annotations (e.g. "uncultured", "metagenome")
# - Replace missing taxonomic ranks with more informative placeholders
#   (e.g. Unidentified_Firmicutes)
# - Preserve original ASV identifiers

ps_prune_r <- subset_taxa(ps, Kingdom!="Eukaryota", Kingdom!="NA")
ps_prune_f <- subset_taxa(ps_prune_r, Order!="Chloroplast" & Family!="Mitochondria")

ps_prune_fp <- prune_samples(sample_sums(ps_prune_f) >= 1000, ps_prune_f)
prune_reads_sample <- as.data.frame(sample_sums(ps_prune_fp))

psF <- prune_samples(
  ! sample_names(ps_prune_fp) %in% c("B1-Plus","B2-Plus"),
  ps_prune_fp)

orig_taxa <- taxa_names(psF)

# Extract the taxonomy table as a data.frame
tax.clean <- as.data.frame(
  tax_table(psF),
  stringsAsFactors = FALSE
)

# Define patterns to strip out
bad_patterns <- c(
  "(?i)\\buncultured.*",
  "(?i)\\bunidentified.*",
  "(?i)\\bmetagenome.*"
)

# Sweep across all taxonomic ranks and remove unwanted bits
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

# Fill down missing ranks
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

# Restore original taxa names and ensure correct ordering
rownames(tax.clean) <- orig_taxa
tax.clean <- tax.clean[orig_taxa, ]

# Assign cleaned taxonomy back into the phyloseq object
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

#check
sample_data(psF)

#rarefy phyloseq object
bacteria_rare = rarefy_even_depth(psF)

# create subset phyloseq objects for downstream analyis
site1_rare = subset_samples(bacteria_rare, site!="2" & site!="3")
site2_rare = subset_samples(bacteria_rare, site!="1" & site!="3")
site3_rare = subset_samples(bacteria_rare, site!="2" & site!="1")

site1_rare_F18F21 <- subset_samples(site1_rare, season == "fall" & year %in% c(2018,2021))

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
control_data <- subset_samples(psF, treatment!="garden")

experimental_F18F22 <- subset_samples(experimental_data, season == "fall" & year %in% c(2018,2022))
control_F18F22 <- subset_samples(control_data, season == "fall" & year %in% c(2018,2022))

experimental_F18F21 <- subset_samples(experimental_data, season == "fall" & year %in% c(2018,2021))
control_F18F21 <- subset_samples(control_data, season == "fall" & year %in% c(2018,2021))

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