library(phyloseq)
library(qiime2R)
library(ggplot2)
library(vegan)
library(patchwork)
library(ALDEx2)
library(dplyr)
library(ggrepel)
library(tidytext)

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

# 1) Extract count matrix and ensure OTUs are rows
otu_tab <- as(otu_table(fungi_only), "matrix")
if (!taxa_are_rows(fungi_only)) {
  otu_tab <- t(otu_tab)
}
otu_df <- as.data.frame(otu_tab)

# 2) Extract taxonomy and collapse ranks
tax_df <- as.data.frame(tax_table(fungi_only), stringsAsFactors = FALSE)
tax_df$taxonomy <- apply(
  tax_df[, c("Kingdom","Phylum","Class","Order","Family","Genus","Species")],
  1, 
  function(x) paste(x[!is.na(x) & x != ""], collapse = ";")
)

# 3) Combine into one data.frame
funguild_in <- cbind(
  `#OTU ID` = rownames(otu_df),
  otu_df,
  taxonomy = tax_df[rownames(otu_df), "taxonomy"]
)

# 4) Write to disk
write.table(
  funguild_in,
  file = "funguild_input.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##### BASH COMMAND TO RUN FUNGUILD
nohup Guilds_v1.1.py \
-otu funguild_input.tsv \
-db fungi \
-u

### RELATIVE ABUNDANCE

# 1. Read in the guild‐counts (rows = OTUs, cols = samples + taxonomy)
guilds <- read.delim("/Users/michaelyoshida/Downloads/funguild_input.guilds.txt",
                     header=TRUE, row.names=1,
                     check.names=FALSE, sep="\t")

# 2. Identify which columns are your numeric sample‐counts
num_cols <- sapply(guilds, is.numeric)

# 3. Compute relative abundances (each column sums to 1)
rel_abund <- sweep(guilds[, num_cols],
                   2,
                   colSums(guilds[, num_cols]),
                   FUN = "/")

# 4. Put taxonomy back on the right
normalized <- cbind(rel_abund, guilds[, !num_cols])

# 5. Write out
write.table(normalized,
            file = "/Users/michaelyoshida/Downloads/funguild_input.guilds.rel.tsv",
            sep  = "\t",
            quote= FALSE,
            col.names = NA)

  # FUNGuild grouped categories boxplots (matches your FAPROTAX style)
  
  library(tidyverse)
  library(readr)
  library(stringr)
  library(forcats)
  library(stringi)
  
  # ===== 1) Read table (rows = guilds, cols = samples; already relative %) =====
  guild_tab <- read.delim(
    "funguild_input.guilds.rel.tsv",
    header      = TRUE,
    sep         = "\t",
    check.names = FALSE,
    row.names   = 1,
    quote       = "",
    comment.char= ""
  )
  
  # If export includes a "Guild" column, move it to rownames safely
  if ("Guild" %in% colnames(guild_tab)) {
    rn <- make.unique(as.character(guild_tab$Guild))
    guild_tab <- guild_tab %>% select(-Guild)
    rownames(guild_tab) <- rn
  }
  
  # ===== 2) Identify sample columns (only 2018/2022) & coerce JUST those to numeric =====
  sample_cols <- grep("^F(2018|2022)[._-]", colnames(guild_tab), value = TRUE)
  if (length(sample_cols) == 0) {
    # fallback: any F#### pattern
    sample_cols <- grep("^F\\d{4}[._-]", colnames(guild_tab), value = TRUE)
  }
  
  # Keep sample columns only (avoid taxonomy/notes), coerce to numeric, NA->0
  guild_num <- guild_tab[, sample_cols, drop = FALSE]
  guild_num[] <- lapply(guild_num, function(x) readr::parse_number(as.character(x)))
  guild_rel <- guild_num %>% mutate(across(everything(), ~ tidyr::replace_na(., 0)))
  
  # Drop guilds that are zero across all selected samples
  guild_rel <- guild_rel[rowSums(guild_rel, na.rm = TRUE) > 0, , drop = FALSE]
  
  # ===== 3) Long format + Season/Treatment parsing =====
  guild_meta <- guild_rel %>%
    rownames_to_column("Guild") %>%
    pivot_longer(all_of(sample_cols), names_to = "SampleID", values_to = "Abundance") %>%
    extract(SampleID, into = c("SeasonCode","RepCode"),
            regex = "(F\\d{4})[._-](.*)", remove = FALSE) %>%
    mutate(
      Season = recode(SeasonCode,
                      F2018 = "Fall 2018",
                      F2022 = "Fall 2022",
                      .default = NA_character_),
      Treatment = case_when(
        str_detect(RepCode, "C\\d+$") ~ "Control",
        str_detect(RepCode, "G\\d+$") ~ "Experimental",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Season), !is.na(Treatment))
  
  # ===== 4) Clean guild names (normalize; drop junk) =====
  normalize_guild <- function(x) {
    x <- as.character(x)
    x <- stri_trans_general(x, "NFKC")                          # Unicode normalize
    x <- str_replace_all(x, "[\\u2212\\u2010-\\u2015]", "-")    # fancy dashes -> "-"
    x <- str_replace_all(x, "[\\p{Cf}\\p{Cc}]", "")             # strip control/format chars
    x <- str_squish(x)
    x
  }
  
  guild_meta_clean <- guild_meta %>%
    mutate(Guild = normalize_guild(Guild)) %>%
    filter(
      !is.na(Guild),
      Guild != "",
      !Guild %in% c("Unassigned", "function not assigned", "NA"),
      str_detect(Guild, "[[:alpha:]]")                           # must contain a letter
    )
  
  # ===== 5) Assign broader groups by keyword =====
  # Order matters (first match wins). Add/tweak rules as needed.
  guild_meta_grouped <- guild_meta_clean %>%
    mutate(Group = case_when(
      str_detect(Guild, "Plant Pathogen")   ~ "Plant Pathogen",
      str_detect(Guild, "Animal Parasite")  ~ "Animal Pathogen",
      str_detect(Guild, "Lichen Parasite")  ~ "Lichen Parasite",
      str_detect(Guild, "Endophyte")        ~ "Endophyte",
      str_detect(Guild, "Saprotroph")       ~ "Saprotroph",
      TRUE                                  ~ "Other"
    ))
  
  # ===== 6) Aggregate abundance by group (per sample) =====
  plot_df <- guild_meta_grouped %>%
    group_by(Season, Treatment, SampleID, Group) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = fct_reorder(Group, Abundance, .fun = mean, .desc = TRUE))
  
  # ===== 7) Plot (FAPROTAX aesthetics; margins adjusted so labels don't clip) =====
  ggplot(plot_df, aes(x = Group, y = Abundance, fill = Treatment)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
    facet_wrap(~ Season, nrow = 2, scales = "fixed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw(base_size = 12) +
    theme(
      plot.margin     = margin(t = 5, r = 5, b = 5, l = 28),
      strip.background= element_rect(fill = "grey90", color = NA),
      strip.text      = element_text(face = "bold"),
      axis.text.x     = element_text(angle = 60, hjust = 1, size = 7),
      axis.title.x    = element_blank(),
      legend.position = "top"
    ) +
    labs(
      title    = "FUNGuild Functional Guild Abundance by Group",
      subtitle = "Fall 2018 vs Fall 2022 • Control vs Experimental",
      y        = "Relative abundance (%)",
      fill     = "Treatment"
    )

####################### ANOVA CALCULATIONS
  
  library(tidyverse)
  library(broom)
  library(stringr)
  
  in_path  <- "/Users/michaelyoshida/Downloads/funguild_input.guilds.rel.tsv"
  out_path <- "/Users/michaelyoshida/Downloads/two_way_anova_FUNGuild.tsv"
  
  # 1) Read normalized table with OTU/ASV IDs as rownames
  normalized <- read.delim(in_path, header = TRUE, check.names = FALSE, row.names = 1)
  
  # 2) Sample columns are numeric
  sample_cols <- normalized %>% select(where(is.numeric)) %>% names()
  
  # 3) Long format at ASV level; keep FUNGuild 'Guild'
  long_asv <- normalized %>%
    rownames_to_column("ASV") %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "RelAbundance")
  
  if (!"Guild" %in% names(long_asv)) stop("FUNGuild 'Guild' column not found in the TSV.")
  long_asv <- long_asv %>%
    mutate(Guild = if_else(is.na(Guild) | Guild == "", "Unassigned", Guild))
  
  # 4) Aggregate to Guild × Sample
  guild_sample <- long_asv %>%
    group_by(Guild, Sample) %>%
    summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
  
  # 5) Derive Season and treatment directly from Sample names
  guild_sample <- guild_sample %>%
    mutate(
      Season = case_when(
        str_detect(Sample, "^F2018") ~ "Fall 2018",
        str_detect(Sample, "^F2022") ~ "Fall 2022",
        TRUE ~ NA_character_
      ),
      treat_letter = str_match(Sample, "^F20(?:18|22)_[^A-Z]*([A-Z])")[,2],
      treatment = case_when(
        treat_letter == "C" ~ "Control",
        treat_letter == "G" ~ "Ginseng",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Season), !is.na(treatment)) %>%
    mutate(
      Season    = factor(Season, levels = c("Fall 2018", "Fall 2022")),
      treatment = factor(treatment, levels = c("Control", "Ginseng"))
    ) %>%
    select(-treat_letter)
  
  # 6) Two-way ANOVA per Guild (unbalanced allowed)
  anova_results <- guild_sample %>%
    group_by(Guild) %>%
    group_modify(~{
      df <- tidyr::drop_na(.x, Season, treatment, RelAbundance)
      if (nrow(df) < 4 || n_distinct(df$Season) < 2 || n_distinct(df$treatment) < 2 || var(df$RelAbundance) == 0) {
        return(tibble(term = NA_character_, df = NA_real_, sumsq = NA_real_,
                      meansq = NA_real_, statistic = NA_real_, p.value = NA_real_))
      }
      fit <- try(aov(RelAbundance ~ Season * treatment, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(tibble(term = NA_character_, df = NA_real_, sumsq = NA_real_,
                      meansq = NA_real_, statistic = NA_real_, p.value = NA_real_))
      }
      broom::tidy(fit) %>% filter(term != "Residuals")
    }) %>%
    ungroup() %>%
    arrange(Guild, term)
  
  # 7) BH/FDR correction PER TERM across guilds
  anova_fdr <- anova_results %>%
    group_by(term) %>%
    mutate(padj_BH = ifelse(!is.na(p.value), p.adjust(p.value, method = "BH"), NA_real_)) %>%
    ungroup()
  
  # 8) Export
  write.table(anova_fdr, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Also export a significant-only table (FDR < 0.05)
  out_sig <- sub("\\.tsv$", "_sig.tsv", out_path)
  anova_fdr %>%
    filter(!is.na(padj_BH), padj_BH < 0.05) %>%
    arrange(term, padj_BH, p.value) %>%
    write.table(file = out_sig, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Saved: ", out_path, " and ", out_sig)
  
#### ANOVA Visual
  # FUNGuild grouped boxplots — ONLY significant guilds (Season FDR<0.05),
  # with multi-guild labels split & abundance apportioned across splits
  
  library(tidyverse)
  library(readr)
  library(stringr)
  library(forcats)
  library(stringi)
  
  # ---- paths ----
  counts_path <- "/Users/michaelyoshida/Downloads/funguild_input.guilds.rel.tsv"
  anova_path  <- "/Users/michaelyoshida/Downloads/two_way_anova_FUNGuild.tsv"   # has padj_BH
  
  # ===== 0) Get significant guilds (Season FDR<0.05) =====
  anova_df <- read.delim(anova_path, header = TRUE, check.names = FALSE)
  if (!"padj_BH" %in% names(anova_df)) stop("Expected 'padj_BH' in ANOVA file. Use the BH-corrected TSV.")
  sig_guilds <- anova_df %>%
    filter(term == "Season", !is.na(padj_BH), padj_BH < 0.05) %>%
    pull(Guild) %>%
    unique()
  
  # ===== 1) Read counts (rows = FUNGuild labels, cols = samples; already relative) =====
  guild_tab <- read.delim(
    counts_path, header = TRUE, sep = "\t",
    check.names = FALSE, row.names = 1, quote = "", comment.char = ""
  )
  
  # If export includes a "Guild" column, move it to rownames safely
  if ("Guild" %in% colnames(guild_tab)) {
    rn <- make.unique(as.character(guild_tab$Guild))
    guild_tab <- guild_tab %>% select(-Guild)
    rownames(guild_tab) <- rn
  }
  
  # ===== 2) Identify 2018/2022 samples & coerce JUST those to numeric =====
  sample_cols <- grep("^F(2018|2022)[._-]", colnames(guild_tab), value = TRUE)
  if (length(sample_cols) == 0) {
    sample_cols <- grep("^F\\d{4}[._-]", colnames(guild_tab), value = TRUE)
  }
  
  guild_num <- guild_tab[, sample_cols, drop = FALSE]
  guild_num[] <- lapply(guild_num, function(x) readr::parse_number(as.character(x)))
  guild_rel <- guild_num %>% mutate(across(everything(), ~ tidyr::replace_na(., 0)))
  guild_rel <- guild_rel[rowSums(guild_rel, na.rm = TRUE) > 0, , drop = FALSE]
  
  # ===== 3) Long + Season/Treatment parsing =====
  guild_meta <- guild_rel %>%
    rownames_to_column("Guild") %>%
    pivot_longer(all_of(sample_cols), names_to = "SampleID", values_to = "Abundance") %>%
    extract(SampleID, into = c("SeasonCode","RepCode"),
            regex = "(F\\d{4})[._-](.*)", remove = FALSE) %>%
    mutate(
      Season = recode(SeasonCode, F2018 = "Fall 2018", F2022 = "Fall 2022", .default = NA_character_),
      Treatment = case_when(
        str_detect(RepCode, "C\\d+$") ~ "Control",
        str_detect(RepCode, "G\\d+$") ~ "Experimental",  # keep your label
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Season), !is.na(Treatment))
  
  # ===== 4) Clean guild names =====
  normalize_guild <- function(x) {
    x <- as.character(x)
    x <- stri_trans_general(x, "NFKC")
    x <- str_replace_all(x, "[\\u2212\\u2010-\\u2015]", "-")
    x <- str_replace_all(x, "[\\p{Cf}\\p{Cc}]", "")
    x <- str_squish(x)
    x
  }
  
  guild_meta_clean <- guild_meta %>%
    mutate(Guild = normalize_guild(Guild)) %>%
    filter(
      !is.na(Guild), Guild != "",
      !Guild %in% c("Unassigned", "function not assigned", "NA"),
      str_detect(Guild, "[[:alpha:]]")
    )
  
  # ===== 4b) Restrict to SIGNIFICANT guilds (match BEFORE splitting) =====
  guild_sig_only <- guild_meta_clean %>%
    filter(Guild %in% sig_guilds)
  
  # ===== 4c) Split multi-guild labels & apportion abundance across splits =====
  # Split on semicolons, commas, or hyphens used as separators, then divide Abundance equally.
  guild_split <- guild_sig_only %>%
    mutate(Guild_raw = Guild) %>%
    separate_rows(Guild, sep = "\\s*;\\s*|\\s*,\\s*|\\s*-\\s*") %>%
    mutate(Guild = str_squish(Guild)) %>%
    filter(Guild != "") %>%
    group_by(SampleID, Season, Treatment, Guild_raw) %>%
    mutate(n_parts = n()) %>%
    ungroup() %>%
    mutate(Abundance = Abundance / n_parts) %>%
    select(-n_parts, -Guild_raw)
  
  # ===== 5) Assign broader groups by keyword =====
  guild_grouped <- guild_split %>%
    mutate(Group = case_when(
      str_detect(Guild, "Plant Pathogen")   ~ "Plant Pathogen",
      str_detect(Guild, "Animal Parasite")  ~ "Animal Pathogen",
      str_detect(Guild, "Lichen Parasite")  ~ "Lichen Parasite",
      str_detect(Guild, "Endophyte")        ~ "Endophyte",
      str_detect(Guild, "Saprotroph")       ~ "Saprotroph",
      TRUE                                  ~ "Other"
    ))
  
  # ===== 6) Aggregate abundance by group (per sample) =====
  plot_df <- guild_grouped %>%
    group_by(Season, Treatment, SampleID, Group) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = fct_reorder(Group, Abundance, .fun = mean, .desc = TRUE))
  
  # ===== 7) Plot =====
  ggplot(plot_df, aes(x = Group, y = Abundance, fill = Treatment)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
    facet_wrap(~ Season, nrow = 2, scales = "fixed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw(base_size = 12) +
    theme(
      plot.margin     = margin(t = 5, r = 5, b = 5, l = 28),
      strip.background= element_rect(fill = "grey90", color = NA),
      strip.text      = element_text(face = "bold"),
      axis.text.x     = element_text(angle = 60, hjust = 1, size = 7),
      axis.title.x    = element_blank(),
      legend.position = "top"
    ) +
    labs(
      title    = "FUNGuild Groups (Significant Season Effects Only)",
      subtitle = "Fall 2018 vs Fall 2022 • Control vs Experimental",
      y        = "Relative abundance (%)",
      fill     = "Treatment"
    )
  
#######################
  
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(ggplot2)
  library(readr)
  
  target_seasons <- c("Fall 2018","Fall 2022")
  
  # --- paths ---
  counts_path <- "/Users/michaelyoshida/Downloads/funguild_input.guilds.rel.tsv"
  anova_path  <- "/Users/michaelyoshida/Downloads/two_way_anova_FUNGuild.tsv"  # must include padj_BH
  
  # --- 0) Which guilds are Season-significant? (FDR<0.05) + build symbol map ---
  anova_all <- read.delim(anova_path, header = TRUE, check.names = FALSE)
  
  sig_info <- anova_all %>%
    filter(term %in% c("Season","treatment","Season:treatment")) %>%
    mutate(sig = !is.na(padj_BH) & padj_BH < 0.05) %>%
    select(Guild, term, sig) %>%
    tidyr::pivot_wider(names_from = term, values_from = sig, values_fill = FALSE) %>%
    mutate(
      symbols = paste0(ifelse(Season, "*", ""),
                       ifelse(treatment, "†", ""),
                       ifelse(`Season:treatment`, "‡", "")),
      symbols = ifelse(symbols == "", "", paste0(" ", symbols))
    )
  
  season_sig <- sig_info %>% filter(Season) %>% pull(Guild) %>% unique()
  
  # --- 1) Read counts (rows = ASVs with a 'Guild' column; cols = samples) ---
  x <- read.delim(counts_path, header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1, quote = "", comment.char = "")
  
  sample_cols <- x %>% select(where(is.numeric)) %>% names()
  base <- x
  if (!"Guild" %in% names(base)) base <- base %>% tibble::rownames_to_column("Guild")
  base <- base %>% select(Guild, all_of(sample_cols))
  
  # --- 2) Long + derive Season; keep only Fall 2018/2022 ---
  long <- base %>%
    pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "RelAbund") %>%
    mutate(
      Season = case_when(
        str_detect(Sample, "^F2018") ~ "Fall 2018",
        str_detect(Sample, "^F2022") ~ "Fall 2022",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Season), Season %in% target_seasons)
  
  # --- 3) Keep ONLY Season-significant guilds (before splitting) ---
  long_sig <- long %>% filter(Guild %in% season_sig)
  
  # --- 4) Split multi-guild labels and APPORTION abundance across parts ---
  # Split on semicolons/commas/hyphens; divide by number of parts per Sample.
  long_split <- long_sig %>%                         # or `long` if not filtering by significance
    mutate(Guild_raw = Guild) %>%
    separate_rows(Guild, sep = "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*|\\s*-\\s*") %>%
    mutate(
      Guild = Guild %>%
        str_squish() %>%
        str_replace_all("[\\[\\]]", "") %>%         # remove [ and ]
        str_replace_all("^\\|+|\\|+$", "")          # trim leading/trailing pipes, just in case
    ) %>%
    filter(Guild != "") %>%
    group_by(Sample, Season, Guild_raw) %>%
    mutate(n_parts = n()) %>%
    ungroup() %>%
    mutate(RelAbund = RelAbund / n_parts)
  
  # Attach symbol string from the raw (pre-split) guild label
  long_split <- long_split %>%
    left_join(sig_info %>% select(Guild_raw = Guild, symbols), by = "Guild_raw") %>%
    select(-n_parts)
  
  # Build one symbol string per split Guild (collapse duplicates after cleaning)
  label_map <- long_split %>%
    distinct(Guild, symbols) %>%
    group_by(Guild) %>%
    summarise(
      symbols = paste0(
        if (any(str_detect(symbols %||% "", "\\*"))) "*" else "",
        if (any(str_detect(symbols %||% "", "†")))  "†" else "",
        if (any(str_detect(symbols %||% "", "‡")))  "‡" else ""
      ),
      .groups = "drop"
    ) %>%
    mutate(symbols = ifelse(symbols == "", "", paste0(" ", symbols)))
  
  # --- 5) Build wide safely (mean across samples per Season), compute change, back to long ---
  wide <- long_split %>%
    group_by(Guild, Season) %>%
    summarize(RelAbund = mean(RelAbund, na.rm = TRUE) * 100, .groups = "drop") %>%  # to %
    pivot_wider(names_from = Season, values_from = RelAbund)
  
  missing <- setdiff(target_seasons, names(wide))
  if (length(missing)) wide[missing] <- NA_real_
  
  summ <- wide %>%
    mutate(change_22_minus_18 = .data[["Fall 2022"]] - .data[["Fall 2018"]]) %>%
    pivot_longer(all_of(target_seasons), names_to = "Season", values_to = "RelAbund") %>%
    filter(!is.na(RelAbund)) %>%
    left_join(label_map, by = "Guild") %>%
    mutate(
      label  = paste0(Guild, coalesce(symbols, "")),
      Season = factor(Season, levels = target_seasons)
    )
  
  # --- 6) Plot (same style as FAPROTAX minus symbols) ---
  ggplot(summ, aes(x = RelAbund, y = fct_reorder(label, change_22_minus_18), fill = Season)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    labs(
      x = "Mean relative abundance (%)",
      y = "Guild",
      title = "FUNGuilds with Significant Season Effect (FDR < 0.05)",
      subtitle = "Labels: * Season, † Treatment, ‡ Season×Treatment (BH corrected)",
      caption  = "Multi-guild entries split; abundance apportioned across parts"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(hjust = 0), legend.position = "bottom")
  
########################### FINAL VISUALIZATION
  
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(ggplot2)
  library(readr)
  
  target_seasons <- c("Fall 2018","Fall 2022")
  
  # --- paths ---
  counts_path <- "/Users/michaelyoshida/Downloads/funguild_input.guilds.rel.tsv"
  anova_path  <- "/Users/michaelyoshida/Downloads/two_way_anova_FUNGuild.tsv"  # must include padj_BH
  
  # --- 0) Which guilds are Season-significant? (FDR<0.05) + build symbol map ---
  anova_all <- read.delim(anova_path, header = TRUE, check.names = FALSE)
  
  sig_info <- anova_all %>%
    filter(term %in% c("Season","treatment","Season:treatment")) %>%
    mutate(sig = !is.na(padj_BH) & padj_BH < 0.05) %>%
    select(Guild, term, sig) %>%
    pivot_wider(names_from = term, values_from = sig, values_fill = FALSE) %>%
    # Symbols for labels: * = Treatment, ‡ = Season×Treatment (omit Season symbol since all are Season-sig)
    mutate(
      symbols = paste0(
        ifelse(treatment, "*", ""),
        ifelse(`Season:treatment`, "‡", "")
      ),
      symbols = ifelse(symbols == "", "", paste0(" ", symbols))
    )
  
  season_sig <- sig_info %>% filter(Season) %>% pull(Guild) %>% unique()
  
  # --- 1) Read counts (rows = ASVs with a 'Guild' column; cols = samples) ---
  x <- read.delim(counts_path, header = TRUE, sep = "\t",
                  check.names = FALSE, row.names = 1, quote = "", comment.char = "")
  
  sample_cols <- x %>% select(where(is.numeric)) %>% names()
  base <- x
  if (!"Guild" %in% names(base)) base <- base %>% tibble::rownames_to_column("Guild")
  base <- base %>% select(Guild, all_of(sample_cols))
  
  # --- 2) Long + derive SeasonCode/RepCode → Season & Treatment; keep only Fall 2018/2022 ---
  long <- base %>%
    pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "RelAbund") %>%
    mutate(
      SeasonCode = str_extract(Sample, "^F(?:2018|2022)"),
      RepCode    = str_extract(Sample, "(?:C|G)\\d+$"),
      Season     = recode(SeasonCode, F2018 = "Fall 2018", F2022 = "Fall 2022", .default = NA_character_),
      Treatment  = case_when(
        str_detect(RepCode, "C\\d+$") ~ "Control",
        str_detect(RepCode, "G\\d+$") ~ "Experimental",
        TRUE                          ~ NA_character_
      )
    ) %>%
    filter(!is.na(Season), Season %in% target_seasons)
  
  # --- 3) Keep ONLY Season-significant guilds (before splitting) ---
  long_sig <- long %>% filter(Guild %in% season_sig)
  
  # --- 4) Split multi-guild labels and APPORTION abundance across parts ---
  long_split <- long_sig %>%
    mutate(Guild_raw = Guild) %>%
    separate_rows(Guild, sep = "\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*|\\s*-\\s*") %>%
    mutate(
      Guild = Guild %>%
        str_squish() %>%
        str_replace_all("[\\[\\]]", "") %>%
        str_replace_all("^\\|+|\\|+$", "")
    ) %>%
    filter(Guild != "") %>%
    group_by(Sample, Season, Treatment, Guild_raw) %>%
    mutate(n_parts = n()) %>%
    ungroup() %>%
    mutate(RelAbund = RelAbund / n_parts)
  
  # Attach symbol string from the raw (pre-split) guild label
  long_split <- long_split %>%
    left_join(sig_info %>% select(Guild_raw = Guild, symbols), by = "Guild_raw") %>%
    select(-n_parts)
  
  # Build one symbol string per cleaned Guild
  label_map <- long_split %>%
    distinct(Guild, symbols) %>%
    group_by(Guild) %>%
    summarize(
      symbols = paste0(
        if (any(str_detect(replace_na(symbols, ""), "\\*"))) "*" else "",
        if (any(str_detect(replace_na(symbols, ""), "‡")))   "‡" else ""
      ),
      .groups = "drop"
    ) %>%
    mutate(symbols = ifelse(symbols == "", "", paste0(" ", symbols)))
  
  # --- 5) Ordering key based on overall change (pooled across Treatment) ---
  order_map <- long_split %>%
    group_by(Guild, Season) %>%
    summarize(RelAbund = mean(RelAbund, na.rm = TRUE) * 100, .groups = "drop") %>%
    pivot_wider(names_from = Season, values_from = RelAbund) %>%
    mutate(change_22_minus_18 = .data[["Fall 2022"]] - .data[["Fall 2018"]]) %>%
    select(Guild, change_22_minus_18)
  
  # --- 6) Summaries per Season × Treatment and plot ---
  summ_treat <- long_split %>%
    group_by(Guild, Season, Treatment) %>%
    summarize(RelAbund = mean(RelAbund, na.rm = TRUE) * 100, .groups = "drop") %>%
    left_join(order_map, by = "Guild") %>%
    left_join(label_map, by = "Guild") %>%
    mutate(
      label     = paste0(Guild, coalesce(symbols, "")),
      Season    = factor(Season, levels = target_seasons),
      Treatment = fct_relevel(Treatment, c("Control","Experimental"))
    )
  
  ggplot(summ_treat, aes(x = RelAbund, y = fct_reorder(label, change_22_minus_18), fill = Treatment)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    facet_wrap(~ Season, nrow = 1, scales = "fixed") +
    labs(
      x = "Mean relative abundance (%)",
      y = "Guild",
      title = "FUNGuilds with Significant Season Effect (FDR < 0.05)",
      subtitle = "Labels next to names: * Treatment, ‡ Season×Treatment (BH corrected).",
      caption  = "Multi-guild entries split; abundance apportioned across parts"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.spacing.x = grid::unit(30, "pt"),  # add space between the two facets
      panel.border    = element_rect(color = "grey85", fill = NA)  # optional: clearer separation
    ) +
    theme(
      axis.text.y      = element_text(hjust = 0),
      legend.position  = "bottom",
      legend.title     = element_blank()
    )
  