library(phyloseq)
library(qiime2R)

setwd("/Users/michaelyoshida/Downloads")

# 1) Import
ps <- qza_to_phyloseq(
  features  = "fwd-table-openref-99.qza",
  taxonomy  = "16S_fwd_openref99_taxonomy.qza",
  metadata  = "forward_read_metadata.tsv"
)

# 2) Prune unwanted taxa & low-depth samples
ps_prune_r  <- subset_taxa(ps, Kingdom != "Eukaryota" & Kingdom != "NA")
ps_prune_f  <- subset_taxa(ps_prune_r, Order != "Chloroplast" & Family != "Mitochondria")
ps_prune_fp <- prune_samples(sample_sums(ps_prune_f) >= 1000, ps_prune_f)

psF <- prune_samples(!sample_names(ps_prune_fp) %in% c("B1-Plus","B2-Plus"), ps_prune_fp)

# 3) Extract OTU & TAX tables (aligned)
otu <- as.data.frame(otu_table(psF))
if (!taxa_are_rows(psF)) otu <- t(otu)

tax <- as.data.frame(tax_table(psF), stringsAsFactors = FALSE)

# 4) Full taxonomy string (original; keep if you still want it)
tax$taxonomy <- apply(tax, 1, function(x) paste(x, collapse = "; "))

# 5) Build a clean FAPROTAX key ------------------------------------
strip_rank <- function(x){
  x <- gsub("^[A-Za-z]__|D_[0-9]+__", "", x)                        # rank prefixes
  x <- gsub("^(Unidentified_|uncultured_|Candidatus_).*", "", x, ignore.case = TRUE)
  x <- gsub("_sensu_stricto.*", "", x)
  x <- gsub("_+$", "", x)
  trimws(x)
}

best_for_faprotax <- function(vec){
  vec <- strip_rank(vec)
  vec <- rev(vec)                          
  nm  <- vec[which(nzchar(vec))[1]]
  if (length(nm) == 0) NA_character_ else nm
}

tax$FAPRO_name <- apply(tax, 1, best_for_faprotax)

tax$FAPRO_name <- gsub("(?i)^(uncultured|unidentified).*", "", tax$FAPRO_name, perl = TRUE)
tax$FAPRO_name[nchar(tax$FAPRO_name) == 0] <- NA

# keep counts only
counts_only <- otu   # copy
# add the clean metadata column
counts_only$FAPRO_name <- tax[rownames(counts_only), "FAPRO_name"]

# order columns: #OTU ID, samples..., FAPRO_name
otu_out <- cbind(`#OTU ID` = rownames(counts_only),
                 counts_only[, c(sample_names(psF), "FAPRO_name")])

write.table(otu_out, "otu_fapro.tsv", sep="\t", quote=FALSE, row.names=FALSE)

####### BASH COMMANDS FOR GENERATING FAPROTAX DATA

source activate qiime2-2022.2

biom convert \
-i otu_fapro.tsv \
-o otu_fapro.h5.biom \
--table-type "OTU table" \
--to-hdf5 \
--process-obs-metadata naive

collapse_table.py \
-i otu_fapro.h5.biom \
-g /usr/local/bin/FAPROTAX.txt \
-n columns_before_collapsing \
--group_leftovers_as "function not assigned" \
--collapse_by_metadata FAPRO_name \
-o faprotax_out.txt \
--force -v

collapse_table.py -i otu_fapro.h5.biom -g /usr/local/bin/FAPROTAX.txt \
--collapse_by_metadata FAPRO_name -o faprotax_out_binary.txt \
-b --force -v

### VISUALIZATIONS

# read in FAPROTAX table (rows = functions, cols = samples)
func <- read.table("/Users/michaelyoshida/Downloads/faprotax_out.txt", header=TRUE, row.names=1, sep="\t")
# convert to % per sample
func_rel <- sweep(func, 2, colSums(func), FUN="/") * 100

library(reshape2); library(ggplot2)

mat_rel <- as.matrix(func_rel)

func_melt <- melt(mat_rel,
                  varnames   = c("Function","Sample"),
                  value.name = "Abundance")


############### FINAL VISUALIZATION (NO ANOVA CALCULATIONS)

library(tidyverse)    # loads dplyr, tidyr, ggplot2, stringr, etc.

# 2. Read in your collapsed FAPROTAX table
func_counts <- read.table(
  "faprotax_out.txt",
  header      = TRUE,
  sep         = "\t",
  row.names   = 1,
  check.names = FALSE
)

# 3. Convert counts → relative abundance (%) per sample
func_rel <- sweep(func_counts, 2, colSums(func_counts), FUN = "/") * 100

# 3a. Drop any function that is zero in *all* samples
func_rel <- func_rel[rowSums(func_rel) > 0, ]

# 4. Reshape to long form and extract metadata
func_meta <- func_rel %>%
  as.data.frame() %>%
  rownames_to_column("Function") %>%
  pivot_longer(
    cols      = -Function,
    names_to  = "SampleID",
    values_to = "Abundance"
  ) %>%
  # split SampleID like "F2018_1C1" → SeasonCode = "F2018", RepCode = "1C1"
  extract(
    col   = SampleID,
    into  = c("SeasonCode", "RepCode"),
    regex = "(F\\d{4})_(.*)"
  ) %>%
  mutate(
    # only recode F2018 and F2022, everything else becomes NA
    Season = recode(
      SeasonCode,
      F2018    = "Fall 2018",
      F2022    = "Fall 2022",
      .default = NA_character_
    ),
    # Control if RepCode ends in C<number>, Experimental if G<number>
    Treatment = case_when(
      str_detect(RepCode, "C\\d+$") ~ "Control",
      str_detect(RepCode, "G\\d+$") ~ "Garden",
      TRUE                          ~ NA_character_
    )
  ) %>%
  # keep only the seasons & treatments we care about
  filter(!is.na(Season), !is.na(Treatment))

# 5. Plot boxplots: one panel per Season, boxes for Control vs Experimental under each Function
ggplot(func_meta,
       aes(x    = Function,
           y    = Abundance,
           fill = Treatment)) +
  geom_boxplot(
    outlier.size = 0.5,
    position     = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ Season,
             nrow   = 2,
             scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 60, hjust = 1, size = 7),
    axis.title.x     = element_blank(),
    legend.position  = "top"
  ) +
  labs(
    title    = "FAPROTAX Functional Group Abundance",
    subtitle = "Fall 2018 vs Fall 2022 • Control vs Experimental",
    y        = "Relative abundance (%)",
    fill     = "Treatment"
  )


############################ ANOVA CALCULATION

library(dplyr)
library(rstatix)


anova_results <- func_meta %>%
  group_by(Function) %>%
  do(tidy(aov(Abundance ~ Season * Treatment, data = .))) %>%
  ungroup()

# 2. Adjust p‐values (BH) *within* each term (Season, Treatment, Season:Treatment)
anova_adj <- anova_results %>%
  group_by(term) %>%
  mutate(
    p.adj      = p.adjust(p.value, method = "BH"),
    p.signif   = symnum(p.adj, 
                        cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                        symbols   = c("***","**","*","ns"))
  ) %>%
  ungroup()

# 3. Pivot to wide so each term’s stats are columns
anova_summary <- anova_adj %>%
  select(Function, term, statistic, p.value, p.adj, p.signif) %>%
  pivot_wider(
    names_from  = term,
    values_from = c(statistic, p.value, p.adj, p.signif),
    names_glue  = "{term}_{.value}"
  )

# 4. Compute per‐Season×Treatment medians (just like before)
medians <- func_meta %>%
  group_by(Season, Function, Treatment) %>%
  summarize(median_abund = median(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from   = Treatment,
    values_from  = median_abund,
    names_prefix = "median_"
  )

# 5. Join medians to ANOVA table
report_final <- anova_summary %>%
  left_join(medians, by = "Function") %>%
  arrange(Season_p.adj)   # or pick whichever term you want to sort by

# 6. Write out
write.table(
  report_final,
  file      = "faprotax_statistics.tsv",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

########################################################################## OUTPUT ANOVA RESULTS W CORRECTION

# 1. Fit two‐way ANOVA for each Function
anova_results <- func_meta %>%
  group_by(Function) %>%
  do(tidy(aov(Abundance ~ Season * Treatment, data = .))) %>%
  ungroup()

# 2. Adjust p‐values (BH) *within* each term (Season, Treatment, Season:Treatment)
anova_adj <- anova_results %>%
  group_by(term) %>%
  mutate(
    p.adj      = p.adjust(p.value, method = "BH"),
    p.signif   = symnum(p.adj, 
                        cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                        symbols   = c("***","**","*","ns"))
  ) %>%
  ungroup()

# 3. Pivot to wide so each term’s stats are columns
anova_summary <- anova_adj %>%
  select(Function, term, statistic, p.value, p.adj, p.signif) %>%
  pivot_wider(
    names_from  = term,
    values_from = c(statistic, p.value, p.adj, p.signif),
    names_glue  = "{term}_{.value}"
  )

# 4. Compute per‐Season×Treatment medians (just like before)
medians <- func_meta %>%
  group_by(Season, Function, Treatment) %>%
  summarize(median_abund = median(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from   = Treatment,
    values_from  = median_abund,
    names_prefix = "median_"
  )

# 5. Join medians to ANOVA table
report_final <- anova_summary %>%
  left_join(medians, by = "Function") %>%
  arrange(Season_p.adj)   # or pick whichever term you want to sort by

# 6. Write out
write.table(
  report_final,
  file      = "faprotax_statistics.tsv",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

sig_time <- report_final %>%
  filter(Season_p.adj < 0.05)

# 2. Among those, flag functions also significant by Treatment,
#    but with no significant Season:Treatment interaction
starred <- sig_time %>%
  filter(Treatment_p.adj < 0.05,
         `Season:Treatment_p.adj` >= 0.05) %>%
  pull(Function)

# 3. Prepare data for plotting: –log10 of Season p-adj
plot_df <- sig_time %>%
  mutate(logp = -log10(Season_p.adj),
         star = if_else(Function %in% starred, "*", ""))

plot_df2 <- plot_df %>%
  # create a new label that appends "*" if starred
  mutate(label = if_else(star == "*",
                         paste0(Function, " *"),
                         Function))

ggplot(plot_df2, aes(x = logp, y = fct_reorder(label, logp))) +
  geom_col() +
  labs(
    x = expression(-log[10]*"(Season p.adj)"),
    y = "Function",
    title = "Functions with Significant Effect Between Fall 2018 and Fall 2022",
    subtitle = "Stars indicate also significant for Treatment (no interaction)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(hjust = 0)  # left-justify names
  )

# 5. Map ASVs to each starred function
asv_map <- tax %>%
  rownames_to_column("ASV") %>%
  select(ASV, FAPRO_name) %>%
  filter(FAPRO_name %in% starred) %>%
  arrange(FAPRO_name)

# 6. Print out the ASV–function mapping
asv_map

######################################## ANOVA VISUALIZATION (NO TREATMENT INCORPORATION)

library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(ggplot2)

target_seasons <- c("Fall 2018","Fall 2022")

# 1) Normalize whitespace in Season (kills non-breaking spaces, etc.)
func_meta <- func_meta %>%
  mutate(Season = str_replace_all(Season, "\\u00A0", " ")) %>%  # NBSP -> space
  mutate(Season = str_squish(Season))

# 2) Build wide safely
wide <- func_meta %>%
  semi_join(sig_time, by = "Function") %>%
  filter(Season %in% target_seasons, !is.na(Abundance)) %>%
  group_by(Function, Season) %>%
  summarize(RelAbund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Season, values_from = RelAbund)

# 3) If a season column is missing, create it as NA
missing <- setdiff(target_seasons, names(wide))
if (length(missing)) wide[missing] <- NA_real_

# 4) Compute change and go back to long
summ <- wide %>%
  mutate(change_22_minus_18 = .data[["Fall 2022"]] - .data[["Fall 2018"]]) %>%
  pivot_longer(all_of(target_seasons), names_to = "Season", values_to = "RelAbund") %>%
  filter(!is.na(RelAbund)) %>%
  mutate(
    label  = if_else(Function %in% starred, paste0(Function, " *"), Function),
    Season = factor(Season, levels = target_seasons)
  )

# 5) Plot
ggplot(summ, aes(x = RelAbund, y = fct_reorder(label, change_22_minus_18), fill = Season)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  labs(
    x = "Mean relative abundance (%)",
    y = "Function",
    title = "Functions with Significant Season Effect (Fall 2018 vs Fall 2022)",
    subtitle = "Stars indicate also significant for Treatment (no Season×Treatment interaction)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(hjust = 0), legend.position = "bottom")

############################# FINAL VERSION, COLOUR BY TREATMENT, FACET BY SEASON

library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)

target_seasons <- c("Fall 2018","Fall 2022")

# 1) Clean Season
func_meta <- func_meta %>%
  mutate(
    Season = str_replace_all(Season, "\u00A0", " "),
    Season = str_squish(Season)
  )

# 2) Treatment metadata
func_meta <- func_meta %>%
  mutate(
    # normalize any pre-existing values
    Treatment = recode(Treatment, "Ginseng" = "Garden"),
    Treatment = na_if(str_squish(Treatment), ""),
    
    # try to infer from Sample (e.g., F2018_C1 / F2022_G3)
    letter_sample = str_match(Sample, "^F20(?:18|22)_[^A-Za-z]*([A-Za-z])")[,2],
    
    # fallback from RepCode like "...C1" or "...G2" if available
    letter_rep = if ("RepCode" %in% names(.)) str_match(RepCode, "([CG])[0-9]+$")[,2] else NA_character_,
    
    letter = toupper(coalesce(letter_sample, letter_rep)),
    
    Treatment_inferred = case_when(
      letter == "C" ~ "Control",
      letter == "G" ~ "Garden",
      TRUE ~ NA_character_
    ),
    
    # fill only where Treatment is missing
    Treatment = coalesce(Treatment, Treatment_inferred)
  ) %>%
  select(-letter_sample, -letter_rep, -letter, -Treatment_inferred)


# 3) Build ordering by global change (avg over treatments)
wide_order <- func_meta %>%
  semi_join(sig_time, by = "Function") %>%
  filter(Season %in% target_seasons, !is.na(Abundance)) %>%
  group_by(Function, Season) %>%
  summarize(RelAbund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Season, values_from = RelAbund)

missing <- setdiff(target_seasons, names(wide_order))
if (length(missing)) wide_order[missing] <- NA_real_

wide_order <- wide_order %>%
  mutate(change_22_minus_18 = .data[["Fall 2022"]] - .data[["Fall 2018"]]) %>%
  select(Function, change_22_minus_18)

# 4) Summaries per Function × Season × Treatment
summ_treat <- func_meta %>%
  semi_join(sig_time, by = "Function") %>%
  filter(Season %in% target_seasons, !is.na(Abundance), !is.na(Treatment)) %>%
  group_by(Function, Season, Treatment) %>%
  summarize(RelAbund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(wide_order, by = "Function") %>%
  mutate(
    label     = if_else(Function %in% starred, paste0(Function, " *"), Function),
    Season    = factor(Season, levels = target_seasons),
    Treatment = factor(Treatment, levels = c("Control","Garden"))
  )

# 5) Plot: facet by Season, color/fill by Treatment
ggplot(summ_treat, aes(x = RelAbund, y = fct_reorder(label, change_22_minus_18))) +
  geom_col(aes(fill = Treatment), position = position_dodge(width = 0.7), width = 0.65) +
  facet_wrap(~ Season, nrow = 1, scales = "fixed") +
  labs(
    x = "Mean relative abundance (%)",
    y = "Function",
    title = "Functions with Significant Season Effect (Fall 2018 vs Fall 2022)",
    subtitle = "Faceted by Season; color = Treatment. * also significant for Treatment (no Season×Treatment interaction)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(hjust = 0), legend.position = "bottom")

