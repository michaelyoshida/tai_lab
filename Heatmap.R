### SCRIPT TO CREATE HEATMAP VISUALIZATIONS
### TWO FUNCTIONS: ONE SUBSETTING BY ALDEX2 RESULTS, ONE SHOWING RELATIVE ABUNDANCES OF TAXONOMY BASED ON PREDETERMINED LIST
### CREATED BY MICHAEL YOSHIDA FOR THE TAI LAB, 2025.

library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(dplyr)

# rank/filter top_n by aldex effect size
bacteria_aldex_genus_heatmap_CH <- function(
    physeq,
    sample_col     = "seasonYear",
    compare_times,                    # e.g., c("2018-fall","2021-fall")
    group_var      = "seasonYear",    # used for column labels
    facet_var      = "treatment",     # facet (column split)
    tax_level      = "Genus",
    sig_cutoff     = 0.05,            # ALDEx2 significance threshold (wi.eBH)
    mc.samples     = 128,
    denom          = "all",
    rand.seed      = 123,
    title          = "Genus-Level Heatmap (ALDEx2-filtered, faceted)",
    collapse_using = c("clean","original"),
    percent        = FALSE,           # show 0-1 (FALSE) or 0-100 (TRUE)
    top_n          = NULL,            # optionally show top N genera by mean rel. abundance
    facet_order    = NULL,            # character vector to order facets
    column_rotate  = 45,
    row_fontsize   = 9,
    col_fontsize   = 9,
    sort_by        = c("mean","abs_effect","alpha")  # NEW
){
  collapse_using <- match.arg(collapse_using)
  sort_by        <- match.arg(sort_by)
  
  # checks to make sure required arguments are provided (past debugging)
  stopifnot(inherits(physeq, "phyloseq"))
  if (is.null(compare_times) || length(compare_times) != 2) {
    stop("Please supply compare_times = c(time1, time2)")
  }
  
  sd_all <- as.data.frame(sample_data(physeq), stringsAsFactors = FALSE)
  if (!all(c(sample_col, group_var, facet_var) %in% colnames(sd_all))) {
    stop("sample_col, group_var, or facet_var not found in sample_data(physeq).")
  }
  if (!tax_level %in% colnames(as.data.frame(tax_table(physeq)))) {
    stop("`tax_level` not found in tax_table(physeq).")
  }
  
  # create a working a copy of phyloseq object for function to manipulate
  phy <- physeq
  
  # subset to the two time points for ALDEx2 
  sel     <- sample_data(phy)[[sample_col]] %in% compare_times
  phy_cmp <- prune_samples(sel, phy)
  
  otu_mat <- as(otu_table(phy_cmp), "matrix")
  if (!taxa_are_rows(phy_cmp)) otu_mat <- t(otu_mat)
  otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]
  groups  <- as.character(sample_data(phy_cmp)[[sample_col]])
  
  # ALDEx2 filtering at OTU level, run CLR transformation, statistical test, analyze effect sizes
  set.seed(rand.seed)
  clr  <- ALDEx2::aldex.clr(otu_mat, conds = groups, mc.samples = mc.samples, denom = denom)
  test <- ALDEx2::aldex.ttest(clr)
  eff  <- ALDEx2::aldex.effect(clr)
  
  # create res (result) object that has combined statistical-test and effect-size outputs
  res  <- cbind(as.data.frame(test), as.data.frame(eff))
  res$taxa <- rownames(res)
  
  # retain significant (by p-value threshold) OTUs
  sig_res  <- subset(res, wi.eBH <= sig_cutoff)
  sig_otus <- sig_res$taxa
  if (length(sig_otus) == 0) stop("No significant OTUs found between those two times.")
  
  # use full dataset so significant taxa can be heatmap visualized through entire timeline
  sd0 <- as.data.frame(sample_data(phy), stringsAsFactors = FALSE)
  
  # create unique labels for subsequent replicate agglomeration
  comb_lab <- paste(sd0[[group_var]], sd0[[facet_var]], sep = " | ")
  comb_df <- data.frame(
    .comb = comb_lab,
    Group = sd0[[group_var]],
    Facet = sd0[[facet_var]],
    row.names = rownames(sd0),
    check.names = FALSE
  )
  
  # calculate mean abundance within each group-by-facet combination (pooling the five replicates)
  phy_agg <- merge_samples(phy, group = comb_df$.comb, fun = mean)
  
  # reconstruct metadata table for pooled replicates 
  merged_keys <- unique(comb_df[, c(".comb","Group","Facet")])
  rownames(merged_keys) <- merged_keys$.comb; merged_keys$.comb <- NULL
  merged_keys[] <- lapply(merged_keys, as.character)
  sample_data(phy_agg) <- phyloseq::sample_data(merged_keys)
  
  # convert aggregated replicate counts into relative abundance
  phy_agg_ra <- transform_sample_counts(phy_agg, function(x) x / sum(x))
  
  # extract relative-abundance matrix in taxa (rows) by sample (column) orientation
  otu_full_mat <- as(otu_table(phy_agg_ra), "matrix")
  if (!taxa_are_rows(phy_agg_ra)) otu_full_mat <- t(otu_full_mat)
  
  # retain only significant OTUs that remain in aggregated replicates dataset
  sig_full_ra  <- otu_full_mat[intersect(sig_otus, rownames(otu_full_mat)), , drop = FALSE]
  if (!nrow(sig_full_ra)) stop("Significant OTUs not present after aggregation (check names).")
  tax_for_labels <- if (collapse_using == "clean") tax.clean else orig_tax_df
  tax_for_labels <- tax_for_labels[rownames(sig_full_ra), , drop = FALSE]
  
  # extract taxnomic labels used to group OTUs
  lvl_vec <- tax_for_labels[[tax_level]]
  
  #debug (generate labels for taxa missing various ranks)
  is_blank <- is.na(lvl_vec) | lvl_vec == ""
  if (any(is_blank)) {
    fam <- tax_for_labels[["Family"]]
    ord <- tax_for_labels[["Order"]]
    lvl_vec[is_blank] <- ifelse(
      !is.na(fam) & nzchar(fam),
      paste0("Unidentified_", fam),
      ifelse(
        !is.na(ord) & nzchar(ord),
        paste0("Unidentified_", ord),
        "Unidentified"
      )
    )
  }
  
  # sum the relative abundances of significant OTUs sharing the same taxonomic label (here it's genus)
  genus_rel_all <- rowsum(sig_full_ra, group = lvl_vec)
  genus_rel_all <- genus_rel_all[rowSums(genus_rel_all, na.rm = TRUE) > 0, , drop = FALSE]
  if (!nrow(genus_rel_all)) stop("All ALDEx2-significant genera have zero abundance after collapse.")
  
  # Convert to % if requested
  mat <- if (percent) genus_rel_all * 100 else genus_rel_all
  mat <- as.matrix(mat); storage.mode(mat) <- "numeric"
  
  # preparing taxa level effect sizes. Retain only significant OTUs to be represented in heatmap
  sig_res_present <- subset(sig_res, taxa %in% rownames(sig_full_ra))
  otu_to_genus <- data.frame(
    taxa  = rownames(sig_full_ra),
    genus = lvl_vec,
    stringsAsFactors = FALSE
  )
  sig_with_genus <- merge(sig_res_present, otu_to_genus, by = "taxa", all.x = TRUE)
  
  #summarize multiple OTU-level effects within same taxon using mean absolute ALDEx2 effect size
  gen_abs_effect <- tapply(
    abs(sig_with_genus$effect),
    sig_with_genus$genus,
    function(v) mean(v, na.rm = TRUE)
  )
  if (is.null(gen_abs_effect)) gen_abs_effect <- setNames(numeric(0), character(0))
  
  # determine top taxas (+ how many to show based on top_n) and row order 
  if (!is.null(top_n)) {
    
    ge_all <- gen_abs_effect[rownames(mat)]
    ge_all[is.na(ge_all)] <- -Inf
    
    keep_idx <- order(-ge_all, rownames(mat))
    keep_idx <- keep_idx[seq_len(min(top_n, length(keep_idx)))]
    mat <- mat[keep_idx, , drop = FALSE]
  }
  
  #calculate mean relative abundance in case of sorting by abundance
  row_means <- rowMeans(mat, na.rm = TRUE)
  
  #sort by alphabetical name
  if (sort_by == "alpha") {
    ord_rows <- order(rownames(mat))
    
    #sort by absolute ALDEx2 effect sizes
  } else if (sort_by == "abs_effect") {
    ge <- gen_abs_effect[rownames(mat)]
    ge[is.na(ge)] <- -Inf
    ord_rows <- order(-ge, -row_means, rownames(mat))
  
    # sort by mean relative abundance 
    } else { # "mean"
    ord_rows <- order(-row_means, rownames(mat))
  }
  mat <- mat[ord_rows, , drop = FALSE]
  
  # retrieve and align metadata with the heatmap matrix columns
  sd2 <- as.data.frame(sample_data(phy_agg_ra), stringsAsFactors = FALSE)
  sd2 <- sd2[colnames(mat), , drop = FALSE]
  
  # format facet names to title case
  cap <- function(x) tools::toTitleCase(tolower(x))
  sd2$Facet <- cap(sd2$Facet)
  
  # create facet order
  if (!is.null(facet_order)) {
    sd2$Facet <- factor(sd2$Facet, levels = cap(facet_order))
  } else {
    sd2$Facet <- factor(sd2$Facet, levels = unique(sd2$Facet))
  }
  
  column_labels <- as.character(sd2$Group)
  col_split     <- sd2$Facet
  
  # Parse year and season
  parse_year <- function(x){
    suppressWarnings(as.integer(stringr::str_extract(x, "(?<!\\d)\\d{4}(?!\\d)")))
  }
  parse_season <- function(x){
    s <- tolower(x)
    dplyr::case_when(
      grepl("winter", s) ~ "winter",
      grepl("spring", s) ~ "spring",
      grepl("summer", s) ~ "summer",
      grepl("fall|autumn", s) ~ "fall",
      TRUE ~ NA_character_
    )
  }
  # to allow for chronological ordering in heatmap, number rank the seasons
  season_rank <- c(winter=1, spring=2, summer=3, fall=4)
  grp_chr <- as.character(sd2$Group)
  yr <- parse_year(grp_chr)
  ss <- parse_season(grp_chr)
  rk <- unname(season_rank[ss])
  
  # order columns first by facet, then by year, then by season
  ord_cols <- order(sd2$Facet, yr, rk, na.last = TRUE)
  sd2 <- sd2[ord_cols, , drop = FALSE]
  mat <- mat[, ord_cols, drop = FALSE]
  column_labels <- as.character(sd2$Group)
  col_split     <- sd2$Facet
  
  # create ComplexHeatmap plot (white → red scale of increasing metric) 
  lo <- 0
  hi <- max(mat, na.rm = TRUE); if (hi == 0) hi <- 1
  col_fun <- circlize::colorRamp2(c(lo, hi), c("white","red"))
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = if (percent) "Rel. Abund. (%)" else "Rel. Abund.",
    col  = col_fun,
    cluster_rows    = FALSE,
    show_row_dend   = FALSE,
    cluster_columns = FALSE,
    column_split    = col_split,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    column_labels     = column_labels,
    column_names_rot  = column_rotate,
    row_names_side    = "left",
    column_names_side = "bottom",
    width  = grid::unit(ncol(mat) * 7,  "mm"),
    height = grid::unit(nrow(mat) * 13, "mm"),
    column_names_gp = grid::gpar(fontsize = col_fontsize + 3),
    row_names_gp    = grid::gpar(fontsize = row_fontsize + 5),
    gap = grid::unit(2, "mm")
  )
  
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  if (!is.null(title)) {
    grid::grid.text(
      title,
      y = grid::unit(1, "npc") - grid::unit(3.5, "mm"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    )
  }
  invisible(ht)
}

bacteria_aldex_genus_heatmap_CH(site1, compare_times = c("2018-fall","2021-fall"), title = "Open Ref Bacteria Garden Site 1 Aldex Significantly Different Taxa", sort_by = "abs_effect", top_n = 10)
bacteria_aldex_genus_heatmap_CH(site1_control, compare_times = c("2018-fall","2021-fall"), title = "Open Ref Bacteria Control Site 1 Aldex Significantly Different Taxa", sort_by = "abs_effect", top_n = 10)

bacteria_aldex_genus_heatmap_CH(ginseng_site1_garden, compare_times = c("2018-fall","2021-fall"), title = "Fungi Garden Site 1 Aldex Significantly Different Taxa", sort_by = "abs_effect", top_n = 10)
bacteria_aldex_genus_heatmap_CH(ginseng_site1_control, compare_times = c("2018-fall","2021-fall"), title = "Fungi Control Site 1 Aldex Significantly Different Taxa", sort_by = "abs_effect", top_n = 10)

#####
genus_heatmap_pathogens_CH <- function(
    physeq,
    group_var        = "seasonYear",     # column in sample_data used as column labels
    facet_var        = "Treatment",      # column to facet (split) columns by
    pathogens,                            # character vector of genus names of interest
    tax_level        = "Genus",
    title            = NULL,
    percent          = FALSE,            # show matrix in 0-1 (FALSE) or 0-100 (TRUE)
    facet_order      = NULL,             # optional: character vector to order facets
    column_rotate    = 45,               # angle for column labels
    row_fontsize     = 9,
    col_fontsize     = 9
) {
  
  # convert taxonomy vlaues to character strings and replace missing values
  tax <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  tax[] <- lapply(tax, as.character)
  tax[is.na(tax)] <- ""
  tax[tax == "__"] <- ""
  
  # identify the taxonomic ranks and fix if needed (debugging)
  fill_cols <- intersect(c("Phylum","Class","Order","Family","Genus","Species"), colnames(tax))
  for (i in seq_len(nrow(tax))) {
    new <- NULL
    if      ("Phylum" %in% colnames(tax) && tax[i, "Phylum"] == "")  new <- paste0("Unidentified_", tax[i, "Kingdom"])
    else if ("Class"  %in% colnames(tax) && tax[i, "Class"]  == "")  new <- paste0("Unidentified_", tax[i, "Phylum"])
    else if ("Order"  %in% colnames(tax) && tax[i, "Order"]  == "")  new <- paste0("Unidentified_", tax[i, "Class"])
    else if ("Family" %in% colnames(tax) && tax[i, "Family"] == "")  new <- paste0("Unidentified_", tax[i, "Order"])
    else if ("Genus"  %in% colnames(tax) && tax[i, "Genus"]  == "")  new <- paste0("Unidentified_", tax[i, "Family"])
    if (!is.null(new)) tax[i, fill_cols] <- new
  }
  rownames(tax) <- taxa_names(physeq)
  tax_table(physeq) <- as.matrix(tax)
  
  # extract sample metadata and again confirm that the grouping variable exists
  sd0 <- as.data.frame(sample_data(physeq), stringsAsFactors = FALSE)
  if (!all(c(group_var, facet_var) %in% colnames(sd0))) {
    stop("`group_var` and/or `facet_var` not found in sample_data(physeq).")
  }
  
  # create unique grouping keys to pool replicates
  comb_lab <- paste(sd0[[group_var]], sd0[[facet_var]], sep = " | ")
  comb_df <- data.frame(
    .comb  = comb_lab,
    Group  = sd0[[group_var]],
    Facet  = sd0[[facet_var]],
    row.names   = rownames(sd0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # aggregate samples by the combination (mean), then relative abundance
  phy_agg <- merge_samples(physeq, group = comb_df$.comb, fun = mean)
  
  merged_keys <- unique(comb_df[, c(".comb","Group","Facet")])
  rownames(merged_keys) <- merged_keys$.comb
  merged_keys$.comb <- NULL
  merged_keys[] <- lapply(merged_keys, function(x) ifelse(is.na(x), "Unknown", as.character(x)))
  sample_data(phy_agg) <- phyloseq::sample_data(merged_keys)
  
  phy_agg_ra <- transform_sample_counts(phy_agg, function(x) x / sum(x))
  
  # collapse to OTUs belonging to same Genus levels
  otu_mat <- as(otu_table(phy_agg_ra), "matrix")
  if (!taxa_are_rows(phy_agg_ra)) otu_mat <- t(otu_mat)
  genus_vec <- as.character(tax_table(phy_agg_ra)[, tax_level])
  genus_vec[genus_vec == "" | is.na(genus_vec)] <- "Unclassified"
  genus_rel <- rowsum(otu_mat, group = genus_vec)
  
  # subset to requested pathogen genera; drop all-zero rows
  found_pathogens <- intersect(rownames(genus_rel), pathogens)
  if (!length(found_pathogens)) stop("None of the specified pathogens were found.")
  mat <- genus_rel[found_pathogens, , drop = FALSE]
  keep_rows <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  if (!nrow(mat)) stop("All selected pathogen genera have zero abundance.")
  
  if (percent) mat <- mat * 100
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  # sort rows by overall mean abundance (DESCENDING) 
  row_means <- rowMeans(mat, na.rm = TRUE)
  ord_rows <- order(-row_means, rownames(mat))   # descending by mean, then name asc
  mat <- mat[ord_rows, , drop = FALSE]
  
  # Column labels, chronological ordering & split
  sd2 <- as.data.frame(sample_data(phy_agg_ra), stringsAsFactors = FALSE)
  sd2 <- sd2[colnames(mat), , drop = FALSE]   # align to matrix columns
  
  # capitalize facet labels for titles
  cap <- function(x) tools::toTitleCase(tolower(x))
  sd2$Facet <- cap(sd2$Facet)
  
  if (!is.null(facet_order)) {
    facet_order <- cap(facet_order)
    sd2$Facet <- factor(sd2$Facet, levels = facet_order)
  } else {
    sd2$Facet <- factor(sd2$Facet, levels = unique(sd2$Facet))
  }
  
  # Parse year and season from group labels (debug: robust to "2018-fall" or "Fall 2018")
  parse_year <- function(x) {
    y <- suppressWarnings(as.integer(stringr::str_extract(x, "(?<!\\d)\\d{4}(?!\\d)")))
    y
  }
  parse_season <- function(x) {
    s <- tolower(x)
    s <- dplyr::case_when(
      grepl("winter", s) ~ "winter",
      grepl("spring", s) ~ "spring",
      grepl("summer", s) ~ "summer",
      grepl("fall|autumn", s) ~ "fall",
      TRUE ~ NA_character_
    )
    s
  }
  # same as previous function
  season_rank <- c(winter = 1, spring = 2, summer = 3, fall = 4)
  grp_chr <- as.character(sd2$Group)
  yr  <- parse_year(grp_chr)
  ss  <- parse_season(grp_chr)
  rk  <- unname(season_rank[ss])
  # For cases like "2018-fall" that don't include the word, infer by suffix:
  missing_rk <- is.na(rk)
  if (any(missing_rk)) {
    ss2 <- tolower(gsub(".*[-_ ]", "", grp_chr))  # try last token
    rk2 <- dplyr::recode(ss2, winter = 1, spring = 2, summer = 3, fall = 4, .default = NA_real_)
    rk[missing_rk & !is.na(rk2)] <- rk2[missing_rk & !is.na(rk2)]
  }
  
  # order columns by facet, then year, then season rank
  ord_cols <- order(sd2$Facet, yr, rk, na.last = TRUE)
  sd2 <- sd2[ord_cols, , drop = FALSE]
  mat <- mat[, ord_cols, drop = FALSE]
  
  column_labels <- as.character(sd2$Group)
  col_split     <- sd2$Facet
  
  # Colors & heatmap specification (same white -> red)
  lo <- 0
  hi <- max(mat, na.rm = TRUE); if (hi == 0) hi <- 1
  col_fun <- circlize::colorRamp2(c(lo, hi), c("white", "red"))
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = if (percent) "Rel. Abund. (%)" else "Rel. Abund.",
    col = col_fun,
    
    # No row clustering/dendrogram
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    row_dend_width = grid::unit(0, "mm"),
    
    # Columns: no clustering; split by facet
    cluster_columns = FALSE,
    column_split = col_split,
    
    # Labels & layout
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_labels = column_labels,
    column_names_rot = column_rotate,
    row_names_side = "left",
    column_names_side = "bottom",
    
    # Size tweaks (wider columns, shorter rows)
    width  = grid::unit(ncol(mat) * 7,  "mm"),
    height = grid::unit(nrow(mat) * 13,  "mm"),
    
    # Readability
    column_names_gp = grid::gpar(fontsize = col_fontsize + 3),
    row_names_gp    = grid::gpar(fontsize = row_fontsize + 5),
    
    # Reduce gap between facet panels (optional)
    gap = grid::unit(2, "mm")
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  if (!is.null(title)) {
    grid::grid.text(title, y = grid::unit(1, "npc") - grid::unit(3.5, "mm"),
                    gp = grid::gpar(fontsize = 18, fontface = "bold"))
  }
  invisible(ht)
}

## fungal pathogens of interest. Import "fungalpathogens" object from "PathogensofInterest.R"
genus_heatmap_pathogens_CH(ginseng_site1, group_var = "seasonYear", pathogens = fungalpathogens, facet_var = "treatment", title = "Fungal Pathogens of Interest, Site 1")
genus_heatmap_pathogens_CH(ginseng_site2, group_var = "seasonYear", pathogens = fungalpathogens, facet_var = "treatment", title = "Fungal Pathogens of Interest, Site 2")
genus_heatmap_pathogens_CH(ginseng_site3, group_var = "seasonYear", pathogens = fungalpathogens, facet_var = "treatment", title = "Fungal Pathogens of Interest, Site 3")

## bacteria pathogens of interest
species <- c(
  "Bacillus megaterium", "Bacillus subtilis", "Pseudomonas putida", 
  "Pseudomonas arsenicooxydans", "Paenibacillus urinalis", 
  "Pseudomonas koreensis", "Brevibacterium frigoritolerans", 
  "Pseudomonas graminis", "Arthrobacter pascens", "Bacillus aryabhattai", 
  "Brevundimonas vesicularis", "Pantoea agglomerans", 
  "Paenarthrobacter nitroguajacolicus", "Microbacterium trichothecenolyticum", 
  "Microbacterium hydrocarbonoxydans"
)

# Extract unique genera
unique_genera <- unique(sub(" .*", "", species))

genus_heatmap_pathogens_CH(site1, group_var = "seasonYear",pathogens = unique_genera, facet_var = "treatment", title = "Immunosuppresive Bacteria of Interest, Site 1")
genus_heatmap_pathogens_CH(site2, group_var = "seasonYear",pathogens = unique_genera, facet_var = "treatment", title = "Immunosuppresive Bacteria of Interest, Site 2")
genus_heatmap_pathogens_CH(site3, group_var = "seasonYear",pathogens = unique_genera, facet_var = "treatment", title = "Immunosuppresive Bacteria of Interest, Site 3")