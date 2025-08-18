#### SCRIPT ON GENERATING STATISTICS FOR ALPHA AND BETA DIVERSITY
#### BY MICHAEL YOSHIDA, FOR THE TAI LAB, 2025


#### ALPHA DIVERSITY

library(phyloseq)
library(ggplot2)
library(multcomp)

ps <- garden_rareF18F22
group_var <- "year"  

# 3) Estimate alpha diversity measures
alpha_df <- estimate_richness(ps, measures =  "Shannon")

# 4) Bring in metadata
meta_df  <- data.frame(sample_data(ps))
# Make sure the rownames match
alpha_df <- cbind(alpha_df, meta_df[rownames(alpha_df), , drop = FALSE])

# 5) Convert grouping column to factor
alpha_df[[group_var]] <- factor(alpha_df[[group_var]])
alpha_df$site          <- factor(alpha_df$site)         # make site a factor too

# 6) ANOVA: e.g. Shannon diversity by group
shannon_aov <- aov(Shannon ~ year * site, data = alpha_df)
summary(shannon_aov)

TukeyHSD(shannon_aov, which = "year:site")



#### BETA DIVERSITY
library(vegan)
library(ggplot2)

ps <- site1_rare

bc_dist <- distance(ps, method = "bray")

# 2) Extract sample metadata
meta <- data.frame(sample_data(ps))
meta$site <- factor(meta$site)

set.seed(123)
adonis_res <- adonis2(bc_dist ~ year + treatment, 
                      data = meta,
                      permutations = 999, by = "margin")
print(adonis_res)

set.seed(123)
adonis_res <- adonis2(bc_dist ~ year + treatment, 
                      data = meta,
                      permutations = 999, by = "margin")
print(adonis_res)

# spread (homogenous vs heterogeneous)
set.seed(123)
disp <- betadisper(bc_dist, meta$collection_year)
disp_test <- permutest(disp, pairwise = TRUE, permutations = 999)
print(disp_test)

# exploring a plotting option
dist_df <- data.frame(
  Year = meta$collection_year,
  Dist = disp$distances
)
ggplot(dist_df, aes(x = Year, y = Dist)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Multivariate Dispersion (distance to centroid) by Year",
    y = "Distance to group centroid"
  )



