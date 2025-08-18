library(readxl) 
library(phyloseq)

pathogen_df <- read_excel("Fungal_Plant_Pathogens_2025.xlsx",
                          sheet = "Fungal Pathogenic Genera")

head(pathogen_df)

fungalpathogens <- pathogen_df[["Fungal Plant Pathogen Genus"]]

# 3) Extract the taxa names from phyloseq object
#    Assumptions: 
#    tax_table row names (taxa_names) correspond to ASV/OTU IDs,
#    want to check “Genus” column in the tax_table.

taxa_tbl <- as.data.frame(tax_table(fungi_only))  
phy_genus <- taxa_tbl$Genus

# 4) Which pathogens are in your phyloseq dataset?
present  <- intersect(pathogens, phy_genus)
missing  <- setdiff(pathogens,   phy_genus)

# 5) Report
cat("Present in phyloseq:\n")
print(present)

cat("\nMissing from phyloseq:\n")
print(missing)



