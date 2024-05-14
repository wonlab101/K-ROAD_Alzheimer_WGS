## change log10PVAL to PVAL for manhattan plot
# by CMY

args = commandArgs(trailingOnly=TRUE)

fn=args[1]
output=args[2]

mt <- read.table(fn, header=TRUE,fill=TRUE)

mt$LOG10P <- as.numeric(as.character(mt$LOG10P))
mt$PVAL <- 10**(-mt$LOG10P)


mt$bonfP_all <- p.adjust(mt$PVAL, "bonferroni") 
mt$fdr_all <- p.adjust(mt$PVAL, "fdr") 

##per category
# Split the table based on ALLELE1 category
split_tables <- split(mt, mt$ALLELE1)

# Initialize empty data frames to store results
bonf_results <- data.frame()
fdr_results <- data.frame()

# Initialize lists to store row and column counts
row_counts <- list()
col_counts <- list()

# Define the total number of tests
total_tests <- 17174

# Perform Bonferroni and FDR correction for each category and store results
for (category in names(split_tables)) {
  sub_table <- split_tables[[category]]

  # Count rows and columns
  row_counts[[category]] <- nrow(sub_table)
  col_counts[[category]] <- ncol(sub_table)
  
  # Bonferroni correction
  sub_table$bonfP_cat <- p.adjust(sub_table$PVAL, "bonferroni", n = total_tests)
  bonf_results <- rbind(bonf_results, sub_table)
  
  # FDR correction
  sub_table$fdr_cat <- p.adjust(sub_table$PVAL, "fdr")
  fdr_results <- rbind(fdr_results, sub_table)
}

# Merge the results back together based on ALLELE1 category
mt_bonf_fdr <- merge(bonf_results, fdr_results, by = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "PVAL", "bonfP_all", "fdr_all"))


# Print row and column counts for each split table
for (category in names(split_tables)) {
  cat("Category:", category, "\n")
  cat("Number of Rows:", row_counts[[category]], "\n")
  cat("Number of Columns:", col_counts[[category]], "\n\n")
}




write.table(mt_bonf_fdr, output, row.names = FALSE, sep = '\t', quote=FALSE) 
