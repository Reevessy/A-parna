# Install the below packages if you have not:
# install.packages(c("parallel", "foreach", "doParallel", "purrr", "dplyr", "tidyr", "filesstrings"))


# Load packages
lapply(c("parallel", "foreach", "doParallel", "purrr", "dplyr", "tidyr", "filesstrings"), 
       require, character.only = TRUE)


# Set the folder containing current R script as working directory
setwd(".")
print(paste("Current working dir:", getwd()))
if (dir.exists("original_files")==FALSE){
  dir.create("original_files")
}
if (dir.exists("cache")==FALSE){
  dir.create("cache")
}


# Set cpu cores for parallel computing
numCores <- detectCores(all.tests = FALSE, logical = TRUE)
print(paste("Parallel computing:", numCores, "cores will be used for data processing"))


# Read raw tsv files
files <- list.files(pattern='*.VCs.tsv')
if (is_empty(files)==FALSE){
  registerDoParallel(numCores)
  foreach (file = files) %dopar% {
    tsv <- read.table(file, fill=TRUE, sep="\t", na.strings=c("","NA"))
    file.move(file, "./original_files/a_tsv_files", overwrite=TRUE)
    filename <- paste0(file, ".csv")
    write.table(tsv, file=filename, sep=",", row.names = FALSE, col.names = FALSE)
  }
}


# Trim and filter the data
files <- list.files(pattern='*.tsv.csv')
registerDoParallel(numCores)
foreach (file = files) %dopar% {
  csv <- read.csv(file)
  file.move(file, "./cache/a_csv_files", overwrite=TRUE)
  csv$Sample <- rep(file,nrow(csv))
  csv$Sample <- gsub('.PASSED.VCs.tsv.csv', '', csv$Sample)
  colnames(csv) <- c("Category", "Variant_Caller", "hg19_chr:hg19_pos", "Ref_Allele", "Alt_Allele", 
                     "COSM_ref", "Gene", "NM#:c_Change", "p_Change", "VAF_percentage", "Alt_Fwd,Alt_Rev/Ref_Fwd,Ref_Rev", 
                     "PoN_Mean_SD", "Qscore", "Strand_bias_test", "Sample")
  csv <- mutate(csv, vaf_diff=abs(VAF_percentage - PoN_Mean_SD))
  csv <- filter(csv, vaf_diff >= 0.2)
  csv <- filter(csv, Qscore >= 20)
  csv <- filter(csv, grepl('PASS', Strand_bias_test))
  csv <- csv[c("Category", "Sample", "Variant_Caller", "hg19_chr:hg19_pos", "Ref_Allele", "Alt_Allele", 
               "COSM_ref", "Gene", "NM#:c_Change", "p_Change", "VAF_percentage", "Alt_Fwd,Alt_Rev/Ref_Fwd,Ref_Rev", 
               "PoN_Mean_SD", "vaf_diff", "Qscore", "Strand_bias_test")]
  filename <- paste0(file, ".sample.csv")
  write.table(csv, file=filename, sep=",", row.names = FALSE)
}


# Export and cleanup
if (dir.exists("# Export")==FALSE){
  dir.create("# Export")
}
files <- list.files(pattern='*.sample.csv')
if (is_empty(files)==FALSE){
  csv <- lapply(files, read.csv)
  csv <- do.call(rbind.data.frame, csv)
  file.move(files, "./cache/b_filtered_csv", overwrite=TRUE)
  colnames(csv) <- c("Category", "Sample", "VariantCaller", "hg19.chr:hg19.pos", "Ref.Allele", "Alt.Allele", 
                     "COSM.ref", "Gene", "NM#:c.Change", "p.Change", "%VAF", "Alt.Fwd,Alt.Rev/Ref.Fwd,Ref.Rev", 
                     "PoN.Mean+3SD", "vaf.diff", "Qscore", "Strand.bias.test")
  filename <- paste0("filtered_variants", ".csv")
  write.table(csv, file=filename, sep=",", row.names = FALSE)
}
files <- list.files(pattern='*.csv')
file.move(files, "./# Export", overwrite=TRUE)
rm("csv", "files", "numCores", "filename")
print("All done! Please check the Export folder.")

