# Install the below packages if you have not:
# install.packages(c("parallel", "foreach", "doParallel", "purrr", "dplyr", "tidyr", "expss", "filesstrings"))


# Load packages
lapply(c("parallel", "foreach", "doParallel", "purrr", "dplyr", "tidyr", "expss", "filesstrings"), 
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
if (is_empty(files)==FALSE){
  registerDoParallel(numCores)
  foreach (file = files) %dopar% {
    csv <- read.csv(file)
    file.move(file, "./cache/a_csv_files", overwrite=TRUE)
    csv$Sample <- rep(file,nrow(csv))
    csv$Sample <- gsub('.PASSED.VCs.tsv.csv', '', csv$Sample)
    colnames(csv) <- c("Category", "Variant_Caller", "hg19_Chr_Pos", "Ref", "Alt", 
                       "COSM_ref", "Gene", "NM#:c_Change", "p_Change", "VAF_percentage", "Alt_Fwd,Alt_Rev/Ref_Fwd,Ref_Rev", 
                       "PoN_Mean_SD", "Qscore", "Strand_bias_test", "Sample")
    csv <- separate(csv,  "Alt_Fwd,Alt_Rev/Ref_Fwd,Ref_Rev",
                    into=c("Alt_Fwd,Alt_Rev", "Ref_Fwd,Ref_Rev"), sep="/")
    csv <- separate(csv, "Alt_Fwd,Alt_Rev",
                    into=c("Alt_Fwd", "Alt_Rev"), sep=",")
    csv <- separate(csv, "Ref_Fwd,Ref_Rev",
                    into=c("Ref_Fwd", "Ref_Rev"), sep=",")
    csv <- separate(csv, "Sample",
                    into=c("Sample#", "Timepoint", "Indexes"), sep="_")
    csv <- mutate(csv, vaf_diff=abs(VAF_percentage - PoN_Mean_SD))
    csv <- filter(csv, vaf_diff >= 0.01)
    csv <- filter(csv, Qscore >= 20)
    csv <- filter(csv, grepl('PASS', Strand_bias_test))
    csv$Variant_hg19 <- paste0(csv$hg19_Chr_Pos, "-", csv$Ref, "-", csv$Alt)
    csv <- csv[c("Category", "Sample#", "Timepoint", "Variant_Caller", "Variant_hg19",  
                 "COSM_ref", "Gene", "NM#:c_Change", "p_Change", "VAF_percentage", "Alt_Fwd", "Alt_Rev", "Ref_Fwd", "Ref_Rev", 
                 "PoN_Mean_SD", "vaf_diff", "Qscore", "Strand_bias_test")]
    filename <- paste0(file, ".sample.csv")
    write.table(csv, file=filename, sep=",", row.names = FALSE)
  }
}


# Combine each sample
files <- list.files(pattern='*.sample.csv')
if (is_empty(files)==FALSE){
  csv <- lapply(files, read.csv)
  csv <- do.call(rbind.data.frame, csv)
  file.move(files, "./cache/b_filtered_csv", overwrite=TRUE)
  colnames(csv) <- c("Category", "Sample", "Timepoint", "Variant_Caller", "Variant_hg19", 
                     "COSM_ref", "Gene", "NM_c_Change", "p_Change", "VAF_percentage", "Alt_Fwd", "Alt_Rev", "Ref_Fwd", "Ref_Rev", 
                     "PoN_Mean_3SD", "vaf_diff", "Qscore", "Strand_bias_test")
  samples <- unique(csv[,"Sample"])
  registerDoParallel(numCores)
  foreach (sample = samples) %dopar% {
    csv1 <- filter(csv, Sample == sample)
    filename <- paste0(sample, ".comb.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
  }
}


# Seek for duplicate and unique variants
files <- list.files(pattern='*.comb.csv')
if (is_empty(files)==FALSE){
  registerDoParallel(numCores)
  foreach (file = files) %dopar% {
    csv <- read.csv(file)
    file.move(file, "./cache/c_comb_csv", overwrite=TRUE)
    csv1 <- filter(csv, Timepoint == 0)
    colnames(csv1) <- c("initial_Category", "Sample", "initial_Timepoint", "initial_Variant_Caller", "Variant_hg19", 
                        "COSM_ref", "Gene", "NM_c_Change", "p_Change", "initial_VAF_percentage", 
                        "initial_Alt_Fwd", "initial_Alt_Rev", "initial_Ref_Fwd", "initial_Ref_Rev", "initial_PoN_Mean_3SD", 
                        "initial_vaf_diff", "initial_Qscore", "initial_Strand_bias_test")
    csv2 <- filter(csv, Timepoint != 0)
    colnames(csv2) <- c("followup_Category", "followup_Sample", "followup_Timepoint", "followup_Variant_Caller", "Variant_hg19", 
                        "followup_COSM_ref", "followup_Gene", "followup_NM_c_Change", "followup_p_Change", "followup_VAF_percentage", 
                        "followup_Alt_Fwd", "followup_Alt_Rev", "followup_Ref_Fwd", "followup_Ref_Rev", "followup_PoN_Mean_3SD", 
                        "followup_vaf_diff", "followup_Qscore", "followup_Strand_bias_test")
    csv3 <- add_columns(csv1, csv2, by="Variant_hg19")
    csv4 <- filter(csv3, !is.na(followup_Sample))
    csv5 <- filter(csv3, is.na(followup_Sample))
    csv4 <- csv4[c("Variant_hg19", "Sample", "COSM_ref", "Gene", "NM_c_Change", "p_Change", 
                  "initial_Timepoint", "initial_Category", "initial_Variant_Caller", "initial_VAF_percentage", "initial_Alt_Fwd", "initial_Alt_Rev", "initial_Ref_Fwd", 
                  "initial_Ref_Rev", "initial_PoN_Mean_3SD", "initial_vaf_diff", 
                  "followup_Timepoint", "followup_Category", "followup_Variant_Caller", "followup_VAF_percentage", "followup_Alt_Fwd", "followup_Alt_Rev", "followup_Ref_Fwd", 
                  "followup_Ref_Rev", "followup_PoN_Mean_3SD", "followup_vaf_diff")]
    csv5 <- csv5[c("Variant_hg19", "Sample", "COSM_ref", "Gene", "NM_c_Change", "p_Change", 
                   "initial_Timepoint", "initial_Category", "initial_Variant_Caller", "initial_VAF_percentage", "initial_Alt_Fwd", "initial_Alt_Rev", "initial_Ref_Fwd", 
                   "initial_Ref_Rev", "initial_PoN_Mean_3SD", "initial_vaf_diff")]
    filename1 <- paste0(file, ".dup.csv")
    filename2 <- paste0(file, ".ini.csv")
    write.table(csv4, file=filename1, sep=",", row.names = FALSE)
    write.table(csv5, file=filename2, sep=",", row.names = FALSE)
  }
}


# Export and cleanup
if (dir.exists("# Export")==FALSE){
  dir.create("# Export")
}
files <- list.files(pattern='*.dup.csv')
if (is_empty(files)==FALSE){
  file.move(files, "./# Export/a_initial_and_followup", overwrite=TRUE)
}
files <- list.files(pattern='*.ini.csv')
if (is_empty(files)==FALSE){
  file.move(files, "./# Export/b_initial_only", overwrite=TRUE)
}
rm("csv", "files", "numCores", "samples")
print("All done! Please check the Export folder.")

