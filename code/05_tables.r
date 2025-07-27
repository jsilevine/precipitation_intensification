##---------------------------------------------------------------
## 05_tables.r -- generate summary tables for model outputs
##---------------------------------------------------------------

library(brms)

mod_files <- list.files("data/model_objects")

for (f in mod_files) {

  mod <- readRDS(paste0("data/model_objects/", f))
  summary_table <- summary(mod)$fixed
  summary_table <- summary_table[,c(1,3:ncol(summary_table))]

  write.csv(round(summary_table, 3), paste0("tables/", gsub(".rds", ".csv", f)))

}
