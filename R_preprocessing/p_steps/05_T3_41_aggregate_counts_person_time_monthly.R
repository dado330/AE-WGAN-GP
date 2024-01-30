##%######################################################%##
#                                                          #
####          AGGREGATE DATASETS TO ADD TOTALS          ####
####            FOR SOME STRATA AND CLEANING            ####
#                                                          #
##%######################################################%##

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(diroutput, "D4_counts_persontime_monthly", suffix[[subpop]], ".RData"))
  PT_monthly <- get(paste0("D4_counts_persontime_monthly", suffix[[subpop]]))
  rm(list = paste0("D4_counts_persontime_monthly", suffix[[subpop]]))
  
  # Find the columns for counts and PT
  cols_to_sums <- names(PT_monthly)[grepl("^Persontime|_b$", names(PT_monthly))]
  
  # Sums by sex to get data for both sexes together
  PT_monthly_sex <- PT_monthly[, lapply(.SD, sum, na.rm=TRUE),
                               by = c("month", "Ageband", "COVID19"),
                               .SDcols = cols_to_sums][, sex := "total"]
  
  PT_monthly <- rbindlist(list(PT_monthly, PT_monthly_sex), use.names = T)
  
  # Sums by Ageband to get data for all ages
  PT_monthly_age <- PT_monthly[, lapply(.SD, sum, na.rm=TRUE),
                               by = c("sex", "month", "COVID19"),
                               .SDcols = cols_to_sums][, Ageband := "total"]
  
  PT_monthly <- rbindlist(list(PT_monthly, PT_monthly_age), use.names = T)
  
  PT_monthly[, c("year", "month") := tstrsplit(month, "-")]
  PT_monthly[, month := month.name[as.integer(month)]]
  
  colstokeep <- c('sex','year')
  colstokeep <- c(colstokeep, setdiff(colnames(PT_monthly), colstokeep))
  PT_monthly <- PT_monthly[,..colstokeep]
  
  nameoutput <- paste0("D4_counts_persontime_monthly_aggregated", suffix[[subpop]])
  assign(nameoutput, PT_monthly)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
}



