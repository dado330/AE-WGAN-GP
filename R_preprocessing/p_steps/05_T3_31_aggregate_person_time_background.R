##%######################################################%##
#                                                          #
####  AGGREGATE DATASETS TO ADD TOTALS FOR SOME STRATA  ####
#                                                          #
##%######################################################%##


for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(diroutput, "D4_persontime_background", suffix[[subpop]], ".RData"))
  PT_background <- get(paste0("D4_persontime_background", suffix[[subpop]]))
  rm(list = paste0("D4_persontime_background", suffix[[subpop]]))
  
  # Find the columns for counts and PT
  cols_to_sums <- names(PT_background)[grepl("^Persontime|_b$", names(PT_background))]
  
  # Sums by sex to get data for both sexes together
  PT_background_sex <- PT_background[, lapply(.SD, sum, na.rm=TRUE),
                                     by = c("year", "Ageband", "COVID19"),
                                     .SDcols = cols_to_sums][, sex := "total"]
  
  PT_background <- rbindlist(list(PT_background, PT_background_sex), use.names = T)
  
  # Sums by Ageband to get data for all ages
  PT_background_age <- PT_background[, lapply(.SD, sum, na.rm=TRUE),
                                     by = c("sex", "year", "COVID19"),
                                     .SDcols = cols_to_sums][, Ageband := "total"]
  
  PT_background <- rbindlist(list(PT_background, PT_background_age), use.names = T)
  
  # Sums by year 2019/2020 for table 6
  PT_background_year <- PT_background[year %in% c("2019", "2020"), ][, lapply(.SD, sum, na.rm=TRUE),
                                                               by = c("sex", "Ageband", "COVID19"),
                                                               .SDcols = cols_to_sums][, year := "2019/2020"]
  
  PT_background <- rbindlist(list(PT_background, PT_background_year), use.names = T)
  
  nameoutput <- paste0("D4_persontime_background_aggregated", suffix[[subpop]])
  assign(nameoutput, PT_background)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
}



