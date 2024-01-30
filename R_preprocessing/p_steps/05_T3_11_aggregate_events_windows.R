##%######################################################%##
#                                                          #
####  AGGREGATE DATASETS TO ADD TOTALS FOR SOME STRATA  ####
#                                                          #
##%######################################################%##

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(diroutput, "D4_count_events_windows", suffix[[subpop]], ".RData"))
  persontime_windows <- get(paste0("D4_count_events_windows", suffix[[subpop]]))
  rm(list = paste0("D4_count_events_windows", suffix[[subpop]]))
  
  cols_to_sums <- names(persontime_windows)[grepl("^Persontime|_b$", names(persontime_windows))]
  
  persontime_windows <- persontime_windows[, lapply(.SD, sum, na.rm=TRUE),
                                           by = c("sex", "dose", "type_vax", "period", "Ageband"),
                                           .SDcols = cols_to_sums]
  
  nameoutput <- paste0("D4_count_events_windows_aggregated", suffix[[subpop]])
  assign(nameoutput, persontime_windows)
  save(nameoutput, file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".RData"), list = nameoutput)
  
  fwrite(get(nameoutput), file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".csv"))
}
