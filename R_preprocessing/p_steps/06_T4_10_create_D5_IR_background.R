##%######################################################%##
#                                                          #
####              CALCULATE BACKGROUND IRS              ####
#                                                          #
##%######################################################%##

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  df_to_cycle <- paste0(c("D4_persontime_background_aggregated", "D4_counts_persontime_monthly_aggregated"), suffix[[subpop]])
  df_export_to_cycle <- c("D5_IR_background", "D5_IR_monthly")
  for (i in seq_along(df_to_cycle)) {
    print(df_to_cycle[[i]])
    
    load(paste0(diroutput, df_to_cycle[[i]] ,".RData"))
    persontime_windows <- get(df_to_cycle[[i]])
    rm(list = df_to_cycle[[i]])
    
    for (ev in c(OUTCOME_variables, CONTROL_variables, "DEATH")) {
      name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
      name_count <- paste0(ev,"_b")
      name_pt <- paste0("Persontime_",ev)
      persontime_windows[, (name_cols) := exactPoiCI(persontime_windows, name_count, name_pt)]
    }
    
    nameoutput <- df_export_to_cycle[[i]]
    assign(nameoutput, persontime_windows)
    save(nameoutput, file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".RData"), list = nameoutput)
    
    fwrite(get(nameoutput), file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".csv"))
    
  }
}
