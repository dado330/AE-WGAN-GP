##%######################################################%##
#                                                          #
####  CREATE TIME DEPENDENT DATASETS FOR ALL VARIABLE   ####
#### NEEDED IN COHORT ANALYSIS EXCEPT NUMBER_CONDITIONS ####
#                                                          #
##%######################################################%##


# pick the coresponding conceptset(s), make some cleaning (restriction to study population, restriction to study period, restiction to appopopriate meanings if necessary), make each record last 90 days (medications) or 365 days (diagnosis), merge all the time periods (using CreateSpells), store

# use variable_definition[[element]] to retrieve all the conceptsets asscoiated to one of the elements; see 3_30 for a similar thing

for (subpop in subpopulations_non_empty) {
  
  print(subpop)
  
  # Import the study population
  name_D4_study_population <- paste0("D4_study_population", suffix[[subpop]])
  load(paste0(diroutput, name_D4_study_population, ".RData"))
  study_population <- get(name_D4_study_population)[, .(person_id, study_entry_date, study_exit_date)]
  rm(list = name_D4_study_population)
  
  # Cycle for all TD datasets
  for (COV_COHORT in c(list_of_covariates_for_cohort, list_of_comedication_for_cohort,
                       list_of_prior_history_for_cohort, list_of_pregnancy_variables)) {
    print(paste("NOW CALCULATING THE COVARIATE:", COV_COHORT))
    
    # ELEMENTS are the variable that creates the covariates
    ELEMENTS <- elements_for_TD_variables[[COV_COHORT]]
    
    # SUBELEMENTS are the conceptsets of the ELEMENTS
    SUBELEMENTS <- unlist(lapply(ELEMENTS, function(ELEMENT) variable_definition[[ELEMENT]]))
    
    # Keep only the available conceptsets
    SUBELEMENTS <- SUBELEMENTS[SUBELEMENTS %in% sub('\\.RData$', '', list.files(dirconceptsets))]
    
    # Empty array with path of saved intermediate files
    SUBELEMENTS_path <- c()
    
    for (SUBELEMENT in SUBELEMENTS) {
      
      print(paste("Covariate", COV_COHORT, "element", SUBELEMENT))
      
      # Load conceptset
      subelement_path <- paste0(dirconceptsets, SUBELEMENT, ".RData")
      CONCEPT <- rbindlist(lapply(list(subelement_path),
                                  function(x) get(load(x)[[1]])))
      
      # delete records that are not observed in this whole subpopulation
      if (this_datasource_has_subpopulations){
        CONCEPT <- CONCEPT[(eval(parse(text = select_in_subpopulationsEVENTS[[subpop]]))), ]
      }
      
      # Create dataset with study_entry_date and date of the event
      CONCEPT <- unique(CONCEPT[, .(person_id, date)])
      CONCEPT <- merge(CONCEPT, study_population[, .(person_id, study_entry_date)], by = "person_id")
      
      # Keep only events with non missing date and after the lookback time wrt the study_entry_date
      lookback <- data.table::fcase(
        SUBELEMENT == "DP_ANTIBIO" && COV_COHORT == "DP_ANTIBIO_30", 365,
        SUBELEMENT == "DP_ANTIBIO" && COV_COHORT == "DP_ANTIBIO_14", 14,
        SUBELEMENT %in% c("DP_ANTIVIR", "DP_ANTIBIO"), 14,
        SUBELEMENT == "DP_INFLUENZAVAC" && COV_COHORT == "DP_INFLUENZAVAC_365", 365,
        SUBELEMENT == "DP_INFLUENZAVAC" && COV_COHORT == "DP_INFLUENZAVAC_90", 90,
        SUBELEMENT %in% DP_variables, 90,
        SUBELEMENT %in% recurrent_OUTCOME_variables, 30,
        default = 365
      )
      
      CONCEPT <- CONCEPT[!is.na(date) & date >= study_entry_date - lookback, ]
      
      # If CONCEPT empty skip to the next
      if (nrow(CONCEPT) == 0) next
      
      # remove study entry and create interval wrt the date of events using the lookback period 
      CONCEPT[, study_entry_date := NULL][, date_end := date + days(lookback)]
      
      # Keep events until study_end
      CONCEPT <- CONCEPT[!is.na(date) & date <= study_end, ]
      
      # CONCEPT_end <- copy(CONCEPT)[, date := date + days(lookback + 1)]
      # 
      # CONCEPT <- rbindlist(list(CONCEPT[, type := "start"], CONCEPT_end[, type := "end"]))
      # rm(CONCEPT_end)
      
      # Collapse overlapping intervals
      CONCEPT <- CreateSpells(
        dataset = CONCEPT,
        id = "person_id" ,
        start_date = "date",
        end_date = "date_end",
        quiet = T
      )
      CONCEPT[, num_spell := NULL]
      
      # Save intermediate datasets
      name_export_df <- paste0("temp", "~", COV_COHORT, "~", SUBELEMENT)
      assign(name_export_df, CONCEPT)
      SUBELEMENTS_path <- append(SUBELEMENTS_path, name_export_df)
      
      save(name_export_df, file = paste0(dirTD, "/", name_export_df, ".RData"), list = name_export_df)
      
      rm(CONCEPT)
      rm(list = name_export_df)
    }
    
    # Add exception for pregnancy
    if (COV_COHORT %in% list_of_pregnancy_variables & !skip_pregnancy) {
      SUBELEMENT <- "D3_pregnancy_final"
      print(paste("Covariate", COV_COHORT, "element", SUBELEMENT))
      
      # Load conceptset
      subelement_path <- paste0(dirpregnancy, SUBELEMENT, ".RData")
      CONCEPT <- rbindlist(lapply(list(subelement_path),
                                  function(x) get(load(x)[[1]])))
      
      # Create dataset with study_entry_date and date of the event
      CONCEPT <- CONCEPT[, .(person_id, pregnancy_start_date, pregnancy_end_date)]
      CONCEPT <- merge(CONCEPT, study_population[, .(person_id, study_entry_date)], by = "person_id")
      
      # Keep only events with non missing dates and pregnancy_end_date after the study_entry_date
      CONCEPT <- CONCEPT[!is.na(pregnancy_start_date) & !is.na(pregnancy_end_date) & pregnancy_end_date >= study_entry_date, ]
      
      # remove study entry
      CONCEPT[, study_entry_date := NULL]
      
      # Keep events until study_end
      CONCEPT <- CONCEPT[pregnancy_start_date <= study_end, ]
      
      # If CONCEPT empty skip to the next
      if (nrow(CONCEPT) == 0) next
      
      setnames(CONCEPT, c("pregnancy_start_date", "pregnancy_end_date"), c("entry_spell_category", "exit_spell_category"))
      
      # Save intermediate datasets
      name_export_df <- paste0("temp", "~", COV_COHORT, "~", SUBELEMENT)
      assign(name_export_df, CONCEPT)
      SUBELEMENTS_path <- append(SUBELEMENTS_path, name_export_df)
      
      save(name_export_df, file = paste0(dirTD, "/", name_export_df, ".RData"), list = name_export_df)
      
      rm(CONCEPT)
      rm(list = name_export_df)
    }
    
    # If all conceptsets are empty or with unusable events the save an empty datasets
    if (is.null(SUBELEMENTS_path)) {
      tmp <- study_population[, .(person_id, study_entry_date)][, value_of_variable := 0]
      setnames(tmp, "study_entry_date", "date")
      name_export_df <- paste0("D3_TD_variable_", COV_COHORT, suffix[[subpop]])
      assign(name_export_df, tmp)
      save(name_export_df, file = paste0(dirTD, "/", name_export_df, ".RData"), list = name_export_df)
      
      rm(tmp)
      rm(list = name_export_df)
      
      next
    }
    
    final_COVARIATE <- list()
    
    # Load intermediates dataset and append them
    for (single_file in SUBELEMENTS_path) {
      load(file = paste0(dirTD, "/", single_file, ".RData"))
      final_COVARIATE <- append(final_COVARIATE, list(get(single_file)))
      file.remove(paste0(dirTD, "/", single_file, ".RData"))
      objects_to_remove <- c(single_file)
      rm(list = objects_to_remove)
    }
    final_COVARIATE <- rbindlist(final_COVARIATE, fill = T)
    
    # Collapse overlapping intervals between conceptset this time
    final_COVARIATE <- CreateSpells(
      dataset = final_COVARIATE,
      id = "person_id" ,
      start_date = "entry_spell_category",
      end_date = "exit_spell_category",
      quiet = T
    )
    final_COVARIATE[, num_spell := NULL]
    
    # Copy the dataset and exit_spell_category add 1 day to get the first day without events
    final_COVARIATE_end <- copy(final_COVARIATE)[, entry_spell_category := NULL]
    final_COVARIATE_end[, exit_spell_category := exit_spell_category + 1]
    setnames(final_COVARIATE_end, "exit_spell_category", "date")
    
    # In dataset remove exit_spell_category
    final_COVARIATE[, exit_spell_category := NULL]
    setnames(final_COVARIATE, "entry_spell_category", "date")
    
    # Add the population in both dataset because we need to remove all events happening after the study exit
    final_COVARIATE <- merge(final_COVARIATE, study_population[, .(person_id, study_exit_date)], by = "person_id")
    final_COVARIATE <- final_COVARIATE[date <= study_exit_date, ]
    final_COVARIATE[, study_exit_date := NULL]
    
    final_COVARIATE_end <- merge(final_COVARIATE_end, study_population[, .(person_id, study_exit_date)], by = "person_id")
    final_COVARIATE_end <- final_COVARIATE_end[date <= study_exit_date, ]
    final_COVARIATE_end[, study_exit_date := NULL]
    
    # Create dataset which contains only persons excluded by the prebious datasets
    to_add_rows <- merge(study_population[, .(person_id, study_entry_date)], final_COVARIATE, all.x = T, by = "person_id")
    to_add_rows[, flag_add_row := fifelse(is.na(date) | date > study_entry_date, T, F)]
    
    to_add_rows <- to_add_rows[, .(study_entry_date, flag_add_row = min(flag_add_row)), by = "person_id"]
    to_add_rows <- to_add_rows[flag_add_row == 1, ][, flag_add_row := NULL]
    setnames(to_add_rows, "study_entry_date", "date")
    
    # Events happening before the study entry are set at the study entry
    final_COVARIATE <- merge(final_COVARIATE, study_population[, .(person_id, study_entry_date)], by = "person_id")
    final_COVARIATE[date < study_entry_date, date := study_entry_date]
    final_COVARIATE[, study_entry_date := NULL]
    
    # Combine all datasets. value_of_variable is 1 when person has event.
    final_COVARIATE <- rbindlist(list(final_COVARIATE[, value_of_variable := 1],
                                      final_COVARIATE_end[, value_of_variable := 0],
                                      to_add_rows[, value_of_variable := 0]))
    rm(final_COVARIATE_end, to_add_rows)
    
    # Export final dataset
    name_export_df <- paste0("D3_TD_variable_", COV_COHORT, suffix[[subpop]])
    assign(name_export_df, final_COVARIATE)
    save(name_export_df, file = paste0(dirTD, "/", name_export_df, ".RData"), list = name_export_df)
    
    rm(final_COVARIATE)
    rm(list = name_export_df)
  }
}

