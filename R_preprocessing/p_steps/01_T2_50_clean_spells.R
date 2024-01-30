##%######################################################%##
#                                                          #
####        CLEAN THE SPELLS AND SELECT THE ONE         ####
####         USED IN THE STUDY FOR EACH PERSON          ####
#                                                          #
##%######################################################%##


load(paste0(dirtemp,"D3_PERSONS.RData"))

D3_clean_spells <- list()
for (subpop in subpopulations_non_empty) {
  
  load(paste0(dirtemp,"D3_output_spells_category.RData"))
  
  if (this_datasource_has_subpopulations == TRUE){
    person_spell <- merge(D3_output_spells_category[[subpop]], D3_PERSONS, all.x = T, by = "person_id")
  } else {
    person_spell <- merge(D3_output_spells_category, D3_PERSONS, all.x = T, by = "person_id")
  }
  
  person_spell <- person_spell[, .(person_id, birth_date, death_date, entry_spell_category_crude = entry_spell_category,
                                   exit_spell_category_crude = exit_spell_category, op_meaning, num_spell)]

  person_spell[, entry_spell_category := data.table::fifelse(birth_date < entry_spell_category_crude - 60 | birth_date > instance_creation,
                                                             entry_spell_category_crude,
                                                             birth_date)]
  person_spell[, exit_spell_category := pmin(exit_spell_category_crude, death_date, na.rm = T)]

  person_spell[, op_start_date_cleaned := data.table::fifelse(entry_spell_category != entry_spell_category_crude, 0, 1)]
  person_spell[, op_end_date_cleaned := data.table::fifelse(exit_spell_category != exit_spell_category_crude, 0, 1)]
  person_spell[, starts_at_birth := data.table::fifelse(entry_spell_category == birth_date, 1, 0)]
  person_spell[, starts_after_ending := data.table::fifelse(exit_spell_category < entry_spell_category, 1, 0)]
  
  # find spells that do not overlap the study period
  person_spell[, no_overlap_study_period := fifelse(
    entry_spell_category > study_end | exit_spell_category < study_start, 1, 0)]
  
  # find spells that are shorter than 365 days
  person_spell[, less_than_365_days_and_not_starts_at_birth := fifelse(
    correct_difftime(pmin(exit_spell_category, study_end), entry_spell_category) <= 365 & starts_at_birth == 0, 1, 0)]
  
  # Import doses dataset and create doses criteria
  load(paste0(dirtemp,"D3_vaccines_curated.RData"))
  
  # Create important index dates, select only variables of interest and merge with spells
  setorderv(D3_vaccines_curated, c("person_id", "date_curated"))
  D3_vaccines_curated[, first_vax_date := date_curated[1], by = person_id]
  D3_vaccines_curated <- unique(D3_vaccines_curated[, .(person_id, first_vax_date)])
  person_spell <- merge(person_spell, D3_vaccines_curated, all.x = T, by = "person_id")
  
  # find spells which contains the first vaccination
  person_spell[, spell_without_vax1 := fifelse(
    entry_spell_category > first_vax_date | exit_spell_category < first_vax_date, 1, 0, na = 1)]
  
  # find if spells with the first vaccination has a distance between first vax and entry spell of less or equal than 365
  person_spell[, has_vax1_before_365_days := fifelse(
    spell_without_vax1 == 0 & entry_spell_category >= first_vax_date - 365, 1, 0, na = 0)]
  person_spell[, first_vax_date := NULL]
  
  ## find the study spell of each person
  # select rows not excluded by criteria starts_after_ending, no_overlap_study_period,
  # less_than_365_days_and_not_starts_at_birth and has_vax1_before_365_days
  person_spell[starts_after_ending == 0 & no_overlap_study_period == 0 &
                 less_than_365_days_and_not_starts_at_birth == 0 & has_vax1_before_365_days == 0,
               flag := 0]
  
  # if a person still have a spell with vax1 then it becomes the study spell
  person_spell[flag == 0 & spell_without_vax1 == 0, is_the_study_spell := 1]
  person_spell[person_id %in% person_spell[is_the_study_spell == 1, person_id], flag := 1]
  
  # On the spells with still flag equal to 0 take the one with the minimum exit_spell_category after study_start
  person_spell[flag == 0 & exit_spell_category >= study_start,
               min_exit_spell_category := min(exit_spell_category), by = person_id]
  
  person_spell[exit_spell_category == min_exit_spell_category, is_the_study_spell := 1]
  person_spell[, c("min_exit_spell_category", "flag") := NULL]
  
  person_spell[is.na(is_the_study_spell), is_the_study_spell := 0]
  
  if (this_datasource_has_subpopulations == TRUE){
    D3_clean_spells[[subpop]] <- person_spell
  } else {
    D3_clean_spells <- person_spell
  }
}

save(D3_clean_spells, file = paste0(dirtemp, "D3_clean_spells.RData"))
rm(person_spell, D3_clean_spells)
