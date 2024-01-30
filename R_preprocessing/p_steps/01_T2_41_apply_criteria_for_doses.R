##%######################################################%##
#                                                          #
####  APPLY EXCLUSION CRITERIA FOR VACCINATION RECORDS  ####
#                                                          #
##%######################################################%##


# input: D3_clean_vaccines
# output: Flowchart_QC_criteria, D3_vaccines_curated

load(paste0(dirtemp, "D3_clean_vaccines.RData"))

# Crate the flowchart and filter the record of the doses
D3_vaccines_curated <- CreateFlowChart(
  dataset = D3_clean_vaccines,
  listcriteria = c("duplicated_records", "manufacturer_not_in_study", "missing_date", "date_before_start_vax",
                   "distance_btw_1_2_doses", "distance_btw_2_3_doses", "distance_btw_3_4_doses", "dose_after_4"),
  flowchartname = "Flowchart_criteria_for_doses")

# Save the flowchart
fwrite(Flowchart_criteria_for_doses, paste0(direxpsubpop[[subpop]], "Flowchart_criteria_for_doses.csv"))

# Clean and save the final dataset for thye doses
D3_vaccines_curated <- D3_vaccines_curated[, .(person_id, date_curated, dose_curated, manufacturer_curated)]
save(D3_vaccines_curated, file = paste0(dirtemp, "D3_vaccines_curated.RData"))

rm(D3_vaccines_curated, D3_clean_vaccines, Flowchart_criteria_for_doses)