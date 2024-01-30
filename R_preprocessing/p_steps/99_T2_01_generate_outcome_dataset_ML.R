library(data.table)

events_ALL_OUTCOMES <- list()
for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"))
  temp_events_ALL_OUTCOMES <- append(events_ALL_OUTCOMES, list(get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))))
  rm(list = paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
}

events_ALL_OUTCOMES <- data.table::rbindlist(temp_events_ALL_OUTCOMES)
rm(temp_events_ALL_OUTCOMES)

events_ALL_OUTCOMES <- unique(events_ALL_OUTCOMES[, .(person_id, date, type_outcome)])
events_ALL_OUTCOMES[, flag := 1]
events_ALL_OUTCOMES <- data.table::dcast(events_ALL_OUTCOMES,
                                         person_id + date ~ type_outcome, value.var = "flag", fill = 0)

events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, c("person_id", "date") := NULL]

cols_to_remove <- c("V_CHILBLAIN_AESI", "V_MICROANGIO_AESI", "V_RENOVASCULAR_AESI", "N_TRIGEMINALNEURALGIA_AESI",
                    "O_DEATHSUDDEN_AESI", "M_OSTEOMYELITIS_AESI", "Im_SJOGRENS_AESI", "D_LIVERACUTE_AESI",
                    "C_PERICARD_AESI", "D_DIVERTICULITIS_AESI", "D_LIVERCIRRHOSIS_AESI", "Im_ANAPHYLAXISNCO_AESI",
                    "Im_ANAPHYLAXIS_AESI", "Im_KAWASAKI_AESI", "N_CONVULSION_AESI", "N_CVST_AESI", "N_MENINGOENC_AESI",
                    "N_MISKAW_AESI", "R_ARDS_AESI", "SO_CONJUNCTIVITIS_AESI")
events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, (cols_to_remove) := NULL]

test <- events_ALL_OUTCOMES[, lapply(.SD, sum)]
test_1 <- events_ALL_OUTCOMES[, rowSums(.SD) == 1]
test_1 <- events_ALL_OUTCOMES[events_ALL_OUTCOMES[, rowSums(.SD) == 1], ]
test_1 <- test_1[, lapply(.SD, sum)]
colnames(test)

events_ALL_OUTCOMES <- events_ALL_OUTCOMES[events_ALL_OUTCOMES[, rowSums(.SD) >= 1], ]
events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, count := seq(.N)]
events_ALL_OUTCOMES <- split(events_ALL_OUTCOMES, by = "count", keep.by = FALSE)

events_ALL_OUTCOMES <- sapply(events_ALL_OUTCOMES, as.matrix, simplify = FALSE, USE.NAMES = FALSE)
events_ALL_OUTCOMES[] <- lapply(events_ALL_OUTCOMES, function(x) { attributes(x) <- NULL; x })
names(events_ALL_OUTCOMES) <- NULL

np <- reticulate::import("numpy")
events_ALL_OUTCOMES<- reticulate::r_to_py(events_ALL_OUTCOMES)
np$save(paste0(dirML, "D3_events_ALL_OUTCOMES_ML") , events_ALL_OUTCOMES)
