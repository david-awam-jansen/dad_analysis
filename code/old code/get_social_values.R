Sys.setenv(TZ = 'UTC')
list.of.packages <- list("foreach", "doSNOW", "parallel", "tidyverse",
												 "lubridate", "dbplyr", "purrrlyr",# "RPostgreSQL",
												 "zoo")
new.packages <- list.of.packages[!(list.of.packages %in%
																	 	installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages,function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_social_index_analysis_", latest_version_date,".RData"))
parents_l <- read_csv("./data/parents_l.csv")

ea_dataset_less <- read_csv(file = paste("./data/", "ea_dataset_less_restricted_", latest_version_date, ".csv", sep=""))

dyadic_index_messages <- tibble(row_nr = iyol_sub$row_nr, message = NA)

source('./code/social_indexes/sociality-indices.R')

iyol_sub <- iyol_sub %>%
	inner_join(select(ea_dataset_less, sname)) %>%
	filter(days_present > 59) %>%
	filter(age_start_yrs < 4) %>%
	arrange(sname, start, grp) %>%
	mutate(row_nr = row_number())

load(paste0("./data/data_set_for_social_index_analysis_", latest_version_date,".RData"))
write_csv(ea_dataset_less,"./data/ea_dataset_less.csv")
write_csv(iyol_sub,"./data/iyol_sub.csv")

# sci_data <- sci(my_iyol = iyol_sub , biograph_l = biograph_l,
#  		members_l = members_l, focals_l = focals_l, females_l = females_l,
#  		interactions_l = grooming_l,
#  		min_res_days = 60, parallel = TRUE, ncores = 48,
#  		directional = FALSE)
#
# write_rds(x = sci_data,
# 					file = paste0("./data/sci_data", latest_version_date, ".rds"))

dsi_data <- dyadic_index(my_iyol = iyol_sub[1:4, ], biograph_l = biograph_l,
		members_l = members_l, focals_l = focals_l, females_l = females_l,
		interactions_l = grooming_l,
		min_cores_days = 60, parallel = TRUE, ncores = 4, directional = FALSE)

write_rds(dsi_data, paste0("./data/dsi_data_", latest_version_date, ".rds"))




