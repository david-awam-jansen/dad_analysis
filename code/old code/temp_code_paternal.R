## library

library(broom)
library(multidplyr)
library(purrr)
library(tidyverse)

## variables

zero_daily_count <- 1/365.25
log_zero_daily_count <- log2(zero_daily_count)

## data
dsi_raw <- readRDS('../R01_grant/data/dsi_pop_1_days_20191204.RDS')
biograph_l <- readRDS('../R01_grant/data/biograph_l_20191204.RDS')

dads <- #tbl(babase, 'parents') %>%  collect() %>%
	parents_l %>%
	select(sname = kid, dad)

moms <- #tbl(babase, 'parents') %>%  collect() %>%
	parents_l %>%
	select(sname = kid, mom)

adult_members <- tbl(babase, "members") %>%
	select(grp, date, sname) %>%
	filter(grp < 3) %>%
	inner_join(select(tbl(babase, "biograph"), sname, sex)) %>%
	left_join(select(tbl(babase, "maturedates"), sname, matured)) %>%
	left_join(select(tbl(babase, "rankdates"), sname, ranked)) %>%
	mutate(is_adult = case_when(sex == "F" & matured <= date ~ TRUE,
															sex == "M" & ranked <= date ~ TRUE)) %>%
	filter(is_adult == TRUE) %>%
	collect()


## functions
get_F_dyads <- function(df) {
	df %>%
		filter(dyad_type == "AF-JUV")
}

add_paternal <- function(focal, focal_grp, focal_start, df) {
	## change name of grp column and make the juvenile the partner in all dyads
	## This is to create 'fake' dryads later
	df <- df %>%  rename(sname_grp = grp)
	df <- df %>%
		filter(sname_sex_class == "JUV") %>%
		rename_at(vars(starts_with("sname")), ~str_replace(., "sname", "sname_temp")) %>%
		rename_at(vars(starts_with("partner")), ~str_replace(., "partner", "sname")) %>%
		rename_at(vars(starts_with("sname_temp")), ~str_replace(., "sname_temp", "partner")) %>%
		bind_rows(df %>%
								filter(sname_sex_class != "JUV"))


	df <- df %>%
		left_join(select(moms, sname, sname_mom = mom), by = "sname") %>%
		left_join(select(moms, partner = sname, partner_mom = mom), by = "partner") %>%
		left_join(select(dads, sname, sname_dad = dad), by = "sname") %>%
		left_join(select(dads, partner = sname, partner_dad = dad), by = "partner") %>%
		mutate(paternal_groom = case_when(
			partner == sname_mom & sname_sex_class == "AF" ~ "maternal_check",
			sname == partner_mom  & sname_sex_class == "AF" ~ "maternal",
			sname != partner_mom  & sname_sex_class == "AF" ~ "non_maternal",
			partner == sname_dad & sname_sex_class == "AM" ~ "paternal_check",
			sname == partner_dad  & sname_sex_class == "AM" ~ "paternal",
			sname != partner_dad  & sname_sex_class == "AM" ~ "non_paternal",
			TRUE ~ "other")) %>%
		filter(!str_detect(paternal_groom, "other"))
}

make_NA <- function(x) NA

add_missing_paternal <- function(focal, focal_grp, focal_start, df) {
	## Make sure that every juvenile will have a paternal and not paternal value
	## with at least 1 groom

	mean_AF_log2OE <- df %>%
		filter(dyad_type == "AF-JUV") %>%
		group_by(partner, partner_grp) %>%
		summarise(mean_AF_log2OE = log2(mean(OE)), .groups = "drop")

	df %>%
		bind_rows(add_all_options) %>%
		filter(i_adj > 0) %>%
		ungroup() %>%
		complete(paternal_groom,
						 nesting(partner_grp, sname_grp),
						 fill=list(partner = 'XXX',
						 					sname = 'XXX')) %>%
		group_by(partner_grp, sname_grp) %>%
		complete(partner, nesting(paternal_groom),
						 fill = list(sname = 'XXX'
						 						, i_adj = zero_daily_count
						 						, log2_i_adj = log2(zero_daily_count)
						 						, n_focals = .1
						 )) %>%
		mutate(dyad_type = if_else(str_detect(paternal_groom, "maternal"),
															 "AF-JUV",
															 "AM-JUV")) %>%
		# #bind_rows(add_all_opttions) %>%
		# filter(i_adj > 0) %>%
		# group_by(dyad_type) %>%
		# complete(paternal_groom, nesting(dyad_type, partner, partner_grp, sname_grp)) %>%
		# group_by(partner_grp, sname_grp) %>%
		# complete(partner, nesting(paternal_groom, dyad_type),
		#          fill = list(sname = 'XXX'
		#                      , i_adj = zero_daily_count
		#                      , log2_i_adj = log2(zero_daily_count)
		#                      , n_focals = .1
		#                      #, focal = focal
	#                      #, focal_grp = focal_grp
	#          )) %>%
	left_join(mean_AF_log2OE, by = c("partner", "partner_grp")) %>%
		select(-partner_dad, -partner_mom) %>%
		left_join(select(moms, partner = sname, partner_mom = mom), by = "partner") %>%
		left_join(select(dads, partner = sname, partner_dad = dad), by = "partner") %>%
		filter((str_detect(paternal_groom, "maternal") & !is.na(partner_mom)) |
					 	(str_detect(paternal_groom, "paternal") & !is.na(partner_dad))) %>%
		mutate(log2OE = if_else(sname == 'XXX', mean_AF_log2OE, log2OE))  %>%
		mutate(focal_id = paste(focal, focal_grp, focal_start, sep="-")) %>%
		mutate(focal_check = case_when(
			(sname == focal | partner == focal) & partner_grp == focal_grp ~ TRUE,
			TRUE ~ FALSE))
}

focal_filter <- function(df) {
	df %>%
		filter(focal_check == TRUE)
}

get_sub_regression <- function(df) {
	if(unique(df$dyad_type) == "AF-JUV") {
		df <- df %>%
			filter(!is.na(partner_mom))
	} else if (unique(df$dyad_type) == "AM-JUV") {
		df <- df %>%
			filter(!is.na(partner_dad))
	}

	if(all(is.na(df$i_adj))|all(df$i_adj == 0)) {
		NULL
	} else {
		df %>%
			filter(df$i_adj > 0 & n_focals > 0) %>%
			lm(formula = log2_i_adj ~ log2OE, na.action=na.exclude)
	}
}

zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

get_paternal_universal_zscored <- function(df) {

	sv <- universal_paternal_values
	names(sv) <- c("dyad_type", "B0", "B1")
	df <- df %>%inner_join(sv, by = "dyad_type")

	df <- df %>%
		mutate(i_adj = if_else(i_adj == 0, zero_daily_count, i_adj),
					 log2_i_adj = log2(i_adj)) %>%
		filter(i_adj > 0) %>%
		group_by(paternal_groom) %>%
		mutate(res_value =
					 	log2_i_adj - (B0 + (B1 * log2OE))) %>%
		ungroup() %>%
		filter(res_value > -9999) %>%
		group_by(dyad_type) %>%
		mutate_at(vars(starts_with('res_value')),.funs = list(zscore = ~zscore(.))) %>%
		select(focal_check, sname, partner, dyad_type, paternal_groom, log2OE, i_adj, log2_i_adj, res_value, zscore)

	df
}

create_paternal_universal_zscored <- function(df) {

	df_overall <- df %>%
		filter(focal_check == TRUE) %>%
		group_by(dyad_type) %>%
		arrange(-res_value) %>%
		slice(1:3) %>%
		summarise(DSI = mean(zscore)) %>%
		dplyr::mutate(dyad_type = case_when(
			dyad_type == "AF-JUV" ~ "DSI_F",
			dyad_type == "AM-JUV" ~ "DSI_M",
			TRUE ~ "check")) %>%
		pivot_wider(names_from = dyad_type, values_from = DSI)

	df_paternal <- df %>%
		filter(focal_check == TRUE) %>%
		group_by(paternal_groom) %>%
		arrange(-res_value) %>%
		slice(1:3) %>%
		summarise(DSI = mean(zscore)) %>%
		dplyr::mutate(paternal_groom = case_when(
			paternal_groom == "maternal" ~ "DSI_maternal",
			paternal_groom == "non_maternal" ~ "DSI_Fme",
			paternal_groom == "paternal" ~ "DSI_paternal",
			paternal_groom == "non_paternal" ~ "DSI_Mde",
			TRUE ~ "check")) %>%
		pivot_wider(names_from = paternal_groom, values_from = DSI)

	df_overall %>%
		bind_cols(df_paternal)
}

get_all_adults <- function(selected_sex, focal_grp, focal_start, focal_end) {
	adult_members %>%
		filter(sex == selected_sex &
					 	grp == focal_grp &
					 	date > focal_start &
					 	date < focal_end) %>%
		select(sname) %>%
		distinct() %>%
		pull()
}

get_potential_partners <- function(selected_sex, focal_grp, df) {
	if(selected_sex == "F") {
		selected_dyad = "AF-JUV"
	} else if(selected_sex == "M") {
		selected_dyad = "AM-JUV"
	} else {
		print("A valid sex needs to be selected")
	}

	df %>%
		ungroup() %>%
		filter(dyad_type == selected_dyad) %>%
		filter(partner_grp == focal_grp) %>%
		filter(sname != 'XXX') %>%
		select(sname) %>%
		distinct() %>%
		pull()
}

get_partners <- function(selected_sex, focal, focal_grp, df) {
	if(selected_sex == "F") {
		selected_dyad = "AF-JUV"
	} else if(selected_sex == "M") {
		selected_dyad = "AM-JUV"
	} else {
		print("A valid sex needs to be selected")
	}

	df %>%
		ungroup() %>%
		filter(dyad_type == selected_dyad) %>%
		filter(focal_check) %>%
		#select(partner, sname, paternal_groom, contains("i_adj"))
		filter(i_adj > 0) %>%
		filter(sname != 'XXX') %>%
		select(sname) %>%
		distinct() %>%
		pull()
}

check_parent <- function(selected_parent, df) {
	selected_parent %in% df
}

check_which <- purrr::possibly(which, otherwise = 99)

get_paternal_position <- function(row, df) {
	#print(row)

	df %>%
		filter(focal_check) %>%
		filter(zscore > -999) %>%
		group_by(dyad_type) %>%
		arrange(dyad_type, -zscore) %>%
		select(focal_check, partner, sname, dyad_type, paternal_groom,zscore) %>%
		#mutate(position = percent_rank(zscore)) %>%
		mutate(position = cume_dist(zscore)) %>%
		filter(paternal_groom %in% c("maternal", "paternal")) %>%
		ungroup() %>%
		mutate(dyad_type = if_else(dyad_type == "AF-JUV", "mom_position", "dad_position")) %>% select(dyad_type, position, position) %>%
		pivot_wider(names_from = dyad_type, values_from = position)
}


## first steps
dsi <- dsi_raw %>%
	mutate(nr_dyads = map_dbl(.x = subset, nrow)) %>%
	filter(nr_dyads > 0) %>%
	mutate(F_dyads = map(.x = subset, .f = get_F_dyads)) %>%
	mutate(nr_F_dyads = map_dbl(.x = F_dyads, nrow)) %>%
	filter(nr_F_dyads > 0) %>%
	select(-di) %>%
	left_join(dads, by = "sname")

## We should exclude the additional lines that were added for maternal SCI

dsi <- dsi %>%
	#filter(sname %in%  xdata_females_with_social$sname) %>%
	inner_join(select(biograph_l, sname, birth)) %>%
	mutate(age = as.numeric(round(((start - birth)/365.25)))) %>%
	filter(age < 4)

dsi_final <- dsi %>%
	mutate(paternal_dyads = pmap(.l = list(sname, grp, start, subset),
															 .f = add_paternal))

add_all_options <- dsi_final %>%
	slice(1) %>%
	select(paternal_dyads) %>%
	unnest(cols= c(paternal_dyads)) %>%
	slice(1:4) %>%
	mutate_all(, .funs = make_NA) %>%
	mutate(dyad_type = rep(c("AF-JUV", "AM-JUV"), each = 2)) %>%
	mutate(paternal_groom = c("maternal",
														"non_maternal",
														"paternal",
														"non_paternal")) %>%
	mutate(i_adj = zero_daily_count)

dsi_final <- dsi_final %>%
	mutate(paternal_data = pmap(.l = list(sname, grp, start, paternal_dyads),
															.f = add_missing_paternal)) %>%
	mutate(focal_data = map(.x = paternal_data,
													.f = focal_filter))

count_dyads <- function(df){
	df %>%
		filter(focal_check) %>%
		group_by(paternal_groom) %>%
		summarise(cases = n()) %>%
		pivot_wider(names_from = paternal_groom, values_from = cases)
}

dsi_final %>%
	mutate(check = map(.x = paternal_data, count_dyads)) -> temp

dsi_final_full <- dsi_final

save(dsi_final, dads, moms, file = "dsi_final_full_30JUN20.RData")

my_results_only_focal <- dsi_final %>%
	select(sex_class, sex, focal_data) %>%
	unnest(cols = c(focal_data))

universal_paternal_values <-  my_results_only_focal %>%
	mutate(group = dyad_type) %>%
	filter(sname != 'XXX') %>%
	group_by(group) %>%
	nest()  %>%
	mutate(regression = purrr::map(data,get_sub_regression)) %>%
	mutate(model_results = regression %>%  map(.,.f = broom::tidy)) %>%
	unnest(model_results) %>%
	select(dyad_type = group,  term, estimate) %>%
	pivot_wider(names_from = term, values_from = estimate)

dsi_final <- dsi_final %>%
	select(-subset, -F_dyads, -focal_data)

save(dsi_final, universal_paternal_values, dads, moms, file = "dsi_final_6FEB21.RData")
rm(list=ls(pattern="dsi"))

load("dsi_final_6FEB21.RData")

my_results_step1 <- dsi_final %>%
	ungroup() %>%
	group_by(sname, start, end, grp) %>%
	mutate(zscored_resvalues =
				 	map(.x = paternal_data,
				 			.f = get_paternal_universal_zscored))

write_rds(my_results_step1, "./data/my_results_step1.rds")

my_results_step2 <- my_results_step1 %>%
	mutate(zscored_DSI =
				 	map(.x = zscored_resvalues,
				 			.f = create_paternal_universal_zscored))

write_csv(my_results_step2, "./data/my_results_step2.csv")

my_results_zscored <- my_results_step2 %>%
	unnest(keep_empty = TRUE, cols = c(zscored_DSI))

my_results_zscored <- my_results_zscored %>%
	inner_join(select(tbl(babase, "parents"), sname = kid, mom) %>%  collect()) %>%
	ungroup() %>%
	group_by(sname, grp, start, end) %>%
	mutate(all_AF = pmap(.l = list("F", grp, start, end), .f = get_all_adults)) %>%
	mutate(all_AM = pmap(.l = list("M", grp, start, end), .f = get_all_adults)) %>%
	mutate(potential_females = pmap(.l = list("F", grp, paternal_dyads),
																	.f = get_potential_partners)) %>%
	mutate(potential_males = pmap(.l = list("M", grp, paternal_dyads),
																.f = get_potential_partners)) %>%
	mutate(AF_partners = pmap(.l = list("F", sname, grp, zscored_resvalues),
														.f = get_partners)) %>%
	mutate(AM_partners = pmap(.l = list("M", sname, grp, zscored_resvalues),
														.f = get_partners)) %>%
	mutate(mom_present = map2_lgl(.x = mom, .y = all_AF, .f=check_parent)) %>%
	mutate(mom_groomed = map2_lgl(.x = mom, .y = AF_partners, .f=check_parent)) %>%
	mutate(dad_present = map2_lgl(.x = dad, .y = all_AM, .f=check_parent)) %>%
	mutate(dad_groomed = map2_lgl(.x = dad, .y = AM_partners, .f=check_parent)) %>%
	mutate(nr_AM = map_int(.x = potential_males, n_distinct),
				 nr_AM_partners = map_int(.x = AM_partners, n_distinct),
				 groom_type = case_when(
				 	dad_present == TRUE & dad_groomed == TRUE ~ "groomed by dad",
				 	dad_present == TRUE & dad_groomed == FALSE ~ "not groomed by dad that is present",
				 	dad_present == FALSE & dad_groomed == FALSE ~ "dad not present"))

my_results_zscored <- my_results_zscored %>%
	ungroup() %>%
	bind_cols(tibble(rownumber = seq(from = 1, to = dim(my_results_zscored)[1]))) %>%
	mutate(paternal_positions = map2(.x = rownumber, .y = zscored_resvalues,
																	 .f = get_paternal_position)) %>%
	unnest(cols = c(paternal_positions), keep_empty = TRUE)

save(my_results_zscored, file = "my_results_zscored_6FEB21.RData")
load("my_results_zscored_30JUL20.RData")

df_temp <- my_results_zscored %>%
	select(sname, grp, start, end, days_present, sex, birth,
				 starts_with("DSI_")) %>%
	write_rds("./data/DSI_paternal_values_new_version.RDS")

DSI_paternal_values_new <- read_rds("./data/DSI_paternal_values_new_version.RDS")
DSI_paternal_values %>%
	pivot_longer()


xdata_ea %>%
	inner_join(dad) %>%
	filter(sex == "Females" &
				 	included_cases == TRUE &
				 	!is.na(dad) &
				 	statage > 4) %>%
	select(sname) %>%
	pull() -> selected_females

selected_data <- my_results_zscored %>%
	filter(sname %in% selected_females)

rm("my_results_zscored")

## To show  to Beth
selected_data %>%
	ggplot(aes(dad_present, fill = dad_groomed)) +
	geom_bar() +
	facet_wrap(age ~ .)

selected_data %>%
	ggplot(aes(nr_AM_partners, fill = dad_groomed)) +
	geom_bar() +
	facet_wrap(age ~ .) +
	theme(legend.position = 'bottom')

selected_data %>%
	filter(age == 0 &
				 	nr_AM_partners &
				 	dad_groomed == FALSE) %>%
	slice(1) %>%
	inner_join(dad %>%  select(sname, focal_dad = dad)) %>%
	select(focal = sname, focal_dad, zscored_resvalues) %>%
	unnest(cols = c(zscored_resvalues)) %>%
	filter(focal_check) %>%
	mutate(explain = case_when(
		str_detect(paternal_groom, "paternal") & sname == 'XXX' ~ 'made up paternal dyad',
		str_detect(paternal_groom, "non_paternal") & sname != 'XXX' ~ 'who are these males',
		TRUE ~ ''))







