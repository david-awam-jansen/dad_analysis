source('./code/functions_for_data_prep.R')

## data
## if data needs to be uploaded source('./code/extract_babase_data_for_data_prep')
latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_dad_analysis_", latest_version_date, ".RData"))
anubis_estimates <- read_csv("./data/genome_wide_anubis_estimates_for_David_23Mar2021.csv")

## social data
dyadic_data <- read_rds("./data/selected_dad_grooming.RDS")

male_values <- dyadic_data %>%
	mutate(DSI_values = map(.x = zscored_resvalues, .f = get_DSI_values))

social_indexes_males <- male_values %>%
	inner_join(select(biograph_l, sname, birth)) %>%
	mutate(age = as.numeric((start - birth)/365.25)) %>%
	select(sname, grp, start, end, age, DSI_values) %>%
	unnest(cols = c(DSI_values))

juvenile_social_indexes_males <- social_indexes_males %>%
	filter(age < 4) %>%
	pivot_longer(names_to = "index", values_to = "value", DSI_paternal:DSI_Mde_top) %>%
	group_by(index) %>%
	group_by(age, index) %>%
	mutate(value.scaled = get_zscore(value)) %>%
	group_by(sname, index) %>%
	summarise(jmean = mean(value.scaled, na.rm = TRUE),
						nr_years = sum(!is.na(value.scaled)),
						.groups = 'drop') %>%
	filter(nr_years >= 2) %>%
	select(-nr_years) %>%
	ungroup() %>%
	pivot_wider(names_from = index, values_from = jmean, names_prefix = "j")


## Early adversity dataset
## whole juvenile period
xdata_ea_raw <- read_csv('./data/ea_dataset_less_restricted_20200302.csv')

xdata_ea_raw %>%
	filter(sname %in% selected_females)

xdata_ea_raw <- read_csv('../DSI_paternal_female_survival/data/ea_dataset_less_restricted_20200302.csv')

xdata_ea <- xdata_ea_raw %>%
	inner_join(select(parents_l, sname = kid, mom, dad)) %>%
	mutate(statage = as.numeric((statdate - birth) / 365.25)) %>%
	mutate(included_cases = complete.cases(
		sname, birth, matgrp, bstatus, statage,
		mom, dad,
		maternal_loss, maternal_SCI_F,
		maternal_rank, sibling, density, drought)) %>%
	mutate(survival_status_at_age_4 = case_when(
		ymd(birth) + years(4) < ymd(statdate) ~ 0,
		ymd(birth) + years(4) > ymd(statdate) &
			status %in% c(0, 1) ~ 1,
		ymd(birth) + years(4) > ymd(statdate) &
			status %in% c(2, 3) ~ 0),
		survival_status_at_age_4_desc = case_when(
			ymd(birth) + years(4) < ymd(statdate) ~ 0,
			ymd(birth) + years(4) > ymd(statdate) &
				status %in% c(0, 1) ~ 1,
			ymd(birth) + years(4) > ymd(statdate) &
				status %in% c(2, 3) ~ 999),
		survival_status_at_age_4_desc = factor(
			survival_status_at_age_4_desc,
			labels = c("Alive", "Died", "Censored")
		),
		age_at_age_4 = ifelse(statage <= 4, statage, 4),
		adult_survival_status = ifelse(status == 1, 1, 0),
		## create binary variables using Tung eta al cutoffs
		maternal_rank_binary = maternal_rank >= 12,
		maternal_SCI_binary = maternal_SCI_F <= as.numeric(
			quantile(xdata_ea_raw$maternal_SCI_F, na.rm = TRUE, probs = 0.25)),
		density_binary = density > 35,
		sex = factor(sex, labels = c("Females", "Males"))) %>%
	group_by(sname) %>%
	mutate(cumulative_adversity = sum(maternal_loss + maternal_rank_binary +
																			maternal_SCI_binary + sibling +  density_binary + drought)) %>%
	ungroup()

xdata_females_with_social_temp <- xdata_ea %>%
	filter(sex == 'Females')  %>%
	filter(age_at_age_4 >= 4) %>%
	inner_join(select(parents_l, sname = kid, mom, dad)) %>%
	filter(!is.na(dad) & !is.na(mom)) %>%
	inner_join(juvenile_social_indexes_males, by = "sname")

xdata_females_with_social_temp <- xdata_females_with_social_temp %>%
	mutate(end_juvenile = birth + years(4) - days(1)) %>%
	mutate(dad_overlap = pmap(.l = list(sname, dad, birth, end_juvenile),
														.f = get_overlap)) %>%
	unnest(cols = c(dad_overlap)) %>%
	rename_with(stringr::str_replace,
							pattern = "partner", replacement = "dad")

xdata_females_with_social <- xdata_females_with_social_temp %>%
	arrange(dad_overlap_years, sname) %>%
	mutate(ordered_sname = seq(1:nrow(xdata_females_with_social_temp)))


# xdata_females_with_social %>%
# select(jDSI_paternal, jDSI_Mde_top, jDSI_Mtop) %>%
# 	pivot_longer(names_to = "index", values_to = "value", jDSI_paternal:jDSI_Mtop) %>%
# 	mutate(index = case_when(index == "jDSI_Mtop" ~ "top male DSI",
# 													 index == "jDSI_Mde_top" ~ "non-dad male DSI",
# 													 index == "jDSI_paternal" ~ "DSI with dad")) %>%
# 	ggplot() +
# 	geom_density(aes(value, color = index), size = 3) +
# 	xlab("Dyadic bond strength") +
# 	cowplot::theme_cowplot() +
# 	theme(legend.position = "bottom",
# 				legend.title = element_blank()) +
# 	guides(color = guide_legend(ncol = 1))

temp_data_overlap_plot <- xdata_females_with_social_temp %>%
	mutate(class = round(dad_overlap_years)) %>%
	group_by(class) %>%
	summarise(cases =n ()) %>%
	mutate(class_max = cumsum(cases)) %>%
	mutate(percentage = round(class_max/max(class_max) * 100,0))

## annual time varying
xdata_ea_annual_temp <- read_csv("./data/social_index_values.csv") %>%
	inner_join(select(xdata_females_with_social, sname))

xdata_ea_annual <-  xdata_ea_annual_temp %>%
	mutate(dad_overlap = pmap(.l = list(sname, dad, start, end, per_group = grp),
														.f = get_overlap)) %>%
	unnest(cols = c(dad_overlap)) %>%
	rename_with(stringr::str_replace,
							pattern = "partner", replacement = "dad")


## silk figure data
get_AM_count <- function(focal, focal_grp, start, end, focal_dad) {
	start <- lubridate::ymd(start)
	end <- lubridate::ymd(end)

	members_l %>%
		inner_join(select(rankdates_l, sname, ranked), by = "sname") %>%
		inner_join(members_l %>%
							 	filter(sname == focal
							 				 & grp == focal_grp
							 				 & date >= start
							 				 & date <= end
							 	) %>%
							 	select(grp, date),
							 by = c("grp", "date")) %>%
		filter(date > ranked) %>%
		group_by(date) %>%
		summarise(nr_adult_males = n()
							, nr_adult_males_not_dad = sum(sname != focal_dad)
							, dad_presence = sum(sname == focal_dad)
							, .groups = 'drop') %>%
		summarise(nr_days = n()
							, mean_nr_adult_males = mean(nr_adult_males)
							, mean_nr_adult_males_not_dad = mean(nr_adult_males_not_dad)
							, mean_dad_presence = mean(dad_presence)
							, .groups = 'drop')
}

get_grooming_partners <- function(df) {
	df %>%
		filter(focal_check
					 & sname != 'XXX'
					 & str_detect(paternal_groom, "paternal")
		) %>%
		group_by(paternal_groom) %>%
		summarize(N = n(),
							.groups = 'drop') %>%
		pivot_wider(names_from = 'paternal_groom', values_from = 'N')
}

Silk_fig4_step1 <- male_values %>%
	select(sname, grp, start, end, dad) %>%
	mutate(AM_count = pmap(.l = list(sname, grp, start, end, dad), .f = get_AM_count)) %>%
	unnest(cols = c(AM_count))

Silk_grooming_step1 <- dyadic_data %>%
	mutate(grooming = map(.x = zscored_resvalues, .f = get_grooming_partners))

Silk_grooming_step2 <- Silk_grooming_step1 %>%
	unnest(cols = c(grooming), keep_empty = TRUE) %>%
	select(sname, age, grp, start, end, paternal, non_paternal) %>%
	mutate_at(vars(paternal:non_paternal), ~replace(., is.na(.), 0))

Silk_fig4_step2 <- Silk_fig4_step1 %>%
	inner_join(Silk_grooming_step2) %>%
	select(sname, grp, age, nr_days, contains('mean'), paternal, non_paternal) %>%
	pivot_longer(names_to = 'index', values_to = 'value', mean_nr_adult_males:non_paternal) %>%
	group_by(sname, age, index) %>%
	summarise(value = weighted.mean(x = value, w = nr_days, na.rm = TRUE), .groups='drop')

Silk_fig4_step2_wide <- Silk_fig4_step2 %>%
	pivot_wider(names_from = index, values_from = value) %>%
	mutate(AM_grooming = paternal + non_paternal)

save(Silk_fig4_step1,
		 Silk_fig4_step2,
		 Silk_fig4_step2_wide,
		 # Silk_fig4_step3,
		 # Silk_fig4_step3b,
		 # Silk_fig4_step4,
		 Silk_groomimng_step1,
		 Silk_grooming_step2,
		 file = "./data/Silk_figure_data.Rdata")

dad_top_groomer <- male_values %>%
	mutate(top_dyads = map(.x = zscored_resvalues, .f= get_top_grooms),
				 has_real_partners = map_lgl(.x = top_dyads, .f = has_real_partners),
				 nr_dyads = map_dbl(.x = top_dyads, .f = check_nr_dyads)) %>%
	mutate(is_dad_top = map_lgl(.x = top_dyads, .f = get_dad_top)) %>%
	unnest(cols = c(DSI_values)) %>%
	select(sname, grp, age, DSI_Mde_top, is_dad_top, has_real_partners, dad_present) %>%
	mutate(dad_status =
				 	case_when(dad_present == FALSE ~ "Dad not present",
				 						has_real_partners == FALSE ~ "No adult male grooming partners",
				 						is_dad_top == TRUE ~  "Dad is the top grooming partner",
				 						is_dad_top == FALSE ~  "Dad is the not top partner")) %>%
	mutate(dad_status =
				 	forcats::fct_relevel(dad_status,
				 											 "Dad is the top grooming partner", after =0L))


save(xdata_females_with_social,
		 xdata_ea_annual,
		 male_values,
		 social_indexes_males,
		 juvenile_social_indexes_males,
		 dad_top_groomer,
		 actor_actees_l,
		 file = "processed_data_for_dad_analysis.RData")

## Prepare dataset for model 1 and model2
## per year group include all ranked males present for at least 30 days

load('./data/dsi_final_full_30JUN20.RData')
get_mean_observer_effort <-function(df) {

	if(nrow(df) > 1) {
		df %>%
			ungroup() %>%
			select(mean_AF_log2OE) %>%
			distinct() %>%
			pull()
	} else {
		NA
	}
}

observer_effort_values <-dsi_final %>%
	ungroup() %>%
	filter(sname %in% xdata_females_with_social$sname) %>%
	mutate(observer_effort = map_dbl(.x = focal_data, .f = get_mean_observer_effort)) %>%
	select(focal = sname, focal_grp = grp, start, end, observer_effort)
rm("dsi_final")

who_grooms_step1 <- social_indexes_males %>%
	select(focal = sname, focal_grp = grp, start, end, contains('DSI')) %>%
	inner_join(select(xdata_females_with_social,
										focal = sname, mom, dad, maternal_loss, maternal_rank, maternal_SCI_F,
										density, sibling, drought, cumulative_adversity), by = "focal") %>%
	inner_join(select(biograph_l, focal = sname, focal_birth = birth), by = "focal") %>%
	inner_join(select(biograph_l, mom = sname, mom_birth = birth), by = "mom") %>%
	mutate(kid_age = as.numeric(start - focal_birth)/365.25) %>%
	mutate(mom_age = as.numeric(start - mom_birth)/365.25) %>%  # mom age
	mutate(nr_focal_days = pmap_dbl(
		.l = list(focal, focal_grp, start, end),
		.f = get_focal_days)) %>%
	ungroup()

who_grooms_step2 <- who_grooms_step1 %>%
	inner_join(observer_effort_values) %>%
	ungroup()

who_grooms_step3 <- who_grooms_step2 %>%
	mutate(all_males = pmap(.l = list(focal_grp, start, end), .f = get_all_males)) %>%
	unnest(cols = c(all_males)) %>%
	left_join(male_values %>%
							select(focal = sname, focal_grp = grp, start, end, zscored_resvalues) %>%
							unnest(cols = zscored_resvalues) %>%
							filter(focal_check == TRUE & str_detect(paternal_groom, "paternal")) %>%
							select(focal, focal_grp, start, sname, zscore),
						by = c("focal", "focal_grp", "start", "sname")) %>%
	rename(AMales = sname) %>%
	group_by(focal, focal_grp, start) %>%
	mutate(nr_males = n_distinct(AMales),
				 nr_grooms = sum(!is.na(zscore)),
				 does_groom = !is.na(zscore)) %>%
	group_by(AMales) %>%
	mutate(is_dad = AMales == dad) %>%
	inner_join(select(biograph_l, AMales = sname, AMales_birth = birth)) %>%
	mutate(AMales_age = as.numeric((start - AMales_birth)/365.25)) %>%
	ungroup()

who_grooms_step4 <- who_grooms_step3 %>%
	mutate(AMales_rank = pmap(.l = list(AMales, focal_grp, start, end),
														.f = get_mean_rank)) %>%
	unnest(cols = c(AMales_rank), keep_empty = TRUE)

who_grooms_step5 <- who_grooms_step4 %>%
	left_join(select(potential_dads_l, focal = kid, AMales = pdad,
									 pdad_status = status,
									 estrous_presence, estrous_me, estrous_c),
						by = c('AMales', "focal")) %>%
	mutate(estrous_presence =if_else(is.na(estrous_presence), 0, estrous_presence),
				 estrous_c =if_else(is.na(estrous_c), 0, estrous_c),
				 estrous_me =if_else(is.na(estrous_me), 0, estrous_me)) %>%
	mutate(daily_d_days_rate = pmap_dbl(.l = list(AMales, focal_grp, start, end),
																			.f = get_d_days_rate)) %>%
	mutate(nr_pdads = map_dbl(.x = focal, .f = get_nr_pdads)) %>%
	mutate(offspring_years = pmap_dbl(
		.l = list(focal, AMales, focal_grp, start, end),
		.f = get_offspring_years))

who_grooms_step5 %>%
	select(-zscore) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

who_grooms_step6 <- who_grooms_step5 %>%
	ungroup() %>%
	mutate(kids = pmap(.l = list(focal, focal_birth, mom, AMales),
										 .f = get_kids)) %>%
	unnest(cols = c(kids), keep_empty = TRUE)

who_grooms_step6 %>%
	select(-zscore) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

who_grooms_step7 <- who_grooms_step6 %>%
	left_join(select(hybridgene_scores_l, AMales = sname, hybrid_score = score)) %>%
	left_join(select(anubis_estimates, AMales = sname, anubis_admix, notes) %>%
							filter(is.na(notes) & !is.na(anubis_admix)))

who_grooms_step8 <- who_grooms_step7 %>%
	select(focal, focal_birth, focal_grp, mom, dad, AMales, start, end, nr_days
				 ## ages
				 , kid_age, mom_age, AMales_age,

				 , observer_effort
				 ## early adversity
				 , maternal_loss, maternal_rank, maternal_SCI_F, density, sibling, drought, cumulative_adversity

				 ## male details
				 , is_dad, mean_ordrank, mean_proprank,
				 , contains('estrous'), , daily_d_days_rate, offspring_years, nr_pdads
				 , previous_kid_with_mom, next_kid_with_mom ## see bootstrap

				 ## for subset
				 , hybrid_score, anubis_admix

				 ## responses
				 , does_groom,  zscore) %>%
	filter(!is.na(mean_proprank)) %>%  ## There are a few issues with male ranks
	filter(nr_days >= 30) ## Males have to be in group for at least 30 days

## Not all previous or next siblings have know dads.
## We initially assigned a .5 value to the binary variable if a malke was the dad
## We then decided to use bootstrapping to use 'randonmly' assigned paterity.
## Only males that were potential dads can bve assigned paternity

nr_draws = 1000

who_grooms_step_previous<- who_grooms_step7 %>%
	filter(previous_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "previous") %>%
	mutate(get_previous = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
														 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_previous), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(previous_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																			.f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, previous_random_draws)

who_grooms_step_next <- who_grooms_step7 %>%
	filter(next_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "next") %>%
	mutate(get_next = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
												 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_next), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(next_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																	.f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, next_random_draws)

who_grooms_step_combined <- who_grooms_step7 %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	left_join(who_grooms_step_previous) %>%
	left_join(who_grooms_step_next)

who_grooms_bootstrap <- who_grooms_step_combined %>%
	mutate(data = pmap(.l = list(data, previous_random_draws, next_random_draws),
										 .f = assign_random_draws)) %>%
	unnest(cols = c(data)) %>%
	select(-previous_random_draws, -next_random_draws) %>%
	filter(!is.na(mean_proprank)) %>%
	filter(nr_days >= 30)

save(xdata_females_with_social,
		 file = "./data/data_for_paper_3MAY21.RData")

write_csv(who_grooms_bootstrap, "./data/who_grooms_bootstrap.csv")

# move who_grooms_bootstrap.csv to cluster and run cluster_bootstrap_previous_next_code.R
# and get the results back here

result_files <- list.files(path = './data/bootstrap_results/bootstrap_models')

for(ii in 1:length(result_files)) {
	xdata <- read_delim(file = paste0('./data/bootstrap_results/bootstrap_models/', result_files[ii]),
											delim = " ", col_names = FALSE)

	if(ncol(xdata) == 10) names(xdata) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "p.value", "VIF", "model_run")
	if(ncol(xdata) == 11) names(xdata) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "df", "p.value", "VIF", "model_run")

	xdata <- xdata %>%
		filter(effect == 'fixed') %>%
		select(-effect, -group, -row, -model_run) %>%
		mutate(term = case_when(str_detect(term , "bootstrap_previous") ~ "previous_kid_with_mom",
														str_detect(term , "bootstrap_next") ~ "next_kid_with_mom",
														TRUE ~ term))

	new_filename <- str_replace(string = result_files[ii], pattern = "_raw.text", replacement = ".csv")
	write_csv(xdata, paste0("./data/bootstrap_results/", new_filename))
	print(paste0(new_filename, " has been saved"))
}


xdata_sum <- xdata %>%
	group_by(term) %>%
	select(-sig) %>%
	summarise_all(.funs = mean) %>%
	mutate(sig = p.value < 0.05)

full_model_tidy <- full_model %>%
	tidy() %>%
	filter(effect == 'fixed') %>%
	select(-effect, -group) %>%
	mutate(sig = p.value < .05)

only_know_dads_data <- who_grooms6 %>%
	filter(previous_kid_with_mom != .5 & next_kid_with_mom != .5)

only_known_dads <- update(full_model, data = only_known_dads_data)

only_known_dads_tidy <- only_known_dads %>%
	tidy() %>%
	filter(effect == 'fixed') %>%
	select(-effect, -group) %>%
	mutate(sig = p.value < .05)


full_model_tidy %>%
	mutate(model = "unk dad set at 0.5") %>%
	bind_rows(xdata %>%
							mutate(model = "individual bootstrap models")) %>%
	bind_rows(xdata_sum %>%
							mutate(model = "mean estimate 0f bootstrap")) %>%
	bind_rows(only_known_dads_tidy %>%
							mutate(model = "only known paternities"))  %>%
	ggplot() +
	geom_point(aes(x = estimate, y = term, color = sig, shape = model, size = model)) +
	scale_shape_manual(values = c(20, 19, 17, 15)) +
	scale_size_manual(values = c(1,3, 3, 3)) +
	scale_color_manual(values = c("grey50", "black")) +
	cowplot::theme_cowplot() +
	theme(legend.position = "bottom") +
	guides(size =  FALSE,
				 alpha = FALSE,
				 shape = guide_legend(nrow = 4))

full_model_tidy %>%
	mutate(model = "unk dad set at 0.5") %>%
	bind_rows(xdata %>%
							mutate(model = "individual bootstrap models")) %>%
	bind_rows(xdata_sum %>%
							mutate(model = "mean estimate 0f bootstrap")) %>%
	bind_rows(only_known_dads_tidy %>%
							mutate(model = "only known paternities"))  %>%
	filter(str_detect(term, "kid_with_mom")) %>%
	ggplot() +
	geom_point(aes(x = estimate, y = term,  color = sig, shape = model, size = model, group = model), position=position_dodge(.9)) +
	scale_shape_manual(values = c(20, 19, 17, 15)) +
	scale_size_manual(values = c(1,3, 3, 3)) +
	scale_alpha_manual(values = c(.2, 1)) +
	scale_color_manual(values = c("grey50", "black")) +
	cowplot::theme_cowplot() +
	theme(legend.position = "bottom") +
	guides(size =  FALSE,
				 alpha = FALSE,
				 shape = guide_legend(nrow = 4)) +
	facet_wrap(~term, scales = "free_y")

full_model_results %>%
	kable(format = "html", booktabs = T, longtable = T, digits = 3, escape = F) %>%
	kable_styling(font_size = 12)






