library(RPostgreSQL)
library(tidyverse)

babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "jansen",
	dbname = "babase",
	password = "Bab00n3455")

biograph_l <- tbl(babase, "biograph") %>%  collect()

group_history <- tbl(babase, "groups_history") %>%
	select(grp = gid, permanent, impermanent, cease_to_exist, last_reg_census) %>%
	group_by(grp) %>%
	mutate(last_date = pmin(impermanent, cease_to_exist, last_reg_census, "2020-07-01")) %>%
	ungroup()

behave_gaps <- tbl(babase, "behave_gaps")

################################################################################

members <- tbl(babase, 'members') %>%
	filter(grp < 3) %>%
	select(membid, grp, sname, date)

members <- members %>%
	inner_join(group_history) %>%
	filter(date >= permanent & date <= last_date)

members_bg <- members %>%
	dplyr::left_join(tbl(babase, "behave_gaps"), by = "grp") %>%
	dplyr::filter(date >= gap_start & date <= gap_end) %>%
	select(membid, bgid)

members <- members %>%
	anti_join(members_bg, by = 'membid')

## To get it use it in analysis and/or save it you have to 'collect; the data
members_l <- members %>%
	collect()
################################################################################
selected_behave_gaps <- read_csv('../R01_grant/data/selected_behave_gaps.csv')
excluded_rank <- selected_behave_gaps %>%
	filter(exclude_rank) %>%
	select(bgid) %>%
	pull()

ranks <- tbl(babase, 'proportional_ranks') %>%
	filter(grp < 3)

ranks <- ranks %>%
	inner_join(group_history) %>%
	filter(rnkdate >= permanent & rnkdate <= last_date)

ranks_bg <- ranks %>%
	dplyr::left_join(tbl(babase, "behave_gaps"), by = "grp") %>%
	dplyr::filter(rnkdate >= gap_start & rnkdate <= gap_end) %>%
	select(rnkid, bgid) %>%
	filter(bgid %in% excluded_rank)

ranks <- ranks %>%
	anti_join(ranks_bg, by = 'rnkid')

## To get it use it in analysis and/or save it you have to 'collect; the data
ranks_l <- ranks %>%
	collect()

rankdates_l <- tbl(babase, 'rankdates') %>%  collect()
potential_dads_l <- tbl(babase, 'potential_dads') %>%  collect()
parents_l <- tbl(babase, 'parents') %>%  collect()
hybridgene_scores_l <- tbl(babase, "hybridgene_scores") %>%
	filter(hgaid == 2) %>%
	collect()
#####################################################################
####################################################################
groupsizes <- tbl(babase, "members") %>%
	filter(grp < 3) %>%
	group_by(grp, date) %>%
	summarise(groupsize = n()) %>%
	mutate(rnkdate = rnkdate(date)) %>%
	group_by(grp, rnkdate) %>%
	summarise(groupsize = mean(groupsize, na.rm = TRUE))

groupsize_adult_males <- tbl(babase, "members") %>%
	filter(grp < 3) %>%
	inner_join(tbl(babase, "rankdates")) %>%
	filter(date >= ranked) %>%
	group_by(grp, date) %>%
	summarise(nr_adult_males = n()) %>%
	mutate(rnkdate = rnkdate(date)) %>%
	group_by(grp, rnkdate) %>%
	summarise(nr_adult_males = mean(nr_adult_males, na.rm = TRUE))

weather <- tbl(babase, "wreadings") %>%
	inner_join(tbl(babase, "raingauges")) %>%
	inner_join(tbl(babase, "tempmins")) %>%
	inner_join(tbl(babase, "tempmaxs")) %>%
	mutate(date = !!dbplyr::build_sql(con = babase, "CAST(wrdaytime AS DATE)")) %>%
	select(date, rain, tempmin, tempmax) %>%
	group_by(date) %>%
	summarise(rain = mean(rain, na.rm=TRUE),
						tempmin = min(tempmin, na.rm=TRUE),
						tempmax = max(tempmax, na.rm=TRUE)) %>%
	collect() %>%
	arrange(date)

get_mean_temp <- function(date) {
	sample_date  = lubridate::ymd(date)
	start_date = sample_date - months(1)

	weather %>%
		filter(date >= start_date &
					 	date <= sample_date) %>%
		select(tempmax) %>%
		summarise(tempmax = mean(tempmax, na.rm = TRUE)) %>%
		pull()
}

# hormones
# storage time as fecal powder
# storage time in methanol
tbl(babase, dbplyr::in_schema('fecal', 'prep')) %>%
	select(sid, sname, date, meextract) %>%
	mutate(years_collected_to_meth = as.numeric(meextract - date)/365.25) %>%
	inner_join(select(tbl(babase, "biograph"), sname, birth, sex)) %>%
	filter(sex == "M") %>%
	inner_join(tbl(babase, dbplyr::in_schema('fecal', 'results')) %>%
						 	select(sid, e2, gc, t, th, gc_date)) %>%
	mutate(years_meth_to_assay = as.numeric(gc_date - meextract)/365.25) %>%
	# age & age squared
	mutate(age_at_sampling = as.numeric((date - birth)/365.25),
				 age_at_sampling2 = age_at_sampling^2) %>%
	filter(e2 > 0 | gc > 0 | t > 0 | th > 0) %>%
	mutate(rnkdate = rnkdate(date)) %>%
	inner_join(select(tbl(babase, 'members'), grp, date, sname)) %>%
	# group size
	# group size squared
	inner_join(groupsizes) %>%
	mutate(groupsize2 = groupsize ^2) -> temp1

temp1 %>%
	inner_join(groupsize_adult_males) %>%
	mutate(nr_adult_males2 = nr_adult_males ^ 2) %>%
	# proportional rank
	inner_join(tbl(babase, 'proportional_ranks') %>%  filter(rnktype == 'ADM')) %>%
	# categorical variable coding whether the individual was alpha or not
	mutate(is_alpha = ordrank == 1) %>%
	select(everything())-> temp1

temp2 <- temp1 %>%  collect()

temp_weather <- temp2 %>%
	select(date) %>%
	distinct() %>%
	mutate(season = if_else(lubridate::month(date) >= 5 & lubridate::month(date) <= 10, "DRY", "WET"),
				 max_temp = map_dbl(.x = date, .f = get_mean_temp))

temp3 <- temp2 %>%
	inner_join(temp_weather, by = 'date')

temp4 <- temp3 %>%
	select(sid, sname, date,
				 meextract, gc_date,
				 years_collected_to_meth, years_meth_to_assay,
				 age_at_sampling, age_at_sampling2, e2, gc, t, th,
				 rnkdate, grp, groupsize, groupsize2, nr_adult_males, ordrank, is_alpha,
				 season, max_temp) 	%>%
		filter(complete.cases(.))


temp5 <- temp4 %>%
	pivot_longer(cols = c(e2, gc, t, th),
							 names_to = "hormones",
							 values_to = "value") %>%
	filter(value > 0)

temp6 <- temp5 %>%
	group_by(hormones) %>%
	nest() %>%
	mutate(data = map(.x = data, .f = rownames_to_column))

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
get_hormome_models <- function(df) {
	df<- df %>%
		mutate_at(c("groupsize", "groupsize2", "ordrank", "max_temp"), scale2) %>%
		mutate(value = log(value))

	lmerTest::lmer(data = df,
			 value ~ years_collected_to_meth +years_meth_to_assay +
			 	age_at_sampling + age_at_sampling2 +
			 	groupsize + groupsize2 + ordrank +  season + max_temp + is_alpha +
			 	(1|sname) + (1|grp))
}

get_hormome_models2 <- function(df) {
	df<- df %>%
		mutate_at(c("groupsize", "groupsize2", "ordrank", "max_temp"), scale2)

	lmerTest::lmer(data = df,
								 value ~ years_collected_to_meth +years_meth_to_assay +
								 	age_at_sampling + age_at_sampling2 +
								 	groupsize + groupsize2 + ordrank +  season + max_temp + is_alpha +
								 	(1|sname) + (1|grp))
}

add_resid <- function(df, residuals) {
	df %>%
		inner_join(residuals %>%
							 rownames_to_column() %>%
							 select(rowname, contains(".")),
							 by = 'rowname')
}

hormone_models <- temp6 %>%
	mutate(models = map(.x = data, .f = get_hormome_models)) %>%
	mutate(models2 = map(.x = data, .f = get_hormome_models2))

hormone_residuals <- hormone_models %>%
	mutate(residuals = map(models, broom.mixed::augment)) %>%
	mutate(data = map2(.x = data, .y = residuals, add_resid)) %>%
	unnest(data) %>%
	select(hormones, sname, date, grp, value, .resid)

hormone_results <- hormone_models %>%
	mutate(results = map(models, broom.mixed::tidy)) %>%
	select(hormones, results) %>%
	unnest(cols = c(results))
##########################################################################

## extra data based on blog and Silk talk

d_days <- tbl(babase, 'members') %>%
	inner_join(tbl(babase, 'biograph') %>%
						 	filter(sex == 'F') %>%
						 	select(sname), by = "sname") %>%
	inner_join(tbl(babase, 'repstats') %>%
						 	select(sname, date, repstate = state)
						 , by = c("sname", "date")) %>%
		left_join(tbl(babase, 'cycstats'), by = c("sname", "date")) %>%
		left_join(tbl(babase, 'cycpoints_cycles'),
							by = c("sname", "date")) %>%
	filter(code == "D") %>%
	group_by(grp, date) %>%
	summarise(nr_d_days = n()) %>%
	collect()

sibling_names <- tbl(babase, "parents") %>%
	inner_join(select(tbl(babase, "biograph"),
										kid = sname, birth, statdate)) %>%
	select(dad, mom, kid, birth, statdate) %>%
	collect() %>%
	arrange(mom, birth) %>%
	group_by(mom) %>%
	mutate(previous_sibling = lag(kid),
				 next_sibling = lead(kid))



save(biograph_l, members_l, ranks_l, rankdates_l,
		 potential_dads_l, parents_l,
		 hormone_residuals, hormone_results,
		 hybridgene_scores_l,
		 d_days,
		 sibling_names,
		 file = "./data/babase_data_for_DSI_paternal_analysis_with_results.RData")


load(file = "./data/babase_data_for_DSI_paternal_analysis_with_results.RData")
anubis_estimates <- read_csv("./data/genome_wide_anubis_estimates_for_David_23Mar2021.csv")

get_maternal_prop <- function(mom, focal_grp, date) {
	kid_rnkdate <- floor_date(ymd(date), unit = 'month')
	ranks_l %>%
		filter(sname == mom
					 & rnktype == 'ADF'
					 &  grp == focal_grp
					 & rnkdate == kid_rnkdate
		) %>%
		select(proprank) %>%
		pull()
}

members_l_AM <- members_l %>%
	inner_join(rankdates_l) %>%
	filter(date > ranked) %>%
	select(sname, grp, date, ranked)

get_focal_days  <- function(focal, focal_grp, start, end) {
	start_date <- ymd(start)
	end_date <- ymd(end)

	members_l %>%
		filter(sname == focal
					 & date >= start_date
					 & date <= end_date
					 & grp == focal_grp) %>%
		group_by(sname) %>%
		summarise(nr_days = n(), .groups = 'drop') %>%
		pull()
}

get_all_males <- function(focal_grp, start, end) {
	start_date <- ymd(start)
	end_date <- ymd(end)

	members_l_AM %>%
		filter(date >= start_date
					 & date <= end_date
					 & grp == focal_grp) %>%
		group_by(sname) %>%
		summarise(nr_days = n(), .groups = 'drop')
}


get_mean_rank <- function(AMales, focal_grp, start, end) {
	start_date <- ymd(start)
	end_date <- ymd(start)

	ranks_l %>%
		filter(sname == AMales) %>%
		filter(rnktype == 'ADM') %>%
		filter(rnkdate >= floor_date(ymd(start), "month") &
					 	rnkdate <= floor_date(ymd(end), "month")) %>%
		select(ordrank, proprank) %>%
		summarise(mean_ordrank = mean(ordrank),
							mean_proprank = round(mean(proprank), 2),
							.groups = 'drop')
}

get_d_days_rate <- function(AMales, focal_grp, start, end) {
	start <- ymd(start)
	end <- ymd(end)

	members_l_AM %>%
		filter(sname == AMales
					 & 	grp == focal_grp
					 & date >=start
					 & date <= end) %>%
		left_join(d_days, by = c("grp", "date")) %>%
		summarise(daily_d_days_rate = sum(nr_d_days, na.rm = TRUE)/n()) %>%
		pull()
}

get_kids <- function(focal, focal_birth, focal_mom, AMales) {
	df <- parents_l %>%
		arrange(zdate) %>%
		select(kid, mom, dad) %>%
		inner_join(select(biograph_l, kid = sname, kid_birth = birth), by = "kid") %>%
		filter(mom == focal_mom) %>%
		filter(kid != focal) %>%
		mutate(age_diff = as.numeric(kid_birth - ymd(focal_birth)))

	previous_kid_with_mom_df <- df %>%
		filter(age_diff < 0)

	if(nrow(previous_kid_with_mom_df) >= 1) {
		previous_kid_with_mom_df <- previous_kid_with_mom_df %>%
			filter(age_diff == max(age_diff))

		previous_pdads <- potential_dads_l %>%
			filter(kid %in% previous_kid_with_mom_df$kid) %>%
			select(pdad) %>%
			pull()

		previous_kid_with_mom <- tibble(previous_kid_with_mom =
					 	case_when(nrow(previous_kid_with_mom_df) == 0 ~ 0,
					 						!is.na(previous_kid_with_mom_df$dad) & previous_kid_with_mom_df$dad == AMales ~ 1,
					 						!is.na(previous_kid_with_mom_df$dad) & previous_kid_with_mom_df$dad != AMales ~ 0,
					 						is.na(previous_kid_with_mom_df$dad) & AMales %in% c(previous_pdads) ~ 0.5,
					 						is.na(previous_kid_with_mom_df$dad) & !(AMales %in% c(previous_pdads)) ~ 0,
					 						TRUE ~ 99))
	} else {
			previous_kid_with_mom =tibble(previous_kid_with_mom = 0)
			}

	next_kid_with_mom_df <- df %>%
		filter(age_diff > 0)

	if(nrow(next_kid_with_mom_df) >= 1) {

		next_kid_with_mom_df <- next_kid_with_mom_df %>%
			filter(age_diff == min(age_diff))

		next_pdads <- potential_dads_l %>%
			filter(kid %in% next_kid_with_mom_df$kid) %>%
			select(pdad) %>%
			pull()

		next_kid_with_mom <- tibble(next_kid_with_mom =
					 	case_when(!is.na(next_kid_with_mom_df$dad) & next_kid_with_mom_df$dad == AMales ~ 1,
					 						!is.na(next_kid_with_mom_df$dad) & next_kid_with_mom_df$dad != AMales ~ 0,
					 						is.na(next_kid_with_mom_df$dad) & AMales %in% c(next_pdads) ~ 0.5,
					 						is.na(next_kid_with_mom_df$dad) & !(AMales %in% c(next_pdads)) ~ 0,
					 						nrow(next_kid_with_mom_df) == 0 ~ 0,
					 						TRUE ~ 99))
	} else {
				next_kid_with_mom =tibble(next_kid_with_mom = 0)
				}


	 	bind_cols(previous_kid_with_mom,next_kid_with_mom)
}


get_nr_pdads <- function(focal) {
	potential_dads_l %>%
		filter(kid == focal) %>%
		group_by(kid) %>%
		summarise(nr_potential_dads = n()
						 ,.groups = 'drop'
						 , by = "kid") %>%
	select(nr_potential_dads) %>%
	pull()
}

members_juveniles <- members_l %>%
	inner_join(select(biograph_l, sname, birth), by = 'sname') %>%
	mutate(age = as.numeric(date- birth)/365.25) %>%
	filter(age < 4) %>%
	inner_join(select(parents_l, sname = kid, dad))

focal_sname = 'LIZ'
focal_grp = 2.2
focal_dad = 'ROC'
AMales = 'ROC'
start = "1996-10-19"
end = "1997-10-18"



get_offspring_years <- function(focal_sname, AMales, focal_grp, start, end) {
	start <- lubridate::ymd(start)
	end <- lubridate::ymd(end)

	members_juveniles  %>%
		filter(dad == AMales &
					 	grp == focal_grp) %>%
		filter(date >= start & date <= end) %>%
		select(sname, date, birth) %>%
		filter(sname != focal_sname) %>%
		arrange(date) %>%
		nrow()/365.25
}



who_grooms <- social_index_values %>%
	select(focal = sname, focal_grp = grp, start, end) %>%
	inner_join(select(xdata_females_with_social,
										focal = sname, mom, dad, maternal_loss, maternal_rank, maternal_SCI_F,
										density, sibling, drought, cumulative_adversity), by = "focal") %>%
	inner_join(select(biograph_l, focal = sname, focal_birth = birth), by = "focal") %>%
	inner_join(select(biograph_l, mom = sname, mom_birth = birth), by = "mom") %>%
	mutate(kid_age = as.numeric(start - focal_birth)/365.25) %>%
	mutate(mom_age = as.numeric(start - mom_birth)/365.25) %>%  # mom age
	mutate(nr_focal_days = pmap_dbl(
		.l = list(focal, focal_grp, start, end),
		.f = get_focal_days))

who_grooms2 <- who_grooms %>%
	mutate(all_males = pmap(.l = list(focal_grp, start, end), .f = get_all_males)) %>%
	unnest(cols = c(all_males)) %>%
	left_join(social_index_values %>%
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
	mutate(AMales_age = as.numeric((start - AMales_birth)/365.25))

who_grooms3 <- who_grooms2 %>%
	mutate(AMales_rank = pmap(.l = list(AMales, focal_grp, start, end),
														.f = get_mean_rank)) %>%
	unnest(cols = c(AMales_rank), keep_empty = TRUE)

who_grooms4 <- who_grooms3 %>%
	left_join(select(potential_dads_l, focal = kid, AMales = pdad,
									 pdad_status = status,
									 estrous_presence, estrous_me, estrous_c),
						by = c('AMales', "focal")) %>%
	mutate(estrous_presence =if_else(is.na(estrous_presence),
																	 0, estrous_presence)) %>%
	mutate(daily_d_days_rate = pmap_dbl(.l = list(AMales, focal_grp, start, end),
																			.f = get_d_days_rate)) %>%
	mutate(nr_pdads = map_dbl(.x = focal, .f = get_nr_pdads)) %>%
	mutate(offspring_years = pmap_dbl(
	.l = list(focal, AMales, focal_grp, start, end),
	.f = get_offspring_years))

who_grooms5 <- who_grooms4 %>%
	ungroup() %>%
	mutate(kids = pmap(.l = list(focal, focal_birth, mom, AMales),
										 .f = get_kids)) %>%
  unnest(cols = c(kids), keep_empty = TRUE)

who_grooms5 <- who_grooms5 %>%
	left_join(select(hybridgene_scores_l, AMales = sname, hybrid_score = score)) %>%
	left_join(select(anubis_estimates, AMales = sname, anubis_admix, notes) %>%
							filter(is.na(notes) & !is.na(anubis_admix)))

who_grooms6 <- who_grooms5 %>%
	select(focal, focal_birth, focal_grp, mom, dad, AMales, start, end, nr_days
				 ## ages
				 , kid_age, mom_age, AMales_age,
				 ## early adversity
				 , maternal_loss, maternal_rank, maternal_SCI_F, density, sibling, drought, cumulative_adversity

				 ## male details
				 , is_dad, mean_ordrank, mean_proprank,
				 , estrous_presence, , daily_d_days_rate, offspring_years, nr_pdads
				 , previous_kid_with_mom, next_kid_with_mom ## see bootstrap

				 ## for subset
				 , hybrid_score, anubis_admix

				 ## responses
				 , does_groom,  zscore) %>%
	filter(!is.na(mean_proprank)) %>%
	filter(nr_days >= 30)








full_model <- glmer(data = who_grooms_bootstrap
						, formula = does_groom ~ 1
						+ kid_age + mom_age + AMales_age
						+ is_dad
						+ cumulative_adversity
						+ mean_ordrank + estrous_presence + daily_d_days_rate
						+ previous_kid_with_mom + next_kid_with_mom+
							+ (1|focal) +  (1|AMales)
						, weight = log(nr_days)
						, family = binomial
						, control=glmerControl(optimizer = "bobyqa",
																	 optCtrl=list(maxfun=2e5)),
						, na.action = 'na.fail')

library(nlme)



full_modelstan <- rstanarm::stan_glmer(data = who_grooms6
										, formula = does_groom ~ 1
										+ kid_age + mom_age + AMales_age
										+ is_dad
										+ cumulative_adversity
										+ mean_ordrank + estrous_presence + daily_d_days_rate
										+ previous_kid_with_mom + next_kid_with_mom+
											+ (1|focal) +  (1|AMales)
										, offset = log(nr_days)
										, family = binomial
										# , control=glmerControl(optimizer = "bobyqa",
										# 											 optCtrl=list(maxfun=2e5)),
										, na.action = 'na.fail')






full_model_no_weights <- glmer(data = who_grooms6
										, formula = does_groom ~ 1
										+ kid_age + mom_age + AMales_age
										+ is_dad
										+ cumulative_adversity
										+ mean_ordrank + estrous_presence + daily_d_days_rate
										+ previous_kid_with_mom + next_kid_with_mom+
											+ (1|focal) +  (1|AMales)
										, family = binomial
										, control=glmerControl(optimizer = "bobyqa",
																					 optCtrl=list(maxfun=2e5)),
										, na.action = 'na.fail')



who_grooms6_hybrid <- who_grooms6 %>%
	filter(!is.na(hybrid_score))

full_model_hybrid <- update(object = full_model,
														. ~ . + hybrid_score,
														data = who_grooms6_hybrid)


full_model_dredged <- MuMIn::dredge(global.model = full_model
																		, rank = 'AICc')

full_model_best <- MuMIn::get.models(full_model_dredged, subset = delta < 2)


who_grooms6_no_dad <- who_grooms6 %>%
	filter(is_dad == FALSE)

full_model_no_dad <- update(object = full_model,
			 . ~ . - is_dad,
			 data = who_grooms6_no_dad)

who_grooms6_strength <- who_grooms6 %>%
	filter(does_groom == TRUE)

full_model_only_dad <- update(object = full_model,
															. ~ . - is_dad - estrous_presence,
															data = who_grooms6_only_dad)

full_model_groom_strength <- lmer(data = who_grooms6_strength
										, formula = zscore ~ 1
										+ kid_age + mom_age + AMales_age
										+ is_dad
										+ cumulative_adversity
										+ mean_ordrank + estrous_presence + daily_d_days_rate
										+ previous_kid_with_mom + next_kid_with_mom+
											+ (1|focal) +  (1|AMales)
										, weight = log(nr_days)
										, control=lmerControl(optimizer = "bobyqa",
																					 optCtrl=list(maxfun=2e5)),
										, na.action = 'na.fail')

vif.mer <- function (fit) {
	## adapted from rms::vif

	v <- vcov(fit)
	nam <- names(fixef(fit))

	## exclude intercepts
	ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
	if (ns > 0) {
		v <- v[-(1:ns), -(1:ns), drop = FALSE]
		nam <- nam[-(1:ns)]
	}

	d <- diag(v)^0.5
	v <- diag(solve(v/(d %o% d)))
	names(v) <- nam
	v
}


vif.mer(full_model)

write_csv(who_grooms_bootstrap, "./data/who_grooms_bootstrap.csv")

save(who_grooms6,
		 who_grooms6_hybrid,
		 who_grooms6_strength,
		 who_grooms6_no_dad,
		 full_model,
		 full_model_no_weights,
		 full_model_hybrid,
		 full_model_no_dad,
		 full_model_only_dad,
		 full_model_groom_strength,
		 full_model_dredged,
		 full_model_best,
		 file = "./data/who_grooms.RData")

full_model %>%
	tidy() %>%
	filter(effect == 'fixed') %>%
	select(-effect, -group) %>%
	kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
	kable_styling(font_size = 20) %>%
	kable_classic(full_width = F) %>%
	footnote(general = "The response variable is if male groomed wih juvenile (0/1)",
					 general_title = "Note:",
					 footnote_as_chunk = T,
					 title_format = c("italic", "underline"))

xx <-
best_models <- subset(xx, delta <3)

# select models using Royall's 1/8 rule for strength of evidence
# IMPORTANT: Weights have been renormalized!!
best_models <- subset(xx, 1/8 < weight/max(xx$weight))

importance(best_models)

# Model average using all candidate models, always use revised.var = TRUE
MA.ests<-model.avg(best_models, revised.var = TRUE)
MA.ests

best_models %>%
	str()


	tidy() %>%
	kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
	kable_styling(font_size = 20) %>%
	kable_classic(full_width = F) %>%
	footnote(general = "The response variable is if male groomed wih juvenile (0/1)",
					 general_title = "Note:",
					 footnote_as_chunk = T,
					 title_format = c("italic", "underline"))


varying.link <- list(family = alist(logit = binomial("logit"),
																			probit = binomial("probit"), cloglog = binomial("cloglog") ))
library("snow")
# Set up the cluster
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 6), type = clusterType))
clusterExport(clust, "who_grooms6")
library(MuMIn)
pdd <- pdredge(full_model, cluster = clust, varying = varying.link)


