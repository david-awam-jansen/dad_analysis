xx <- read_csv("./data/who_grooms_bootstrap.csv")

xx %>%
	anti_join(select(who_grooms_step6, focal, focal_grp, start, end))
filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

xdata_females_with_social %>%
	filter(sname == 'EMC')


################################################################################
members <- tbl(babase, 'members') %>%
	filter(grp < 3) %>%
	select(membid, grp, sname, date, grpofresidency)

members <- members %>%
	left_join(group_history) %>%
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


xdata_females_with_social %>%
	filter(dad_overlap_l > 3) %>%
	select(sname, dad, birth, dad_overlap_l)

focal = "LAN"
focal_dad = 'IAG'
focal_birth = '2000-01-21'

## babase data



ggplot() +
	geom_segment(data =xdata_temp,
							 aes(x = 0, xend= dad_overlap_days,
							 		y=ordered_sname,
							 		yend = ordered_sname),
							 size = .5) +
	labs(x = "Duration of time spend in same grp as dad (years)",
			 y = "Focal individuals") +
	scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
	scale_y_continuous(breaks=c(12, 54, 93, 139, 213),
										 labels=c('12 (6%)', '54 (25%)', '93 (44%)',
										 				 '139 (65%)', '213 (100%)'),
										 expand = c(0, 0), limits = c(0, NA)) +
	# geom_segment(data = temp4,
	# 						 aes(x=class, xend = class,
	# 						 		y = 0, yend = class_max),
	# 						 linetype = "dashed", color = "red", size = 1) +
	# geom_segment(data = temp4,
	# 						 aes(x=0, xend = class,
	# 						 		y = class_max, yend = class_max),
	# 						 linetype = "dashed", color = "red", size =1) +
	# geom_segment(data = temp2,
	# 						 aes(x=0, xend = class,
	# 				 		     y = class_max, yend = class_max),
# 				         linetype = "dashed", color = "yellow", size =1) +
cowplot::theme_cowplot()




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













# social_index_values <- read_rds('./data/temp_social_index_values.rds')
# xdata_females_with_social <- read_csv('./data/xdata_females_with_social.csv')
#
# get_mean_observer_effort <-function(df) {
#
# 	if(nrow(df) > 1) {
# 	df %>%
# 		ungroup() %>%
# 		select(mean_AF_log2OE) %>%
# 		distinct() %>%
# 		pull()
# 	} else {
# 		NA
# 	}
# }

who_grooms <- social_index_values %>%
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



who_grooms6 %>%
	select(-zscore) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

library(lmerTest)

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

	v %>% as_tibble() %>%
		mutate(term = nam)
}

full_model_age <- glmer(data = who_grooms6
												, formula = does_groom ~ 1
												+ observer_effort
												## ages
												+ kid_age+ mom_age + poly(AMales_age, 2)
												## early adversity
												+ cumulative_adversity
												## male details
												+ is_dad
												#+ mean_ordrank
												+ estrous_c + daily_d_days_rate +
													offspring_years + nr_pdads
												+ previous_kid_with_mom
												+ (1|focal) +  (1|AMales)
												, weights = log(nr_days)
												, family = binomial
												, control=glmerControl(optimizer = "bobyqa",
																							 optCtrl=list(maxfun=2e5))
												, na.action = 'na.fail')

full_model_age %>%
	broom.mixed::tidy() %>%
	left_join(vif.mer(full_model_age), by = 'term') %>%
	kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
	kable_styling(font_size = 20) %>%
	kable_classic(full_width = F)

full_model_rank <- update(object = full_model_age
													, . ~ . - poly(AMales_age, 2) + mean_ordrank
													, data = who_grooms6)

full_model_rank %>%
	broom.mixed::tidy() %>%
	left_join(vif.mer(full_model_age), by = 'term') %>%
	kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
	kable_styling(font_size = 20) %>%
	kable_classic(full_width = F)

anova(full_model,  full_model_age, full_model_rank)







full_model_age_no_dad <- update(object = full_model_age,
																. ~ . -is_dad,
																data = who_grooms6 %>%  filter(is_dad == FALSE))

full_model %>%
	broom.mixed::tidy() %>%
	left_join(vif.mer(full_model), by = 'term') %>%
	kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
	kable_styling(font_size = 20) %>%
	kable_classic(full_width = F)






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

