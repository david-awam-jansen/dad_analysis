## all functions

get_DSI_values <- function(df) {
	step1 <- df %>%
		filter(focal_check == TRUE) %>%
		arrange(-zscore, paternal_groom) %>%
		ungroup()

	DSI_paternal <- step1 %>%
		filter(paternal_groom == "paternal") %>%
		pull(zscore)

	DSI_M <- step1 %>%
		filter(str_detect(paternal_groom, "paternal")) %>%
		top_n(3, wt = zscore) %>%
		summarise(sum(zscore)/n()) %>%
		pull()

	DSI_Mtop <- step1 %>%
		filter(str_detect(paternal_groom, "paternal")) %>%
		slice(1) %>%
		pull(zscore)

	DSI_Mde <- step1 %>%
		filter(str_detect(paternal_groom, "non_paternal")) %>%
		top_n(3, wt = zscore) %>%
		summarise(sum(zscore)/n()) %>%
		pull()

	DSI_Mde_top <- step1 %>%
		filter(str_detect(paternal_groom, "non_paternal")) %>%
		slice(1) %>%
		pull(zscore)

	tibble(DSI_paternal,
				 DSI_M,
				 DSI_Mtop,
				 DSI_Mde,
				 DSI_Mde_top)
}

has_real_partners <- function(df) {
	df %>%
		summarise(nr_partners = sum(grooming_partner != 'XXX') > 0 ) %>%
		pull()
}

get_top_grooms <- function(df) {
	df %>%
		filter(focal_check) %>%
		filter(str_detect(paternal_groom, "paternal")) %>%
		group_by(paternal_groom) %>%
		arrange(zscore) %>%
		slice(1) %>%
		select(paternal_groom, grooming_partner = sname, zscore)  %>%
		ungroup()

}

get_dad_top <- function(df) {

	if(nrow(df) == 0) {

		FALSE
	} else {
		df %>%
			arrange(zscore, paternal_groom) %>%
			slice(1) %>%
			mutate(dad_top = paternal_groom == 'paternal') %>%
			select(dad_top) %>%
			pull()
	}
}

check_nr_dyads <- function(df) nrow(df)

get_zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

get_overlap <- function(focal, focal_partner, start, end, per_group = FALSE) {
	start = ymd(start)
	end = ymd(end)

	if(per_group != FALSE) 	members_data <- members_partner_overlap_l %>% filter(grp == per_group | grp == grpofresidency)
	if(per_group == FALSE) 	members_data <- members_partner_overlap_l


	df <- members_partner_overlap_l %>%
		filter(sname == focal) %>%
		arrange(date) %>%
		filter(between(date, start, end)) %>%
		select(date, focal = sname, focal_grpofresidency = grpofresidency) %>%
		full_join(members_partner_overlap_l %>%
								filter(sname == focal_partner) %>%
								filter(between(date, start, end)) %>%
								select(date, partner = sname, partner_grp = grp, partner_grpofresidency = grpofresidency),
							by = "date")  %>%
		mutate(same_grp = focal_grpofresidency == partner_grpofresidency |
					 	focal_grpofresidency == partner_grp)

	df %>%
		summarise(partner_overlap_years = sum(same_grp, na.rm = TRUE)/365.25,
							focal_days_known_grp = sum(!is.na(focal_grpofresidency))/365.25,
							partner_overlap_proportional = partner_overlap_years/focal_days_known_grp)
}

check_overlap_plot <- function(focal, focal_partner, start, end) {
	start = ymd(start)
	end = ymd(end)

	df <- members_partner_overlap_l %>%
		filter(sname == focal) %>%
		arrange(date) %>%
		filter(between(date, start, end)) %>%
		select(date, focal = sname, focal_grpofresidency = grpofresidency) %>%
		full_join(members_partner_overlap_l %>%
								filter(sname == focal_partner) %>%
								filter(between(date, start, end)) %>%
								select(date, partner = sname, partner_grp = grp, partner_grpofresidency = grpofresidency),
							by = "date")  %>%
		mutate(same_grp = focal_grpofresidency == partner_grpofresidency |
					 	focal_grpofresidency == partner_grp)

	df %>%
		select(date, same_grp, focal_grpofresidency, dad_grpofresidency) %>%
		pivot_longer(names_to = "sname", values_to = "grp", focal_grpofresidency:dad_grpofresidency) %>%
		inner_join(select(groups_l, grp = gid, social_group = name), by = "grp") %>%
		mutate(sname = if_else(str_detect(sname, 'focal'), 'Juvenile', 'Dad')) %>%
		mutate(age = as.numeric(date - focal_birth)/365.25) %>%
		ggplot(aes(x = age, y = sname, color = social_group)) +
		geom_point() +
		geom_vline(xintercept = df1$first_time_not_same_grp, color = 'red', size = 2) +
		geom_vline(xintercept = df1$last_time_not_same_grp, color = 'green', size = 2) +
		cowplot::theme_cowplot() +
		labs(x = "Age of juvenile in years",
				 y ="") +
		theme(legend.position = "bottom")

}

## for who groomed
get_focal_days  <- function(focal, focal_grp, start, end) {
	start_date <- ymd(start)
	end_date <- ymd(end)

	members_partner_overlap_l %>%
		filter(sname == focal
					 & date >= start_date
					 & date <= end_date
					 & (grp == focal_grp |  grpofresidency== focal_grp )) %>%
		group_by(sname) %>%
		summarise(nr_days = n(), .groups = 'drop') %>%
		pull()
}

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

get_all_males <- function(focal_grp, start, end) {
	start_date <- ymd(start)
	end_date <- ymd(end)

	members_AM_l %>%
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

	members_AM_l %>%
		filter(sname == AMales
					 & 	grp == focal_grp
					 & date >=start
					 & date <= end) %>%
		left_join(d_days_l, by = c("grp", "date")) %>%
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

## for bootstrap
get_maternal_sibling <- function(focal_sname, focal_birth, focal_grp, focal_mom, period) {

	df <- parents_l %>%
		inner_join(select(biograph_l, kid = sname, kid_birth = birth), by = "kid") %>%
		filter(mom == focal_mom) %>%
		arrange(kid_birth)

	if(period == "previous") {
		df2 <- df %>%
			filter(kid_birth < focal_birth) %>%
			#filter(kid != focal_sname) %>%
			filter(kid_birth == max(kid_birth))
	} else {
		df2 <- df %>%
			filter(kid_birth >= focal_birth) %>%
			filter(kid != focal_sname) %>%
			filter(kid_birth == min(kid_birth))
	}

	df2 %>%
		inner_join(potential_dads_l %>%
							 	group_by(kid) %>%
							 	summarise(nr_potential_dads = n())
							 ,.groups = 'drop'
							 , by = "kid") %>%
		select(offspring_sname = kid, nr_potential_dads)
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

get_offspring_years <- function(focal_sname, AMales, focal_grp, start, end) {
	start <- lubridate::ymd(start)
	end <- lubridate::ymd(end)

	members_juveniles_l  %>%
		filter(dad == AMales &
					 	grp == focal_grp) %>%
		filter(date >= start & date <= end) %>%
		select(sname, date, birth) %>%
		filter(sname != focal_sname) %>%
		arrange(date) %>%
		nrow()/365.25
}

get_male_present_count <-  function(focal_kid, focal_grp, start, end) {

	start_date <- ymd(start)
	end_date <- ymd(end)

	members_AM_l %>%
		filter(date >= start_date
					 & date <= end_date
					 & grp == focal_grp) %>%
		distinct(sname) %>%
		left_join(select(potential_dads_l, kid, sname = pdad) %>%
								filter(kid == focal_kid),
							by = 'sname') %>%
		summarise(nr_males = n(),
							pdad_present = sum(!is.na(kid)==TRUE),
							.groups = 'drop')
}

get_random_draws <- function(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = 1000, period)  {

	paternity <- c(1, rep(0, times = nr_potential_dads -1))

	xx <-replicate(draws,sample(paternity, pdad_present ,replace = FALSE))

	if(pdad_present > 1) {
		xx <- xx %>%
			as_tibble() %>%
			mutate(AMales = c(pdads_sname))
	} else {
		xx <- t(xx) %>%
			as_tibble() %>%
			mutate(AMales = c(pdads_sname))
	}

	names(xx)[1:draws] <- paste0(period, "_", seq(1:draws))
	return(xx)
}

get_pdads <- function(df, period) {
	if(period == "previous") {
		df %>%
			filter(previous_kid_with_mom == 0.5) %>%
			select(AMales) %>%  pull()
	} else {
		df %>%
			filter(next_kid_with_mom == 0.5) %>%
			select(AMales) %>%  pull()
	}
}

assign_random_draws <- function(df, previous_set, next_set) {
	if(!is.null(previous_set)) {
		df <- df %>%
			left_join(previous_set, by = 'AMales')
	}

	if(!is.null(next_set)) {
		df <- df %>%
			left_join(next_set, by = 'AMales')
	}
	return(df)
}

## for paper
fix_terms <-function(term) {
	term %>%
		str_replace("kid_age", "Age of juvenile") %>%
		str_replace("dad_age", "Age of dad") %>%
		str_replace("mom_age", "Maternal age") %>%
		str_replace("maternal SCI_M", "Maternal SCI_M") %>%
		str_replace("cumulative_adversity", "Cum. adversity of juvenile") %>%
		str_replace("maternal_proprank", "Maternal rank (at birth)") %>%
		str_replace("paternal_sibling_years", "Number of paternal siblings") %>%
		str_replace("offspring_years", "Number of kids male has in group") %>%
		str_replace("full_sibling_years", "Number of full siblings") %>%
		str_replace("dad_rank", "Ordinal rank of dad") %>%
		str_replace("hybrid_score", "Hybrid score of of the male") %>%
		str_replace("anubis_admix", "Anubis admix score of the male") %>%
		str_replace("nr_potential_dads", "Number of potential dads") %>%
		str_replace("maternal_SCI_F_binaryTRUE", "Socially isolated mom") %>%
		str_replace("maternal_lossTRUE", "Maternal_loss") %>%
		str_replace("maternal_rank_binaryTRUE", "Low ranking mom") %>%
		str_replace("density_binaryTRUE", "Born in large group") %>%
		str_replace("siblingTRUE", "Has competing maternal sibling") %>%
		str_replace("droughtTRUE", "Born in drought") %>%
		str_replace("nr_d_days", "Daily rate of d days") %>%
		str_replace("previous_sibling", "Overlap with previous maternal sibling") %>%
		str_replace("next_sibling", "Overlap with next maternal sibling") %>%
		str_replace("focal_age", "Age class of the juvenile") %>%
		str_replace("is_dadTRUE", "Male is the dad") %>%
		str_replace("AMales_age", "Age of male") %>%
		str_replace("mean_ordrank", "Rank of the male") %>%
		str_replace("mean_proprank", "Proportional rank of the male") %>%
		str_replace("previous_kid_with_mom", "Mother had a kid with male") %>%
		str_replace("next_kid_with_mom", "Mother will have a kid with male") %>%
		str_replace("estrous_presence", "Potential dad") %>%
		str_replace("daily_d_days_rate", "Daily rate of d days") %>%
		str_replace("nr_pdads", "Number of potential dads at conception") %>%
		str_replace("jDSI_paternal", "jDSI with dad") %>%
		str_replace("jDSI_Mtop", "top male jDSI") %>%
		str_replace("jDSI_Mde_top", "top non-dad male jDSI") %>%
		str_replace("dad_overlap_l", "years of co-residency with dad")  %>%
		str_replace("observer_effort", "Observer effort") %>%
		str_replace("poly\\(Age of male, 2\\)1", "First polynomial of male age") %>%
		str_replace("poly\\(Age of male, 2\\)2", "Second polynomial of male age") %>%
		str_replace('estrous_c', "Male had consort with mother during fertile period") %>%
		# str_replace('estrous_me', "The number of mounts and ejaculation interactions observerd between the mother and the male during the fertitle period ")
	str_replace('estrous_me', "The number of mounts")
}

get_coxph_model <- function(response) {
	coxph(data=xdata_females_with_social, as.formula(paste0("Surv(statage, adult_survival_status) ~ ", response)))
}

get_coxzph <- function(x) {
	tail(cox.zph(x)$table, 1)[,3]
}

get_confint <- function(model) {
	model %>%
		confint() %>%
		as_tibble() %>%
		setNames(c("conf.low", "conf.high"))
}

make_table <- function(df, fontSize = 15, style = 'html', figure_zoom = 1) {
	t1 <- df %>%
		mutate(HR = round(exp(estimate), 3),
					 conf.low = round(exp(conf.low), 3),
					 conf.high = round(exp(conf.high), 3),
					 sig_stars = case_when(p.value < 0.001  ~ "***",
					 											p.value < 0.01  ~ "***",
					 											p.value < 0.05  ~ "*",
					 											p.value < 0.01  ~ ".",
					 											TRUE ~ ""),
					 p.value = round(p.value, 3),
					 sig_stars = paste("<strong>", "<sup>", sig_stars,"</sup>", "</strong>", "<br>", sep =""),
					 HR = paste(HR, sig_stars, "(", conf.low, "-", conf.high,")", sep = "")) %>%
		select(formula, term, HR, AICc, model_check) %>%
		arrange(AICc) %>%
		mutate(min_AICc = min(AICc)) %>%
		mutate(dAICc = AICc - min_AICc) %>%
		select(formula, term, HR, AICc, dAICc, model_check) %>%
		#mutate(term = str_remove(term, "TRUE")) %>%
		mutate(term = map_chr(.x = term, .f = fix_terms)) %>%
		mutate(formula =
					 	str_replace(formula, "maternal_loss", "maternal loss"),
					 formula =
					 	str_replace(formula, "jDSI_M","xxx"),
					 formula =
					 	str_replace(formula, "xxxtop",
					 							"strength of strongest juvenile dyadic strength with top male
                    (including dad; top male jDSI)"),
					 formula =
					 	str_replace(formula, "xxxde_top",
					 							"strength of strongest juvenile dyadic strength with top male
                    (excluding dad; top non-dad jDSI)"),
					 formula =
					 	str_replace(formula, "xxxde","juvenile dyadic strength with top3 males (excluding dad; jDSI_M)"),
					 formula =
					 	str_replace(formula, "xxx","juvenile dyadic strength with top3 males (jDSI_M)"),

					 formula =
					 	str_replace(formula, "cumulative_adversity","cumulative adversity"),
					 formula =
					 	str_replace(formula, "jDSI_paternal","juvenile dyadic strength with dad (jDSI with dad)"),
					 formula =
					 	str_replace(formula, "jSCI_M","xxx"),
					 formula =
					 	str_replace(formula, "xxxde","juvenile social connectedness to males excluding dad (jSCI_Mde)"),
					 formula =
					 	str_replace(formula, "xxx","juvenile social connectedness to males (jSCI_M)"),
					 formula =
					 	str_replace(formula, "dad_overlap_l","proportion of juvenile period shared with dad in same grp"),
					 formula =
					 	str_replace(formula, "dad_overlap","proportions of juvenile period shared with dad in same grp"),
					 formula = paste0('Adult female survival - ', formula)) %>%
		group_by(term) %>%
		mutate(HR = ifelse(str_detect(HR, "\\*"),
											 paste('<span style="background-color:Yellow;">',
											 			HR, '</span>', sep = ""),
											 HR)) %>%
		pivot_wider(names_from = term, values_from = HR) %>%
		select(!(AICc:model_check), AICc, dAICc, model_check) %>%
		mutate(formula = paste0("* ", formula)) %>%
		rename(model = formula,
					 '&Delta;AICc' = dAICc,
					 coxzph = model_check)

	t1 <- t1 %>%
		kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
		kable_classic(full_width = F, font_size = fontSize) %>%
		add_header_above(c(" ", "Hazard ratio \n (95% confidence interval)" = ncol(t1) - 4, "Model parameters" = 3)) %>%
		footnote(general = paste0("*** < 0.001, ** < 0.01, * < 0.05 & . < 0.01; ",
															"yellow highlighting indicates significance"),
						 general_title = "Note:",
						 footnote_as_chunk = T,
						 title_format = c("italic", "underline"))

	if(style == 'html') return(t1)
	if(style == 'word')	{

				t2 <- t1 %>%
					kableExtra::as_image(zoom = figure_zoom)
				return(t2)
	}
}

make_flextable <-function(df) {

	t1 <- df %>%
	mutate(HR = round(exp(estimate), 3),
				 conf.low = round(exp(conf.low), 3),
				 conf.high = round(exp(conf.high), 3),
				 sig_stars = case_when(p.value < 0.001  ~ "***",
				 											p.value < 0.01  ~ "***",
				 											p.value < 0.05  ~ "*",
				 											p.value < 0.01  ~ ".",
				 											TRUE ~ ""),
				 p.value = round(p.value, 3),
				 sig_stars = paste("^", "**", sig_stars,"**", "^", sep =""),
				 HR = paste(HR, sig_stars, "\n(", conf.low, "-", conf.high,")", sep = "")) %>%
	select(formula, term, HR, AICc, model_check) %>%
	arrange(AICc) %>%
	mutate(min_AICc = min(AICc)) %>%
	mutate(dAICc = AICc - min_AICc) %>%
	select(formula, term, HR, AICc, dAICc, model_check) %>%
	#mutate(term = str_remove(term, "TRUE")) %>%
	mutate(term = map_chr(.x = term, .f = fix_terms)) %>%
	mutate(formula =
				 	str_replace(formula, "maternal_loss", "maternal loss"),
				 formula =
				 	str_replace(formula, "jDSI_M","xxx"),
				 formula =
				 	str_replace(formula, "xxxtop",
				 							"strength of strongest juvenile dyadic strength with top male
                    (including dad; top male jDSI)"),
				 formula =
				 	str_replace(formula, "xxxde_top",
				 							"strength of strongest juvenile dyadic strength with top male
                    (excluding dad; top non-dad jDSI)"),
				 formula =
				 	str_replace(formula, "xxxde","juvenile dyadic strength with top3 males (excluding dad; jDSI_M)"),
				 formula =
				 	str_replace(formula, "xxx","juvenile dyadic strength with top3 males (jDSI_M)"),

				 formula =
				 	str_replace(formula, "cumulative_adversity","cumulative adversity"),
				 formula =
				 	str_replace(formula, "jDSI_paternal","juvenile dyadic strength with dad (jDSI with dad)"),
				 formula =
				 	str_replace(formula, "jSCI_M","xxx"),
				 formula =
				 	str_replace(formula, "xxxde","juvenile social connectedness to males excluding dad (jSCI_Mde)"),
				 formula =
				 	str_replace(formula, "xxx","juvenile social connectedness to males (jSCI_M)"),
				 formula =
				 	str_replace(formula, "dad_overlap_l","proportion of juvenile period shared with dad in same grp"),
				 formula =
				 	str_replace(formula, "dad_overlap","proportions of juvenile period shared with dad in same grp"),
				 formula = paste0('Adult female survival - ', formula)) %>%
	# mutate(HR = ifelse(str_detect(HR, "\\*"),
	# 									 paste('<span style="background-color:Yellow;">',
	# 									 			HR, '</span>', sep = ""),
	# 									 HR)) %>%
	group_by(term) %>%
	pivot_wider(names_from = term, values_from = HR, values_fill = NA_character_) %>%
	select(!(AICc:model_check), AICc, dAICc, model_check) %>%
	mutate(formula = paste0("", formula)) %>%
	rename(model = formula,
				 '&Delta;AICc' = dAICc,
				 coxzph = model_check)

col.index <- matrix(ncol =2, nrow = 0)
#names(col.index) <- c("row", "col")
for(ii in 1:nrow(t1)) {
	cols <- which((str_detect(string = t1[ii, ], pattern = "\\*", negate = FALSE)))
	temp <- as.matrix(tibble(row = rep(ii, times = length(cols)),
													 col = cols))
	col.index <- rbind(col.index, temp)
}

marking = LETTERS[dim(t1)[1] : 1]
models<- t1$model

ft <- t1 %>%
	mutate(model = marking) %>%
	flextable()

for(i in 1:nrow(col.index)) {
	ft <- ft %>%
		bg(i = col.index[i, 1], j = col.index[i, 2], bg = 'yellow')
}

ft <- ft %>%
	add_header_row(values = c("", "Hazard ratio (95% CI)", "Model parameters"), colwidths = c(1, ncol(t1) - 4, 3),
								 top = TRUE) %>%
	colformat_double(j = (ncol(t1)- 2):ncol(t1), digits = 3) %>%
	theme_vanilla()

ii <- dim(t1)[1]
while(ii >= 1) {
	ft <- ft %>%
		add_footer_row(values = paste(marking[ii], "=", c(models[ii])),
									 colwidths = c(ncol(t1)),
									 top = FALSE)
	ii <-ii - 1
}

ft %>%
	colformat_double(digits = 2) %>%
	align(align = "left", part = "all") %>%
	align(align = "center", part = "header") %>%
	fontsize(size = 9, part = "all") %>%
	colformat_md(part = "all") %>%
	fit_to_width(6.5)
#
# 	width(j = 1, width = .1) %>%
#   width(j = (ncol(t1)- 2):ncol(t1), width = .4) %>%
# 	width(j = (2:(ncol(t1)- 3)), width = 1)
}


make_table2 <- function(df, interpretation_text , colored=F, style = 'html', figure_zoom = 1) {
	t1 <- df %>%
		filter(AICc == min(abs(AICc))) %>%
		select(model) %>%
		slice(1) %>%
		mutate(results = map(.x = model, .f = tidy)) %>%
		unnest(cols = c(results)) %>%
		mutate(HR = exp(estimate)) %>%
		select(-model) %>%
		mutate(interpretation = interpretation_text) %>%
		select(term, HR, everything()) %>%
		mutate(term = map_chr(.x = term, .f = fix_terms))

	t1 <- t1 %>%
		kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
		kable_styling(font_size = 15) %>%
		kable_classic(full_width = F)

	if(colored) {
		t1 + row_spec(which(t1$p.value < 0.05), background = "Yellow")
	} else {
		t1
	}

	if(style == 'html') return(t1)
	if(style == 'word')	{

		t2 <- t1 %>%
			kableExtra::as_image(zoom = figure_zoom)
		return(t2)
	}
}

make_flextable2 <- function(df, interpretation_text) {
	t1 <- df %>%
		filter(AICc == min(abs(AICc))) %>%
		select(model) %>%
		slice(1) %>%
		mutate(results = map(.x = model, .f = tidy)) %>%
		unnest(cols = c(results)) %>%
		mutate(HR = exp(estimate)) %>%
		select(-model) %>%
		mutate(interpretation = interpretation_text) %>%
		select(term, HR, everything()) %>%
		mutate(term = map_chr(.x = term, .f = fix_terms))

	t1 %>%
		flextable() %>%
		colformat_double(digits = 2) %>%
		theme_vanilla() %>%
		align(align = "left", part = "all") %>%
		fontsize(size = 9, part = "all") %>%
		autofit() %>%
		fit_to_width(6.5) %>%
		colformat_md(part = "all")


}



make_model_table <- function(model, dataset, explain_text = NA_character_, style = 'html', figure_zoom = 1,highlight = NULL) {
	t1 <- read_csv(paste0('./data/bootstrap_results/', model, "_results.csv")) %>%
		group_by(term) %>%
		summarise_all(.funs = mean) %>%
		mutate(term = map_chr(.x = term, fix_terms)) %>%
		mutate(order = if_else(str_detect(term, "Intercept"), 1, 2)) %>%
		arrange(order, p.value) %>%
		select(-order) %>%
		mutate(interpretation = explain_text) %>%
		kable(format = 'html', booktabs = F, escape=F, digits = 3) %>%
		kable_styling(font_size = 12) %>%
		kable_classic(full_width = F) %>%
		footnote(general = paste0('There are ',
															nrow(dataset), " data points in this analyis with ",
															n_distinct(dataset$AMales), " adult males.")
						 , general_title = ""
						 , footnote_as_chunk = T
						 , title_format = c("italic", "underline"))

	if(!is.null(highlight)) t1 <- t1 %>%  row_spec(row = highlight, background = "Yellow")

	if(style == 'html') return(t1)
	if(style == 'word')	{

		t2 <- t1 %>%
			kableExtra::as_image(zoom = figure_zoom)
		return(t2)
	}
}

FitFlextableToPage <- function(ft, pgwidth = set_pgwidth){

	ft_out <- ft %>% autofit()

	ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
	return(ft_out)
}


make_model_flextable <- function(model, dataset, explain_text = NA_character_, set_width = 6.5) {
	t1 <- read_csv(paste0('./data/bootstrap_results/', model, "_results.csv")) %>%
		group_by(term) %>%
		summarise_all(.funs = mean) %>%
		mutate(term = map_chr(.x = term, fix_terms)) %>%
		mutate(order = if_else(str_detect(term, "Intercept"), 1, 2)) %>%
		arrange(order, p.value) %>%
		select(-order, -VIF) %>%
		mutate(interpretation = explain_text)

	t1 %>%
		flextable() %>%
		colformat_double(digits = 2) %>%
		set_table_properties(layout = "autofit", width = 0) %>%
   	theme_vanilla() %>%
		align(align = "left", part = "all") %>%
		fontsize(size = 9, part = "all") %>%
		rotate(j = 2:6, rotation="btlr",part="header") %>%
		colformat_md(part = "all") %>%
		add_footer_row(values = paste0('There are ',
															nrow(dataset), " data points in this analyis with ",
															n_distinct(dataset$AMales), " adult males."),
									 colwidths = ncol(t1)) %>%
		fit_to_width(set_width)


	#FitFlextableToPage(f1, pgwidth = set_pgwidth)
}








## Overlap plot
## # ggplot() +
# 	geom_segment(data =xdata_ea_temp,
# 							 aes(x = 0, xend= dad_overlap_days,
# 							 		y=ordered_sname,
# 							 		yend = ordered_sname),
# 							 size = .5) +
# 	labs(x = "Duration of time spend in same grp as dad (years)",
# 			 y = "Focal individuals") +
# 	scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
# 	scale_y_continuous(breaks= temp_data_overlap_plot$class_max,
# 										 labels= paste0(temp_data_overlap_plot$class_max, " (",
# 										 							 temp_data_overlap_plot$percentage, '%)'),
# 										 expand = c(0, 0), limits = c(0, NA)) +
# 	geom_segment(data = temp_data_overlap_plot,
# 							 aes(x=class, xend = class,
# 							 		y = 0, yend = class_max),
# 							 linetype = "dashed", color = "red", size = 1) +
# 	geom_segment(data = temp_data_overlap_plot,
# 							 aes(x=0, xend = class,
# 							 		y = class_max, yend = class_max),
# 							 linetype = "dashed", color = "red", size =1) +
# 	# geom_segment(data = temp2,
# 	# 						 aes(x=0, xend = class,
# 	# 				 		     y = class_max, yend = class_max),
# 	# 				         linetype = "dashed", color = "yellow", size =1) +
# 	cowplot::theme_cowplot()
