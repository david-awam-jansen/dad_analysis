setwd("~/DSI_paternal_female_survival")


## functions
get_kid_data <- function(focal_sname, focal_birth, focal_grp, focal_mom, period) {

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

get_male_present_count <-  function(focal_kid, focal_grp, start, end) {

	start_date <- ymd(start)
	end_date <- ymd(end)

	members_l_AM %>%
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

nr_draws = 1000

step_previous<- who_grooms5 %>%
	filter(previous_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "previous") %>%
	mutate(get_previous = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
														 .f = get_kid_data)) %>%
	unnest(cols= c(get_previous), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														 .f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(previous_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																 .f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, previous_random_draws)

step_next <- who_grooms5 %>%
	filter(next_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "next") %>%
	mutate(get_next = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
														 .f = get_kid_data)) %>%
	unnest(cols= c(get_next), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(next_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																			.f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, next_random_draws)

step_combined <- who_grooms5 %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	left_join(step_previous) %>%
  left_join(step_next)

who_grooms_bootstrap <- step_combined %>%
  mutate(data = pmap(.l = list(data, previous_random_draws, next_random_draws),
  									 .f = assign_random_draws)) %>%
	unnest(cols = c(data)) %>%
	select(-previous_random_draws, -next_random_draws) %>%
	filter(!is.na(mean_proprank)) %>%
	filter(nr_days >= 30)

write_csv(who_grooms_bootstrap, "./data/who_grooms_bootstrap.csv")

# move who_grooms_bootstrap.csv to cluster and run cluster_bootstrap_previous_next_code.R
# and get the results back here

getwd()
result_files <- list.files(path = './data/bootstrap_results/bootstrap_models')
result_files <- result_files[str_detect(result_files, "full") |str_detect(result_files, "best")]


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








write_csv(x = xdata, file = 'full_model_results.csv',
					)



