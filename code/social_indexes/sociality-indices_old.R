.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("The sociality indices code was last updated on ", Sys.time()))
}

print(paste0("The sociality indices code was last updated on ", Sys.time()))

## general functions
zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
make_NA <- function(x) NA


## SCI support functions
get_mem_dates <- function(my_sub, members_l, df, sel = NULL) {

	mem_dates <- my_sub %>%
		dplyr::ungroup() %>%
		dplyr::inner_join(df, by = c("grp")) %>%
		dplyr::filter(date >= start & date <= end)

	# Remove all rows for dates when the particular animal wasn't present in grp
	remove_rows <- mem_dates %>%
		dplyr::anti_join(members_l, by = c("sname", "grp", "sex_class", "date"))

	# Take set difference and calculate summary
	mem_dates <- mem_dates %>%
		dplyr::setdiff(remove_rows) %>%
		dplyr::select(sname, grp, sex_class, date, !!sel)

	return(mem_dates)
}

get_interaction_dates <- function(my_sub, members_l, df,  my_sex_var,
																	my_role, my_sex=c("F", "M"), my_age) {

	groom_dates <- my_sub %>%
		dplyr::ungroup() %>%
		dplyr::inner_join(df, by = c("sname" = my_role)) %>%
		dplyr::filter(date >= start & date <= end & !!(my_sex_var) %in% my_sex)

	if(my_role == "actor") {
		groom_dates <- groom_dates %>%  filter(is_actee_adult == my_age)
	} else {
		groom_dates <- groom_dates %>%  filter(is_actor_adult == my_age)
	}

	# Remove all rows for dates when the particular animal wasn't present in grp
	remove_rows <- groom_dates %>%
		dplyr::anti_join(members_l, by = c("sname", "grp", "date"))

	# Take set difference and calculate summary
	groom_dates <- groom_dates %>%
		dplyr::setdiff(remove_rows) %>%
		dplyr::select(sname, grp, date, iid)

	return(groom_dates)
}

## SCI functions
get_sci_subset <- function(df, biograph_l, members_l,
													 focals_l, females_l, interactions_l,
													 min_res_days, directional) {

	zero_daily_count <- 1/365.25
	log_zero_daily_count <- log2(zero_daily_count)
	my_sname <- df$sname
	my_sex_class <- df$sex_class

	my_members <- members_l %>%
		filter(sex_class == my_sex_class)

	# Get all members of same sex as the focal animal during relevant time period
	my_subset <- my_members %>%
		dplyr::inner_join(select(df, -sname, -grp), by = c("sex_class")) %>%
		dplyr::filter(date >= start & date <= end) %>%
		dplyr::group_by(sname, grp, sex_class) %>%
		dplyr::summarise(days_present = n(),
										 start = min(date),
										 end = max(date),
										 .groups = "drop") %>%
		ungroup()

	# To allow animals to be non-adults only if focal is not adults
	# Should not  affect results for adults
	if(df$sex_class %in% c("AM", "AF")) {
		my_interactions <- interactions_l %>%
			dplyr::filter((is_actor_adult & is_actee_adult) )
	} else {
		my_interactions <- interactions_l %>%
			dplyr::filter(!(is_actor_adult & is_actee_adult))
	}
	my_interactions <- my_interactions %>%
		dplyr::filter(actor %in% my_subset$sname |
										actee %in% my_subset$sname) %>%
		dplyr::filter(date >= df$start & date <= df$end)

	moms <- biograph_l %>%
		select(sname, pid)  %>%
		mutate(mom = str_sub(pid,1,3))

	dads <- biograph_l %>%
		inner_join(select(parents_l, sname = kid, dad), by = "sname") %>%
		mutate(dad = if_else(dad == "\\N", NA_character_, dad)) %>%
		select(sname, dad)

	my_interactions_mom_excluded <- my_interactions %>%
		filter(actor_sex_class == "AF" | actee_sex_class == "AF") %>%
		left_join(select(moms, actor = sname, actor_mom = mom), by = "actor") %>%
		left_join(select(moms, actee = sname, actee_mom = mom), by = "actee") %>%
		filter(!(actor_sex_class == "JUV" & is.na(actor_mom)) |
					 	!(actee_sex_class == "JUV" & is.na(actee_mom))) %>%
		mutate(maternal_grooms = actor == actee_mom | actee == actor_mom) %>%
		filter(maternal_grooms == FALSE)

	my_interactions_dad_excluded <- my_interactions %>%
		filter(actor_sex_class == "AM" | actee_sex_class == "AM") %>%
		left_join(select(dads, actor = sname, actor_dad = dad), by = "actor") %>%
		left_join(select(dads, actee = sname, actee_dad = dad), by = "actee") %>%
		filter(actor_sex_class == "JUV" | actee_sex_class == "JUV") %>%
		filter(!(actor_sex_class == "JUV" & is.na(actor_dad)) |
					 	!(actee_sex_class == "JUV" & is.na(actee_dad))) %>%
		mutate(paternal_grooms =
					 	case_when(actor_sex_class == "JUV" ~ actee == actor_dad,
					 						actee_sex_class == "JUV" ~ actor == actee_dad)) %>%
		filter(paternal_grooms == FALSE)

	## Focal counts
	# Get all focals during relevant time period in grp
	my_focals <- get_mem_dates(my_subset, my_members, focals_l, sel = quo(sum))

	## Observation days
	# Focal animal was present and at least one focal sample was collected
	obs_days <- my_focals %>%
		group_by(grp, sname, sex_class) %>%
		summarise(days_observed = n(), .groups = "drop")

	my_subset <- my_subset %>%
		left_join(obs_days, by = c("sname", "grp", "sex_class"))

	my_focals <- my_focals %>%
		dplyr::group_by(grp, sname, sex_class) %>%
		dplyr::summarise(n_focals = sum(sum), .groups = "drop")

	## Female counts
	my_females <- get_mem_dates(my_subset, my_members, females_l,
															sel = quo(nr_females)) %>%
		dplyr::group_by(grp, sname, sex_class) %>%
		dplyr::summarise(mean_f_count = mean(nr_females), .groups = "drop")

	# Join back to my_subset to add n_focals column
	my_subset <- my_subset %>%
		dplyr::left_join(my_focals, by = c("grp", "sname", "sex_class")) %>%
		dplyr::left_join(my_females, by = c("grp", "sname", "sex_class"))

	if (nrow(my_subset) == 0 | nrow(my_focals) == 0 | nrow(my_females) == 0) {
		return(tibble::as.tibble(NULL))
	}

	# Filter and calculate variables
	my_subset <- my_subset %>%
		dplyr::filter(days_present >= min_res_days & mean_f_count > 0) %>%
		dplyr::mutate(OE = (n_focals / mean_f_count) / days_present,
									log2OE = log2(OE)) %>%
		dplyr::filter(!is.na(OE))

	## Interactions given to females by each actor of focal's sex
	gg_f <- get_interaction_dates(my_subset, members_l, my_interactions,
																quo(actee_sex), "actor", "F", TRUE) %>%
		dplyr::group_by(grp, sname) %>%
		dplyr::summarise(ItoF = n(), .groups = "drop")

	## Interactions received from females by each actee of focal's sex
	gr_f <- get_interaction_dates(my_subset, members_l, my_interactions,
																quo(actor_sex), "actee", "F", TRUE) %>%
		dplyr::group_by(grp, sname) %>%
		dplyr::summarise(IfromF = n(),  .groups = "drop")

	# Calculate variables for interactions with males only if:
	# - the interactions are grooming AND the focal animal is female OR
	# - the interactions are anything but grooming
	include_males <- my_interactions$act[[1]] != "G" | (my_interactions$act[[1]] == "G" & df$sex_class != "AM")

	if (include_males) {
		## Interactions given to males by each actor of focal's sex
		gg_m <- get_interaction_dates(my_subset, members_l, my_interactions,
																	quo(actee_sex), "actor", "M", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(ItoM = n(), .groups = "drop")

		## Interactions received from males by each actee of focal's sex
		gr_m <- get_interaction_dates(my_subset, members_l, my_interactions,
																	quo(actor_sex), "actee", "M", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(IfromM = n(),  .groups = "drop")
	}

	if(df$sex_class == "JUV") {
		## Interactions given to juveniles by each actor of focal's sex
		gg_j <- get_interaction_dates(my_subset, members_l, my_interactions,
																	quo(actee_sex), "actor", my_age = FALSE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(ItoJ = n(),  .groups = "drop")

		## Interactions received from juveniles by each actee of focal's sex
		gr_j <- get_interaction_dates(my_subset, members_l, my_interactions,
																	quo(actor_sex), "actee", my_age = FALSE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(IfromJ = n(),  .groups = "drop")

		## Interactions given to juveniles by females that are not their mom
		gg_fme <- get_interaction_dates(my_subset, members_l, my_interactions_mom_excluded,
																		quo(actee_sex), "actor", "F", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(ItoFme = n(),  .groups = "drop")

		## Interactions received from females by that are not their mom
		gr_fme <- get_interaction_dates(my_subset, members_l, my_interactions_mom_excluded,
																		quo(actor_sex), "actee", "F", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(IfromFme = n(),  .groups = "drop")

		## Interactions given to juveniles by females that are not their mom
		gg_mde <- get_interaction_dates(my_subset, members_l, my_interactions_dad_excluded,
																		quo(actee_sex), "actor", "M", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(ItoMde = n(),  .groups = "drop")

		## Interactions received from females by that are not their mom
		gr_mde <- get_interaction_dates(my_subset, members_l, my_interactions_dad_excluded,
																		quo(actor_sex), "actee", "M", TRUE) %>%
			dplyr::group_by(grp, sname) %>%
			dplyr::summarise(IfromMde = n(),  .groups = "drop")
	}


	my_subset <- my_subset %>%
		dplyr::left_join(gg_f, by = c("grp", "sname")) %>%
		dplyr::left_join(gr_f, by = c("grp", "sname"))

	my_subset <- my_subset %>%
		tidyr::replace_na(list(ItoF = 0, IfromF = 0))

	my_subset <- my_subset %>%
		mutate(ItF = ItoF + IfromF)

	if (include_males) {
		my_subset <- my_subset %>%
			dplyr::left_join(gg_m, by = c("grp", "sname")) %>%
			dplyr::left_join(gr_m, by = c("grp", "sname"))

		my_subset <- my_subset %>%
			tidyr::replace_na(list(ItoM = 0, IfromM = 0))

		my_subset <- my_subset %>%
			mutate(ItM = ItoM + IfromM)

	}

	if (df$sex_class == 'JUV') {
		my_subset <- my_subset %>%
			dplyr::left_join(gg_j, by = c("grp", "sname")) %>%
			dplyr::left_join(gr_j, by = c("grp", "sname")) %>%
			dplyr::left_join(gg_fme, by = c("grp", "sname")) %>%
			dplyr::left_join(gr_fme, by = c("grp", "sname")) %>%
			dplyr::left_join(gg_mde, by = c("grp", "sname")) %>%
			dplyr::left_join(gr_mde, by = c("grp", "sname"))

		my_subset <- my_subset %>%
			tidyr::replace_na(list(ItoJ = 0, IfromJ = 0)) %>%
			tidyr::replace_na(list(ItoFme = 0, IfromFme = 0)) %>%
			tidyr::replace_na(list(ItoMde = 0, IfromMde = 0))

		my_subset <- my_subset %>%
			mutate(ItJ = ItoJ + IfromJ) %>%
			mutate(ItFme = ItoFme + IfromFme) %>%
			mutate(ItMde = ItoMde + IfromMde)
	}



	# Calculate variables, first for interactions with females only
	my_subset <- my_subset %>%
		dplyr::mutate(ItoF_daily = ItoF / days_present,
									log2ItoF_daily = dplyr::case_when(
										ItoF == 0 ~ log_zero_daily_count,
										TRUE ~ log2(ItoF_daily)),
									IfromF_daily = IfromF / days_present,
									log2IfromF_daily = dplyr::case_when(
										IfromF == 0 ~ log_zero_daily_count,
										TRUE ~ log2(IfromF_daily)),
									ItF_daily = ItF / days_present,
									log2ItF_daily = dplyr::case_when(
										ItF == 0 ~ log_zero_daily_count,
										TRUE ~ log2(ItF_daily)))

	if (include_males) {
		my_subset <- my_subset %>%
			dplyr::mutate(ItoM_daily = ItoM / days_present,
										log2ItoM_daily = dplyr::case_when(
											ItoM == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItoM_daily)),
										IfromM_daily = IfromM / days_present,
										log2IfromM_daily = dplyr::case_when(
											IfromM == 0 ~ log_zero_daily_count,
											TRUE ~ log2(IfromM_daily)),
										ItM_daily = ItM / days_present,
										log2ItM_daily = dplyr::case_when(
											ItM == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItM_daily)))
	}

	if (df$sex_class == "JUV") {
		my_subset <- my_subset %>%
			dplyr::mutate(ItoJ_daily = ItoJ / days_present,
										log2ItoJ_daily = dplyr::case_when(
											ItoJ == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItoJ_daily)),
										IfromJ_daily = IfromJ / days_present,
										log2IfromJ_daily = dplyr::case_when(
											IfromJ == 0 ~ log_zero_daily_count,
											TRUE ~ log2(IfromJ_daily)),
										ItJ_daily = ItJ / days_present,
										log2ItJ_daily = dplyr::case_when(
											ItJ == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItJ_daily))) %>%
			dplyr::mutate(ItoFme_daily = ItoFme / days_present,
										log2ItoFme_daily = dplyr::case_when(
											ItoFme == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItoFme_daily)),
										IfromFme_daily = IfromFme / days_present,
										log2IfromFme_daily = dplyr::case_when(
											IfromFme == 0 ~ log_zero_daily_count,
											TRUE ~ log2(IfromFme_daily)),
										ItFme_daily = ItFme / days_present,
										log2ItFme_daily = dplyr::case_when(
											ItFme == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItFme_daily))) %>%
			dplyr::mutate(ItoMde_daily = ItoMde / days_present,
										log2ItoMde_daily = dplyr::case_when(
											ItoMde == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItoMde_daily)),
										IfromMde_daily = IfromMde / days_present,
										log2IfromMde_daily = dplyr::case_when(
											IfromMde == 0 ~ log_zero_daily_count,
											TRUE ~ log2(IfromMde_daily)),
										ItMde_daily = ItMde / days_present,
										log2ItMde_daily = dplyr::case_when(
											ItMde == 0 ~ log_zero_daily_count,
											TRUE ~ log2(ItMde_daily)))
	}

	my_subset$SCI_F_Dir <- as.numeric(residuals(lm(data = my_subset, log2ItoF_daily ~ log2OE)))
	my_subset$SCI_F_Rec <- as.numeric(residuals(lm(data = my_subset, log2IfromF_daily ~ log2OE)))
	my_subset$SCI_F_tot <- as.numeric(residuals(lm(data = my_subset, log2ItF_daily ~ log2OE)))

	if (include_males) {
		my_subset$SCI_M_Dir <- as.numeric(residuals(lm(data = my_subset, log2ItoM_daily ~ log2OE)))
		my_subset$SCI_M_Rec <- as.numeric(residuals(lm(data = my_subset, log2IfromM_daily ~ log2OE)))
		my_subset$SCI_M_tot <- as.numeric(residuals(lm(data = my_subset, log2ItM_daily ~ log2OE)))
	}

	if (df$sex_class == "JUV") {
		my_subset$SCI_J_Dir <- as.numeric(residuals(lm(data = my_subset, log2ItoJ_daily ~ log2OE)))
		my_subset$SCI_J_Rec <- as.numeric(residuals(lm(data = my_subset, log2IfromJ_daily ~ log2OE)))
		my_subset$SCI_J_tot <- as.numeric(residuals(lm(data = my_subset, log2ItJ_daily ~ log2OE)))
		my_subset$SCI_Fme_Dir <- as.numeric(residuals(lm(data = my_subset, log2ItoFme_daily ~ log2OE)))
		my_subset$SCI_Fme_Rec <- as.numeric(residuals(lm(data = my_subset,
																										 log2IfromFme_daily ~ log2OE)))
		my_subset$SCI_Fme_tot <- as.numeric(residuals(lm(data = my_subset, log2ItFme_daily ~ log2OE)))
		my_subset$SCI_Mde_Dir <- as.numeric(residuals(lm(data = my_subset, log2ItoMde_daily ~ log2OE)))
		my_subset$SCI_Mde_Rec <- as.numeric(residuals(lm(data = my_subset,
																										 log2IfromMde_daily ~ log2OE)))
		my_subset$SCI_Mde_tot <- as.numeric(residuals(lm(data = my_subset, log2ItMde_daily ~ log2OE)))
	}

	if (!directional) {
		my_subset$SCI_F <- (my_subset$SCI_F_Dir + my_subset$SCI_F_Rec) / 2
		if (include_males) {
			my_subset$SCI_M <- (my_subset$SCI_M_Dir + my_subset$SCI_M_Rec) / 2
		}
		if (df$sex_class == "JUV") {
			my_subset$SCI_J <- (my_subset$SCI_J_Dir + my_subset$SCI_J_Rec) / 2
			my_subset$SCI_Fme <- (my_subset$SCI_Fme_Dir + my_subset$SCI_Fme_Rec) / 2
			my_subset$SCI_Mde <- (my_subset$SCI_Mde_Dir + my_subset$SCI_Mde_Rec) / 2
		}
	}

	return(my_subset)
}


sci <- function(my_iyol, biograph_l, members_l, focals_l, females_l, interactions_l,
								min_res_days = 60, parallel = FALSE, ncores = NULL,
								directional = FALSE) {
	ptm <- proc.time()

	# Return an empty tibble if the subset is empty
	if (is.null(my_iyol) |
			!all(names(my_iyol) %in% c("sname", "grp", "start", "end", "days_present", "sex",
																 "birth", "first_start_date", "matured", "ranked","statdate",
																 "birth_dates",  "midpoint", "age_start_yrs", "age_class",
																 "obs_date", "sex_class")) | min_res_days < 0) {
		stop("Problem with input data. Use the 'make_iyol' or 'make_target_df' function to create the input.")
	}

	if (parallel) {
		avail_cores <- detectCores() -1
		if (!is.null(ncores)) {
			if (ncores > avail_cores) {
				message(paste0("Ignoring 'ncores' argument because only ", avail_cores,
											 " cores are available."))
				ncores <- avail_cores
			}
		} else {
			message(paste0("Using all available cores: ", avail_cores,
										 ". Use 'ncores' to specify number of cores to use."))
			ncores <- avail_cores
		}

		cl <- makeCluster(ncores, outfile = "out.log")
		registerDoSNOW(cl)
		clusterExport(cl, list("get_sci_subset",
													 "get_mem_dates",
													 "get_interaction_dates",
													 "parents_l"));
		pb <- txtProgressBar(min = 0, max = nrow(my_iyol), style = 3)
		progress <- function(n) setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
		subset <- foreach(i = 1:nrow(my_iyol), .options.snow = opts,
											.packages = c('tidyverse')) %dopar% {
												get_sci_subset(my_iyol[i, ], biograph_l, members_l, focals_l,
																			 females_l, interactions_l, min_res_days,
																			 directional)
											}
		close(pb)
		stopCluster(cl)
		my_iyol <- add_column(my_iyol, subset)
	} else {
		if (!is.null(ncores)) {
			message("Ignoring 'ncores' argument because 'parallel' set to FALSE.")
		}
		my_iyol$subset <- list(NULL)
		pb <- txtProgressBar(min = 0, max = nrow(my_iyol), style = 3) # Progress bar
		for (i in 1:nrow(my_iyol)) {
			my_iyol[i, ]$subset <- list(get_sci_subset(my_iyol[i, ],
																								 biograph_l,
																								 members_l,
																								 focals_l,
																								 females_l,
																								 interactions_l,
																								 min_res_days = 60,
																								 directional = FALSE))
			setTxtProgressBar(pb, i)
			close(pb)
		}
	}

	sci_focal <- my_iyol %>%
		unnest(cols = c(subset), names_repair = ~ make.unique(.x)) %>%
		mutate(focal = (sname == sname.1 & grp == grp.1)) %>%
		filter(focal) %>%
		select(sname, grp, start, end, contains("SCI_"))

	res <- left_join(my_iyol, sci_focal, by = c("sname", "grp", "start", "end"))

	attr(res, "directional") <- directional

	tdiff <- (proc.time() - ptm)["elapsed"] / 60
	message(paste0("Elapsed time: ", round(tdiff, 3), " minutes (",
								 round(tdiff / 60, 3), ") hours."))

	return(res)
}


## DSI support functions
get_dyadic_subset <- function(df, biograph_l, members_l, focals_l, females_l,
															interactions_l, min_cores_days, within_grp,
															directional) {
	## some functions
	get_resident_dates <- function(focal_sname, focal_grp, focal_sex_class) {
		my_members %>%
			dplyr::filter(sname == focal_sname
										& grp == focal_grp
										& sex_class == focal_sex_class) %>%
			select(date) %>%
			pull()
	}

	get_resident_overlap <- function(list1, list2) intersect(list1, list2)

	get_actor_interactions <- function(focal_sname, focal_grp, focal_sex_class) {
		my_interactions %>%
			arrange(actor, actee) %>%
			dplyr::filter(actor == focal_sname
										& actor_grp == focal_grp
										& actor_sex_class == focal_sex_class)
	}

	get_actee_interactions <- function(focal_sname, focal_grp, focal_sex_class) {
		my_interactions %>%
			arrange(actee, actor) %>%
			dplyr::filter(actee == focal_sname
										& actee_grp == focal_grp
										& actee_sex_class == focal_sex_class)
	}

	# Return total count of focals during co-residence dates
	get_focal_counts <- function(coresidence_dates, focal_grp) {

		res <- my_focals %>%
			dplyr::filter(date %in% coresidence_dates & grp == focal_grp)

		return(sum(res$sum))
	}

	# Return average number of females present in grp during co-residence dates
	get_female_counts <- function(coresidence_dates, focal_grp) {

		res <- my_females %>%
			dplyr::filter(date %in% coresidence_dates & grp == focal_grp)

		return(mean(res$nr_females))
	}

	# Return grooming by actor to actee during co-residence dates
	get_interactions <- function(my_actor, my_actee,
															 focal_grp, partner_grp,
															 focal_sex_class, partner_sex_class) {
		res <- my_interactions %>%
			dplyr::filter(actor == my_actor
										& actor_grp == focal_grp
										& actor_sex_class == focal_sex_class
										& actee == my_actee
										& actee_grp == partner_grp
										& actee_sex_class == partner_sex_class)

		return(nrow(res))
	}

	get_i_given <- function(df, partner, partner_grp, partner_sex_class) {
		sum(which(df$actee == partner
							& df$actee_grp == partner_grp
							& df$actee_sex_class == partner_sex_class))
	}

	get_i_received <- function(df, partner, partner_grp, partner_sex_class) {
		sum(which(df$actor == partner
							& df$actor_grp == partner_grp
							& df$actor_sex_class == partner_sex_class))
	}

	## create subset
	my_grp <- df$grp
	my_sname <- df$sname
	my_start <- df$start
	my_end <- df$end
	my_sex_class <- df$sex_class
	my_row <- df$row_nr

	print(my_row)

	if (my_sex_class == "SM") {
		message("Social indexes are not calculated for sub adult males")
		return(tibble::as_tibble(NULL))
	}

	# Put some subsets in environment for faster performance
	if (within_grp) {
		message("This has not yet been coded for juvenile version")
	} else {
		my_members <- dplyr::filter(members_l, date >= my_start & date <= my_end)
		my_focals <- dplyr::filter(focals_l, date >= my_start & date <= my_end)
		my_females <- dplyr::filter(females_l, date >= my_start & date <= my_end)
		my_interactions <- dplyr::filter(interactions_l, date >= my_start & date <= my_end)

		# This block filters out all interactions with juveniles if focal is adult
		if(df$sex_class %in% c("AF", "AM")) {
			my_members <- my_members %>%
				dplyr::filter(sex_class %in% c("AF", "AM"))
		}

		# This block filters out all interactions with juveniles if focal is adult
		if(df$sex_class %in% c("AF", "AM")) {
			my_interactions <- my_interactions %>%
				dplyr::filter((is_actor_adult & is_actee_adult) )
		} else {
			# This block filters out all interactions with only adults if focal is juvenile
			my_interactions <- my_interactions %>%
				dplyr::filter(!(is_actor_adult & is_actee_adult))
		}

		if (any(map_int(list(my_members, my_females, my_focals, my_interactions),
										nrow) == 0)) {
			return(tibble::as_tibble(NULL))
		}

		# For each animal in each group during this time, calculate:
		# number of days present, first date, last date
		my_subset <- my_members %>%
			dplyr::rename(sname_sex = sex) %>%
			dplyr::rename(sname_sex_class = sex_class) %>%
			dplyr::rename(sname_grp = grp) %>%
			dplyr::group_by(sname, sname_grp, sname_sex, sname_sex_class) %>%
			dplyr::summarise(days_present = n(),
											 start = min(date),
											 end = max(date),
											 .groups = 'drop')
	}

	## Adding recidence and grooming dat to my subset
	my_subset <- my_subset %>%
		mutate(sname_resident_dates =
					 	pmap(.l = list(sname, sname_grp, sname_sex_class),
					 			 .f = get_resident_dates)) %>%
		mutate(sname_actor_interactions =
					 	pmap(.l = list(sname, sname_grp, sname_sex_class),
		 				  	 .f = get_actor_interactions)) %>%
		mutate(sname_actee_interactions =
					 	pmap(.l = list(sname, sname_grp, sname_sex_class),
		 						 .f = get_actee_interactions)) %>%
		unite("temp_sname", sname, sname_sex, sname_sex_class, sname_grp,
					remove = FALSE)

	my_subset_list <-
		my_subset %>%
		#select(temp_sname, sname_resident_dates) %>%
		select(temp_sname, which(sapply(.,class)=="list"))

	my_subset <- my_subset %>%  select_if(negate(is_list))

	# Find all distinct members IN POPULATION between start and end dates
  dyads <- my_subset %>%
		dplyr::mutate(temp_partner = list(my_subset$temp_sname)) %>%
		unnest(cols = c(temp_partner)) %>%
		separate(temp_partner, c("partner", "partner_sex",
														 "partner_sex_class", "partner_grp"),
						 sep = "_", remove = TRUE) %>%
		mutate(partner_grp = as.numeric(partner_grp)) %>%
		filter(sname != partner) %>%
		dplyr::filter(sname_grp == partner_grp)

	## change name of grp column and make the juvenile the partner in all dyads
	## This is to create 'fake' dryads later
  dyads <- dyads %>%
  	filter(sname_sex_class == "JUV") %>%
  	rename_at(vars(starts_with("sname")),
  						~str_replace(., "sname", "sname_temp")) %>%
	  rename_at(vars(starts_with("partner")),
		  				~str_replace(., "partner", "sname")) %>%
	  rename_at(vars(starts_with("sname_temp")),
		  				~str_replace(., "sname_temp", "partner")) %>%
	  mutate(swapped = TRUE) %>%
	  bind_rows(dyads %>%
		  					filter(sname_sex_class != "JUV") %>%
			  				mutate(swapped = FALSE)) %>%
  	unite("temp_sname", sname, sname_sex, sname_sex_class, partner_grp,
  				remove = FALSE) %>%
  	unite("temp_partner", partner, partner_sex, partner_sex_class, partner_grp,
  				remove = FALSE)

  # Note that the step above created duplicated dyads
	# e.g., sname A and partner B, sname B and partner A
	# If DSI is symmetric, these can be removed
	# That's true here, so remove duplicate dyads
	my_subset <- dyads %>%
		dplyr::rowwise() %>%
		dplyr::mutate(tmp = paste(sort(c(temp_sname, temp_partner)),
															collapse = '')) %>%
		dplyr::distinct(tmp, .keep_all = TRUE) %>%
		ungroup()

	## to save memory we can remove some items
	rm(list = c("dyads", "my_members", "my_interactions"))

	# Remove male-male dyads for grooming
	if (interactions_l$act[[1]] == "G") {
		my_subset <- my_subset %>%
			dplyr::filter(!(sname_sex_class == "AM" & partner_sex_class == "AM"))
	}

	# Remove all dyads where not at least one of the two is the correct sex_class
	my_subset <- my_subset %>%
		dplyr::filter(sname_sex_class == df$sex_class |
										partner_sex_class == df$sex_class)

	# Remove all dyads with sub adult males
	my_subset <- my_subset %>%
		dplyr::filter(!(sname_sex_class == "SM" | partner_sex_class == "SM"))

  ## Adding resident and grooming data back for sname and partners (recidence only)
	my_subset <- my_subset %>%
		select(swapped, sname, partner, contains("temp"), everything()) %>%
		inner_join(my_subset_list, by = "temp_sname") %>%
		inner_join(my_subset_list %>%
							 	select(temp_sname, sname_resident_dates) %>%
							 	rename_at(vars(contains("sname")),
							 						~str_replace(., "sname", "partner")),
							 , by = "temp_partner") %>%
		select(-tmp, !(contains("temp"))) ## getting rid off temp colums

  ## getting overlap data between dyads
	my_subset <- my_subset %>%
		mutate(resident_overlap = map2(.x = sname_resident_dates,
																	 .y = partner_resident_dates,
																	 .f = get_resident_overlap)) %>%
		dplyr::mutate(coresidence_days = length(resident_overlap)) %>%
		dplyr::ungroup() %>%
		dplyr::filter(coresidence_days >= min_cores_days)

	## Focal counts
	# Get total count of focals during each dyad's co-residence dates
	my_subset <- my_subset %>%
		dplyr::mutate(n_focals =
										purrr::pmap_dbl(.l = list(resident_overlap, sname_grp),
																		.f = get_focal_counts)) %>%
		dplyr::filter(n_focals > 0)

	## Female counts
	# Get average number of females in group during the dyad's co-residence dates
	my_subset <- my_subset %>%
		dplyr::mutate(n_females =
										purrr::pmap_dbl(.l = list(resident_overlap, sname_grp),
																							.f = get_female_counts)) %>%
		dplyr::filter(n_females > 0)

	# Remove coresidence_dates to make object simpler
	my_subset <- dplyr::select(my_subset, -contains("resident"))

	# Filter and calculate variables
	my_subset <- my_subset %>%
		dplyr::mutate(OE = (n_focals / n_females) / coresidence_days,
									log2OE = log2(OE)) %>%
		dplyr::filter(!is.na(OE))

	if (directional) {
		message("This has not yet been coded for juvenile version")
	} else {
		# Exit if no data meet the criteria, exit without proceeding further
		# Return NULL
		if (nrow(my_subset) == 0) {
			return(tibble::as_tibble(NULL))
			} else {
				my_subset <- my_subset %>%
					mutate(i_given = pmap_dbl(.l = list(sname_actor_interactions,	partner,
																							partner_grp, partner_sex_class),
																		.f = get_i_given),
								 i_received = pmap_dbl(.l = list(sname_actee_interactions, partner,
								 																partner_grp, partner_sex_class),
								 											.f = get_i_received),
								 i_total = i_given + i_received) %>%
					mutate(i_adj = i_total / coresidence_days,
								 log2_i_adj = log2(i_adj))

				# Classify dyads by dyad type ("F-F", "F-M", or "M-M"), and nest by dyad type
				# Since this is not directional, "F-M" and "M-F" are combined into one category: F-M
				my_subset <- my_subset %>%
					dplyr::rowwise() %>%
					dplyr::mutate(dyad = paste(sort(c(sname, partner)), collapse = '-'),
												dyad_type = paste(sort(c(sname_sex_class, partner_sex_class)),
																					collapse = '-')) %>%
					dplyr::ungroup() %>%
					dplyr::group_by(dyad_type) %>%
					tidyr::nest()

				# Fit regression separately for the dyad types and get residuals
				my_subset <- my_subset %>%
					dplyr::mutate(data = purrr::map(data, fit_dyadic_regression)) %>%
					tidyr::unnest(cols = c(data))

				# Reorganize columns
				my_subset <- my_subset %>%
					dplyr::select(sname, sname_sex, sname_sex_class,
												partner, partner_sex, partner_sex_class, everything()) %>%
					dplyr::arrange(sname, sname_sex_class, partner, partner_sex_class)
			}
}
	my_subset <- my_subset %>%
		select_if(negate(is.list))


	return(my_subset)
}


#' Fit dyadic index regression on subset of data
#'
#' @param df A subset of data on which to fit a regression of interactions on observer effort.
#'
#' @return The input data with an additional column containing the regression residuals.
#' @export
#'
#' @examples
fit_dyadic_regression <- function(df) {

	if (all(df$i_adj == 0)) {
		return(tbl_df(NULL))
	}

	# There will be lots of zeros in most subsets
	# If present, remove these before fitting regression model
	zero_subset <- dplyr::filter(df, i_adj == 0)

	if (nrow(zero_subset) > 0) {

		# Fit regression to non-zero values
		nonzero_subset <- dplyr::filter(df, i_adj != 0)
		nonzero_subset$res_i_adj <- as.numeric(residuals(lm(data = nonzero_subset, log2_i_adj ~ log2OE)))

		# Assign a value of -Inf for the residuals of zero values (for calculating quantiles)
		zero_subset$res_i_adj <- -Inf

		# Combine with zero and non-zero subsets
		df <- dplyr::bind_rows(zero_subset, nonzero_subset)
	} else {
		df$res_i_adj <- as.numeric(residuals(lm(data = df, log2_i_adj ~ log2OE)))
	}

	return(df)

}








#' Calculate dyadic index variables from individual-year-of-life data
#'
#' @param my_iyol Individual-year-of-life data.
#' @param biograph_l A local copy of the biograph table
#' @param members_l A subset of members table produced by the function 'subset_members'
#' @param focals_l A subset of focals produced by the function 'subset_focals'
#' @param females_l A subset of female counts produced by the function 'subset_females'
#' @param agonism_l A subset of agonism data produced by the function 'subset_agonism'
#' @param min_cores_days The minimum number of coresidence days needed for dyad to be included. Defaults to 60 days.
#' @param within_grp Logical value indicating whether regressions should be fit relative to dyads in entire population (default) or to dyads within-group
#' @param parallel Logical value indicating whether to process in parallel
#' @param ncores Integer value indicating how many cores to use in parallel processing
#' @param directional Logical value indicating whether to preserve directionality
#'
#' @return The input data with an additional list columsn containing the full DSI subset and the focal DSI variables.
#' @export
#'
#' @examples
dyadic_index <- function(my_iyol, biograph_l, members_l, focals_l, females_l,
                         interactions_l, min_cores_days = 60,
                         within_grp = FALSE, parallel = FALSE,
                         ncores = NULL, directional = FALSE) {

  ptm <- proc.time()

  # Return an empty tibble if the subset is empty
  if (is.null(my_iyol) |
      !all(names(my_iyol) %in% c("sname", "grp", "start", "end", "days_present", "sex",
                                 "birth", "first_start_date", "statdate", "birth_dates",
                                 "midpoint", "age_start_yrs", "age_class", "obs_date",
                                 "matured", "ranked", "sex_class", "row_nr")) |
      min_cores_days < 0) {
    stop("Problem with input data. Use the 'make_iyol' or 'make_target_df' function to create the input.")
  }

  if (parallel) {
    avail_cores <- detectCores() -1
    if (!is.null(ncores)) {
      if (ncores > avail_cores) {
        message(paste0("Ignoring 'ncores' argument because only ", avail_cores,
                       " cores are available."))
        ncores <- avail_cores
      }
    } else {
      message(paste0("Using all available cores: ", avail_cores,
                     ". Use 'ncores' to specify number of cores to use."))
      ncores <- avail_cores
    }

    cl <- makeCluster(ncores, outfile = "")
    registerDoSNOW(cl)
    pb <- txtProgressBar(min = 0, max = nrow(my_iyol), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    clusterExport(cl, list("get_dyadic_subset",
    											 "fit_dyadic_regression",
    											 "get_mem_dates",
    											 "get_interaction_dates",
    											 "parents_l"));
    opts <- list(progress = progress)
    subset <- foreach(i = 1:nrow(my_iyol), .options.snow = opts,
                      .packages = c('tidyverse')) %dopar% {
                        get_dyadic_subset(my_iyol[i, ], biograph_l,
                                          members_l, focals_l, females_l,
                                          interactions_l, min_cores_days,
                                          within_grp, directional)
                      }
    close(pb)
    stopCluster(cl)
    my_iyol <- add_column(my_iyol, subset)
  } else {
    if (!is.null(ncores)) {
      message("Ignoring 'ncores' argument because 'parallel' set to FALSE.")
    }
    my_iyol$subset <- list(NULL)
    pb <- txtProgressBar(min = 0, max = nrow(my_iyol), style = 3) # Progress bar
    for (i in 1:nrow(my_iyol)) {
      my_iyol[i, ]$subset <- list(get_dyadic_subset(my_iyol[i, ], biograph_l,
                                                    members_l, focals_l, females_l,
                                                    interactions_l, min_cores_days,
                                                    within_grp, directional))
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  my_iyol <- my_iyol %>%
  	select(-subset) %>%
    dplyr::mutate(di = purrr::pmap(list(sname, grp, subset), get_focal_index))

  attr(my_iyol, "directional") <- directional

  tdiff <- (proc.time() - ptm)["elapsed"] / 60
  message(paste0("Elapsed time: ", round(tdiff, 3), " minutes (",
                 round(tdiff / 60, 3), ") hours."))

  return(my_iyol)
}

add_paternal <- function(focal, focal_grp, focal_start, df) {
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



#' Obtain dyadic interaction subsets for the animal's year of life
#'
#' @param df One row individual-year-of-life data
#' @param biograph_l A local copy of the biograph table
#' @param members_l A subset of members table produced by the function 'subset_members'
#' @param focals_l A subset of focals produced by the function 'subset_focals'
#' @param females_l A subset of female counts produced by the function 'subset_females'
#' @param interactions_l A subset of interaction data
#' @param min_cores_days The minimum number of coresidence days needed for dyad to be included. Defaults to 60 days.
#' @param within_grp Logical value indicating whether regressions should be fit relative to dyads in entire population (default) or to dyads within-group
#' @param directional Logical value indicating whether to preserve directionality
#'
#' @return The input row with an additional list column containing the subset
#'
#' @examples
# get_dyadic_subset <- function(df, biograph_l, members_l, focals_l, females_l,
#                               interactions_l, min_cores_days, within_grp,
#                               directional) {
#
#   # Find and return all co-residence dates for focal_sname and partner_sname in my_members
#   get_overlap_dates <- function(focal_sname, partner_sname,
#                                 focal_grp, partner_grp,
#                                 focal_sex_class, partner_sex_class) {
#
#
#     focal_dates <- my_members %>%
#       dplyr::filter(sname == focal_sname &
#                       grp == focal_grp &
#                       sex_class == focal_sex_class) %>%
#       dplyr::pull(date)
#
#     partner_dates <- my_members %>%
#       dplyr::filter(sname == partner_sname &
#                       grp == partner_grp &
#                       sex_class == partner_sex_class) %>%
#       dplyr::pull(date)
#
#     overlap_dates <- dplyr::intersect(focal_dates, partner_dates)
#
#     overlap_dates <- my_focals %>%
#       filter(date %in% overlap_dates & grp == focal_grp) %>%
#       pull(date)
#
#     return(overlap_dates)
#   }
#
#
#
#   # Return total count of focals during co-residence dates
#   get_focal_counts <- function(coresidence_dates, focal_grp) {
#
#     res <- my_focals %>%
#       dplyr::filter(date %in% coresidence_dates & grp == focal_grp)
#
#     return(sum(res$sum))
#   }
#
#   # Return average number of females present in grp during co-residence dates
#   get_female_counts <- function(coresidence_dates, focal_grp) {
#
#     res <- my_females %>%
#       dplyr::filter(date %in% coresidence_dates & grp == focal_grp)
#
#     return(mean(res$nr_females))
#   }
#
#   # Return grooming by actor to actee during co-residence dates
#   get_interactions <- function(my_actor, my_actee,
#                                focal_grp, partner_grp,
#                                focal_sex_class, partner_sex_class) {
#
#     res <- my_interactions %>%
#       dplyr::filter(actor == my_actor &
#                     actor_grp == focal_grp &
#                     actor_sex_class == focal_sex_class &
#                     actee == my_actee &
#                     actee_grp == partner_grp &
#                     actee_sex_class == partner_sex_class)
#
#     return(nrow(res))
#   }
#
#   my_grp <- df$grp
#   my_sname <- df$sname
#   my_start <- df$start
#   my_end <- df$end
#   my_sex_class <- df$sex_class
#
#   if (my_sex_class == "SM") {
#     message("Social indexes are not calculated for sub adult males")
#     return(tibble::as_tibble(NULL))
#   }
#
#     # Put some subsets in environment for faster performance
#   if (within_grp) {
#     message("This has not yet been coded for juvenile version")
#     my_members <- dplyr::filter(members_l, grp == my_grp & date >= my_start & date <= my_end)
#
#     # This line allows for the focal animal only to be represented as a non-adult
#     my_members <- dplyr::filter(my_members, (sex == "F" & (date >= matured | sname == my_sname)) |
#                                   (sex == "M" & (date >= ranked | sname == my_sname)))
#
#     my_focals <- dplyr::filter(focals_l, grp == my_grp & date >= my_start & date <= my_end)
#     my_females <- dplyr::filter(females_l, grp == my_grp & date >= my_start & date <= my_end)
#     my_interactions <- interactions_l %>%
#       dplyr::filter((actor_grp == my_grp | actee_grp == my_grp) & date >= my_start & date <= my_end)
#
#     # This line allows for the focal animal only to be represented as a non-adult
#     my_interactions <- my_interactions %>%
#       dplyr::filter((is_actor_adult & is_actee_adult) |
#                       (is_actor_adult & actee == my_sname) |
#                       (is_actee_adult & actor == my_sname))
#
#     # Find all distinct members WITHIN GROUP between start and end dates
#     # For each animal in each group durign this time, calculate:
#     # number of days present, first date, last date
#     my_subset <- my_members %>%
#       dplyr::rename(sname_sex = sex) %>%
#       dplyr::inner_join(select(df, -sname, -sex), by = c("grp")) %>%
#       dplyr::group_by(sname, grp, sname_sex) %>%
#       dplyr::summarise(days_present = n(),
#                        start = min(date),
#                        end = max(date))
#
#     } else {
#     my_members <- dplyr::filter(members_l, date >= my_start & date <= my_end)
#
#     # This block filters out all interactions with juveniles if focal is adult
#     if(df$sex_class %in% c("AF", "AM")) {
#       my_members <- my_members %>%
#         dplyr::filter(sex_class %in% c("AF", "AM"))
#       }
#
#     # This line allows for the focal animal only to be represented as a non-adult
#     # my_members <- dplyr::filter(my_members,
#     #                             (sex == "F" & (date >= matured | sname == my_sname)) |
#     #                             (sex == "M" & (date >= ranked | sname == my_sname)))
#
#     my_focals <- dplyr::filter(focals_l, date >= my_start & date <= my_end)
#     my_females <- dplyr::filter(females_l, date >= my_start & date <= my_end)
#
#     my_interactions <- dplyr::filter(interactions_l, date >= my_start & date <= my_end)
#
#     0# This block filters out all interactions with juveniles if focal is adult
#     if(df$sex_class %in% c("AF", "AM")) {
#       my_interactions <- my_interactions %>%
#         dplyr::filter((is_actor_adult & is_actee_adult) )
#     } else {
#       # This block filters out all interactions with only adults if focal is juvenile
#       my_interactions <- my_interactions %>%
#         dplyr::filter(!(is_actor_adult & is_actee_adult))
#     }
#
#     # This line allows for the focal animal only to be represented as a non-adult
#     # my_interactions <- my_interactions %>%
#     #   dplyr::filter((is_actor_adult & is_actee_adult) |
#     #                   (is_actor_adult & actee == my_sname) |
#     #                   (is_actee_adult & actor == my_sname))
#
#     if (any(map_int(list(my_members, my_females, my_focals, my_interactions),
#                     nrow) == 0)) {
#       return(tibble::as_tibble(NULL))
#     }
#
#     # Find all distinct members IN POPULATION between start and end dates
#     # For each animal in each group durign this time, calculate:
#     # number of days present, first date, last date
#     my_subset <- my_members %>%
#       dplyr::rename(sname_sex = sex) %>%
#       dplyr::rename(sname_sex_class = sex_class) %>%
#       dplyr::group_by(sname, grp, sname_sex, sname_sex_class) %>%
#       dplyr::summarise(days_present = n(),
#                        start = min(date),
#                        end = max(date),
#       								 .groups = 'drop')
#   }
#
#   # For each of these records, find all possible dyad partners
#   # Store as new list column and unnest to expand
#    my_subset <- my_subset %>%
#      unite("temp", sname, sname_sex_class, remove = FALSE)
#
#    temp <-my_subset %>%
#    	dplyr::mutate(resident_dates = purrr::pmap(list(sname, grp,  sname_sex_class),
#    																								get_temp_dates))
#
#
#
#   #
#   dyads <- my_subset %>%
#     dplyr::mutate(partner_temp = list(my_subset$temp)) %>%
#     unnest(cols = partner_temp) %>%
#     separate(partner_temp, c("partner", "partner_sex_class"), remove = FALSE) %>%
#     #select(-temp) %>%
#     filter(sname != partner)
#
#   # Get sex and grp of partner
#   # Remove dyads not in same groups
#   dyads <- dyads %>%
#     dplyr::inner_join(dplyr::select(my_subset,
#                                     partner = sname,
#                                     partner_sex = sname_sex,
#                                     partner_grp = grp),
#                       by = "partner") %>%
#     dplyr::filter(grp == partner_grp)
#
#   # Note that the step above created duplicated dyads
#   # e.g., sname A and partner B, sname B and partner A
#   # If DSI is symmetric, these can be removed
#   # That's true here, so remove duplicate dyads
#   my_subset <- dyads %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(tmp = paste(grp, sort(c(temp, partner_temp)), collapse = '')) %>%
#     dplyr::distinct(tmp, .keep_all = TRUE) %>%
#     dplyr::select(-tmp, -temp, -partner_temp) %>%
#     dplyr::ungroup()
#
#   # Remove male-male dyads for grooming
#   if (interactions_l$act[[1]] == "G") {
#     my_subset <- my_subset %>%
#       dplyr::filter(!(sname_sex_class == "AM" & partner_sex_class == "AM"))
#   }
#
#   # Remove all dyads where not at least one of the two is the correct sex_class
#   my_subset <- my_subset %>%
#       dplyr::filter(sname_sex_class == df$sex_class | partner_sex_class == df$sex_class)
#
#   # Remove all dyads with sub adult males
#   my_subset <- my_subset %>%
#     dplyr::filter(!(sname_sex_class == "SM" | partner_sex_class == "SM"))
#
#   ## Co-residence dates
#   # Find all dates during which focal and partner co-resided in my_grp
#   # Get a count of these dates
#   my_subset <- my_subset %>%
#     dplyr::mutate(coresidence_dates = purrr::pmap(list(sname, partner,
#                                                        grp, partner_grp,
#                                                        sname_sex_class, partner_sex_class),
#                                                   get_overlap_dates)) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(coresidence_days = length(coresidence_dates)) %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(coresidence_days >= min_cores_days)
#
#   ## Focal counts
#   # Get total count of focals during each dyad's co-residence dates
#   my_subset <- my_subset %>%
#     dplyr::mutate(n_focals = purrr::pmap_dbl(list(coresidence_dates, grp),
#                                              get_focal_counts)) %>%
#     dplyr::filter(n_focals > 0)
#
#   ## Female counts
#   # Get average number of females in group during the dyad's co-residence dates
#   my_subset <- my_subset %>%
#     dplyr::mutate(n_females = purrr::pmap_dbl(list(coresidence_dates, grp),
#                                               get_female_counts)) %>%
#     dplyr::filter(n_females > 0)
#
#   # Remove coresidence_dates to make object simpler
#   my_subset <- dplyr::select(my_subset, -coresidence_dates)
#
#   # Filter and calculate variables
#   my_subset <- my_subset %>%
#     dplyr::mutate(OE = (n_focals / n_females) / coresidence_days,
#                   log2OE = log2(OE)) %>%
#     dplyr::filter(!is.na(OE))
#
#   # Exit if no data meet the criteria, exit without proceeding further
#   # Return NULL
#   if (nrow(my_subset) == 0) {
#     return(tibble::as_tibble(NULL))
#   }
#
#   ## Interactions between dyad during the period of interest
#   # Since grooming dates are 1st of month for many years, it does not work
#   # to restrict this to co-resident dates only. Instead, use start and end.
#   # Co-residence is (usually) implied by the grooming interaction.
#   # i_given is behavior given by sname to partner
#   # i_received is behavior received by sname from partner
#    my_subset <- my_subset %>%
#     dplyr::mutate(i_given = purrr::pmap_dbl(list(sname, partner,
#                                                  grp, partner_grp,
#                                                  sname_sex_class, partner_sex_class),
#                                       get_interactions),
#            i_received = purrr::pmap_dbl(list(partner, sname,
#                                               partner_grp, grp,
#                                               partner_sex_class, sname_sex_class),
#                                          get_interactions),
#            i_total = i_given + i_received)
#
#   if (!directional) {
#
#     # Calculate variables
#     my_subset <- my_subset %>%
#       dplyr::mutate(i_adj = i_total / coresidence_days,
#                     log2_i_adj = log2(i_adj))
#
#     # Classify dyads by dyad type ("F-F", "F-M", or "M-M"), and nest by dyad type
#     # Since this is not directional, "F-M" and "M-F" are combined into one category: F-M
#     my_subset <- my_subset %>%
#       dplyr::rowwise() %>%
#        dplyr::mutate(dyad = paste(sort(c(sname, partner)), collapse = '-'),
#                      dyad_type = paste(sort(c(sname_sex_class, partner_sex_class)), collapse = '-')) %>%
#       dplyr::ungroup() %>%
#       dplyr::group_by(dyad_type) %>%
#       tidyr::nest()
#
#     # Fit regression separately for the two dyad types and get residuals
#     my_subset <- my_subset %>%
#       dplyr::mutate(data = purrr::map(data, fit_dyadic_regression)) %>%
#       tidyr::unnest()
#
#     # Reorganize columns
#     my_subset <- my_subset %>%
#       dplyr::select(sname, sname_sex, sname_sex_class, partner, partner_sex, partner_sex_class, everything()) %>%
#       dplyr::arrange(sname, sname_sex_class, partner, partner_sex_class)
#     } else {
#     # Behaviors for which direction needs to be preserved, e.g., agonism
#     message("This has not yet been coded for the juvenile version")
#
#     # Duplicate each row to have one for sname->partner and one for partner->sname:
#     sub_d1 <- my_subset %>%
#       dplyr::rename_at(vars(contains("partner")), funs(sub('partner', 'anim2', .))) %>%
#       dplyr::rename("anim1" = "sname", "anim1_grp" = "grp", "anim1_sex" = "sname_sex") %>%
#       dplyr::select(contains("anim1"), contains("anim2"), everything())
#
#     # Reverse the relevant columns
#     sub_d2 <- my_subset %>%
#       dplyr::rename_at(vars(contains("partner")), funs(sub('partner', 'anim1', .))) %>%
#       dplyr::rename("anim2" = "sname", "anim2_grp" = "grp", "anim2_sex" = "sname_sex",
#              "i_given" = "i_received", "i_received" = "i_given") %>%
#       dplyr::select(contains("anim1"), contains("anim2"), everything())
#
#     # Combine the two sets
#     my_subset <- dplyr::bind_rows(sub_d1, sub_d2)
#
#     # Calculate variables
#     my_subset <- my_subset %>%
#       dplyr::mutate(i_adj = i_given / coresidence_days,
#                     log2_i_adj = log2(i_adj))
#
#     # Directional dyad type: F-F, F-M, M-M, or M-F
#     # Since this is directional, "F-M" and "M-F" are distinguished
#     # Nest by dyad type
#     my_subset <- my_subset %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(dyad = paste(anim1, anim2, sep = '-'),
#                     dyad_type = paste(anim1_sex, anim2_sex, sep = '-')) %>%
#       dplyr::ungroup() %>%
#       dplyr::group_by(dyad_type) %>%
#       tidyr::nest()
#
#     # Fit regression separately for the different dyad types and get residuals
#     my_subset <- my_subset %>%
#       dplyr::mutate(data = purrr::map(data, fit_dyadic_regression)) %>%
#       tidyr::unnest()
#
#     # Reorganize columns
#     my_subset <- my_subset %>%
#       dplyr::select(sname = anim1, grp = anim1_grp, sname_sex = anim1_sex,
#                     partner = anim2, partner_grp = anim2_grp,
#                     partner_sex = anim2_sex, everything()) %>%
#       dplyr::arrange(sname, partner)
#     }
#   return(my_subset)
# }



#' Calculate dyadic variables for focal animal from a dyadic subset
#'
#' @param my_sname
#' @param my_grp
#' @param my_subset
#'
#' @return The input data with an additional list column containing the dyadic variables.
#'
#' @examples
get_focal_index <- function(my_sname, my_grp, my_subset) {

  # Return an empty tibble if the subset is empty
  if (nrow(my_subset) == 0) {
    return(tibble::as_tibble(NULL))
  }

  # Calculate 50 and 90 percentiles from full set of residuals
  # This represents all dyads of a given type in the group during the year
  # NOT including the "zero-interaction" values set to -Inf
  percs <- my_subset %>%
    dplyr::group_by(dyad_type) %>%
    dplyr::filter(i_adj > 0) %>%
    dplyr::summarise(perc_50 = quantile(res_i_adj, probs = 0.5),
                     perc_90 = quantile(res_i_adj, probs = 0.9))

  # Add percentile columns to my_subset
  my_subset <- dplyr::left_join(my_subset, percs, by = "dyad_type")

  # The focal subset should contain only dyads that include my_sname
  # This can be in either the sname or partner column
  # Categorize as very strongly bonded, strongly bonded, weakly bonded, or not bonded
  focal_di <- my_subset %>%
    dplyr::filter(grp == my_grp & (sname == my_sname | partner == my_sname)) %>%
    dplyr::mutate(bond_strength = dplyr::case_when(
      res_i_adj >= perc_90 ~ "VeryStronglyBonded",
      res_i_adj >= perc_50 ~ "StronglyBonded",
      res_i_adj >= -9999999 ~ "WeaklyBonded",
      TRUE ~ "NotBonded"))

  return(focal_di)
}


#' Create summary data from dyadic index input
#'
#' @param df Dyadic index input data produced by the 'dyadic_index' function.
#'
#' @return An individual-year-of-life data set with summarized dyadic variables.
#' @export
#'
#' @examples
dyadic_index_summary <- function(df) {

  # Return an empty tibble if the subset is empty
  if (is.null(df) |
      !all(names(df) %in% c("sname", "grp", "start", "end", "days_present", "sex",
                            "birth", "first_start_date", "ranked", "matured",
                            "statdate", "birth_dates",
                            "midpoint", "age_start_yrs", "age_class", "sex_class",
                            #"subset",
      											"di", "obs_date"))) {
    stop("Problem with input data. Use the 'dyadic_index' function to create the input.")
  }

  directional <- attr(df, "directional")

  df$di_sum <- list(NULL)
  pb <- txtProgressBar(min = 0, max = nrow(df), style = 3) # Progress bar
  for (i in 1:nrow(df)) {
    df[i, ]$di_sum <- list(dyadic_row_summary(df$di[[i]], df$sname[[i]], directional))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  df <- df %>%
    #dplyr::select(-subset, -di) %>%
  	dplyr::select(-di) %>%
    tidyr::unnest(cols = c(di_sum))

  di_strength <- df %>%
    dplyr::select(-top_partners, -r_quantity, -r_reciprocity) %>%
    tidyr::unnest(cols = c(r_strength)) %>%
    dplyr::select(-n)

  di_recip <- df %>%
    dplyr::select(-top_partners, -r_quantity, -r_strength) %>%
    tidyr::unnest(cols = c(r_reciprocity)) %>%
    dplyr::select(-n)

  if (directional) {
    di_strength <- di_strength %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_F",
        sex_class == "AM" & dyad_type == "AM-AM" ~ "DSI_M",
        sex_class == "AF" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "DSI_F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AF", "AF-JUV") ~ "DSI_F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AM", "AM-JUV") ~ "DSI_M",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "DSI_J"),
        sex = forcats::fct_recode(sex, Male = "M", Female = "F"),
        DSI_type = paste(DSI_type, direction, sep = "_")) %>%
      dplyr::select(-dyad_type, -direction) %>%
      tidyr::spread(DSI_type, r_strength) %>%
      dplyr::select(sname, grp, start, end, contains("DSI"))

    di_quantity <- df %>%
      dplyr::select(-top_partners, -r_strength, -r_reciprocity) %>%
      tidyr::unnest(cols = c(r_quantity)) %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type %in% c("AF-AM", "AM-AF") ~ "F",
        sex_class == "AM" & dyad_type == "AM-AM" ~ "M",
        sex_class == "AF" & dyad_type %in% c("AF-AM", "AM-AF") ~ "M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AF", "AF-JUV") ~ "F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AM", "AM-JUV") ~ "M",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "J"),
        sex = forcats::fct_recode(sex, Male = "M", Female = "F"),
        DSI_type = paste(DSI_type, direction, sep = "_")) %>%
      dplyr::select(-dyad_type, -direction) %>%
      tidyr::gather(bond_cat, n_bonds, contains("Bonded")) %>%
      tidyr::unite(var, bond_cat, DSI_type) %>%
      tidyr::spread(var, n_bonds, fill = 0) %>%
      dplyr::select(sname, grp, start, end, ends_with("_Dir"), ends_with("_Rec"))

    di_recip <- di_recip %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type %in% c("AF-AM", "AM-AF") ~ "recip_F",
        sex_class == "AM" & dyad_type == "AM-AM" ~ "recip_M",
        sex_class == "AF" & dyad_type %in% c("AF-AM", "AM-AF") ~ "recip_M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "recip_F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AF", "AF-JUV") ~ "recip_F",
        sex_class == "JUV" & dyad_type %in% c("JUV-AM", "AM-JUV") ~ "recip_M",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "recip_J"),
        sex = forcats::fct_recode(sex, Male = "M", Female = "F"),
        DSI_type = paste(DSI_type, direction, sep = "_")) %>%
      dplyr::select(-dyad_type, -direction) %>%
      tidyr::spread(DSI_type, r_reciprocity) %>%
      dplyr::select(sname, grp, start, end, contains("recip"))
  } else {
    di_strength <- di_strength %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type == "AM-AM" ~ "DSI_M",
        sex_class == "AM" & dyad_type == "AF-AM" ~ "DSI_F",
        sex_class == "AF" & dyad_type == "AF-AM" ~ "DSI_M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "DSI_F",
        sex_class == "JUV" & dyad_type == "AM-JUV" ~ "DSI_M",
        sex_class == "JUV" & dyad_type == "AF-JUV" ~ "DSI_F",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "DSI_J"),
        sex = forcats::fct_recode(sex, Male = "M", Female = "F")) %>%
      dplyr::select(-dyad_type) %>%
      tidyr::spread(DSI_type, r_strength) %>%
      dplyr::select(sname, grp, start, end, one_of("DSI_F", "DSI_M", "DSI_J"))

    di_quantity <- df %>%
      dplyr::select(-top_partners, -r_strength, -r_reciprocity) %>%
      tidyr::unnest(cols = c(r_quantity)) %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type == "AM-AM" ~ "M",
        sex_class == "AM" & dyad_type == "AF-AM" ~ "F",
        sex_class == "AF" & dyad_type == "AF-AM" ~ "M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "F",
        sex_class == "JUV" & dyad_type == "AM-JUV" ~ "M",
        sex_class == "JUV" & dyad_type == "AF-JUV" ~ "F",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "J"),
        sex_class = forcats::fct_recode(sex, Male = "M", Female = "F")) %>%
      dplyr::select(-dyad_type) %>%
      tidyr::gather(bond_cat, n_bonds, contains("Bonded")) %>%
      tidyr::unite(var, bond_cat, DSI_type) %>%
      tidyr::spread(var, n_bonds, fill = 0) %>%
      dplyr::select(sname, grp, start, end,
                    ends_with("_M"), ends_with("_F"), ends_with("_J"))

    di_recip <- di_recip %>%
      dplyr::mutate(DSI_type = case_when(
        sex_class == "AM" & dyad_type == "AM-AM" ~ "recip_M",
        sex_class == "AM" & dyad_type == "AF-AM" ~ "recip_F",
        sex_class == "AF" & dyad_type == "AF-AM" ~ "recip_M",
        sex_class == "AF" & dyad_type == "AF-AF" ~ "recip_F",
        sex_class == "JUV" & dyad_type == "AM-JUV" ~ "recip_M",
        sex_class == "JUV" & dyad_type == "AF-JUV" ~ "recip_F",
        sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "recip_J"),
        sex_class = forcats::fct_recode(sex, Male = "M", Female = "F")) %>%
      dplyr::select(-dyad_type) %>%
      tidyr::spread(DSI_type, r_reciprocity) %>%
      dplyr::select(sname, grp, start, end,
                    one_of("recip_F", "recip_M", "recip_J"))
  }

  di_summary <- df %>%
    dplyr::select(-top_partners, -starts_with("r_")) %>%
    dplyr::left_join(di_strength, by = c("sname", "grp", "start", "end")) %>%
    dplyr::left_join(di_quantity, by = c("sname", "grp", "start", "end")) %>%
    dplyr::left_join(di_recip, by = c("sname", "grp", "start", "end"))

  return(di_summary)
}


#' Summarize a single DSI subset.
#' Calculate top partners, bond quantity, bond strength, and bond reciprocity
#' for each individual-year-of-life separately for the focal animal's top 3
#' grooming partners, separately for F-M and F-F dyads (if applicable)
#'
#' @param df A DSI subset
#' @param focal Focal animal
#' @param directional Logical value indicating whether to preserve directionality
#'
#' @return The input data with additional list columns.
#' @export
#'
#' @examples
dyadic_row_summary <- function(df, focal, directional) {

  # Return an empty tibble if the subset is empty
  if (nrow(df) == 0) {
    return(tibble::as_tibble(NULL))
  }

  # Reciprocity is the mean of interaction asymmetry for the top three partners
  if (directional) {

    df <- df %>%
      dplyr::mutate(direction = dplyr::if_else(sname == focal, "Dir", "Rec"))

    # Relationship quantity is the number of bonds in each bond-strength category
    r_quantity <- df %>%
      dplyr::group_by(dyad_type, direction) %>%
      dplyr::count(bond_strength) %>%
      tidyr::spread(bond_strength, n)

    # Top partners are the top three interaction partners
    # This is calculated separately for each dyad type
    top_partners <- df %>%
      dplyr::filter(res_i_adj > -9999) %>%
      dplyr::arrange(dyad_type, direction, desc(res_i_adj)) %>%
      dplyr::group_by(dyad_type, direction) %>%
      dplyr::slice(1:3)

    # Relationship strength is the mean of the index value for the top three partners
    r_strength <- top_partners %>%
      dplyr::filter(res_i_adj > -9999) %>%
      dplyr::group_by(dyad_type, direction) %>%
      dplyr::summarise(r_strength = mean(res_i_adj, na.rm = TRUE),
                       n = n())

    r_reciprocity <- top_partners %>%
      dplyr::mutate(recip = (i_received - i_given) / (i_given + i_received)) %>%
      dplyr::group_by(dyad_type, direction) %>%
      dplyr::summarise(r_reciprocity = mean(recip, na.rm = TRUE),
                       n = n())
  } else {
    # Relationship quantity is the number of bonds in each bond-strength category
    r_quantity <- df %>%
      dplyr::group_by(dyad_type) %>%
      dplyr::count(bond_strength) %>%
      tidyr::spread(bond_strength, n)

    # Top partners are the top three interaction partners
    # This is calculated separately for each dyad type
    top_partners <- df %>%
      dplyr::filter(res_i_adj > -9999) %>%
      dplyr::arrange(dyad_type, desc(res_i_adj)) %>%
      dplyr::group_by(dyad_type) %>%
      dplyr::slice(1:3)

    # Relationship strength is the mean of the index value for the top three partners
    r_strength <- top_partners %>%
      dplyr::filter(res_i_adj > -9999) %>%
      dplyr::group_by(dyad_type) %>%
      dplyr::summarise(r_strength = mean(res_i_adj, na.rm = TRUE),
                       n = n())

    r_reciprocity <- top_partners %>%
      dplyr::mutate(recip = 1 - abs((i_given - i_received) / (i_given + i_received))) %>%
      dplyr::group_by(dyad_type) %>%
      dplyr::summarise(r_reciprocity = mean(recip, na.rm = TRUE),
                       n = n())
  }

  res <- tibble(top_partners = list(top_partners),
                r_quantity = list(r_quantity),
                r_strength = list(r_strength),
                r_reciprocity = list(r_reciprocity))

  return(res)
}

# get_dyadic_subset2 <- function(df, biograph_l, members_l, focals_l, females_l,
# 															 interactions_l, min_cores_days, within_grp,
# 															 directional) {
#
# 	print("Using the updated algorith")
# 	my_grp <- df$grp
# 	my_sname <- df$sname
# 	my_start <- df$start
# 	my_end <- df$end
# 	my_sex_class <- df$sex_class
#
# 	if (my_sex_class == "SM") {
# 		message("Social indexes are not calculated for sub adult males")
# 		return(tibble::as_tibble(NULL))
# 	}
#
# 	# Put some subsets in environment for faster performance
# 	if (within_grp) {
# 		message("This has not yet been coded for juvenile version")
# 		my_members <- dplyr::filter(members_l, grp == my_grp & date >= my_start & date <= my_end)
#
# 		# This line allows for the focal animal only to be represented as a non-adult
# 		my_members <- dplyr::filter(my_members, (sex == "F" & (date >= matured | sname == my_sname)) |
# 																	(sex == "M" & (date >= ranked | sname == my_sname)))
#
# 		my_focals <- dplyr::filter(focals_l, grp == my_grp & date >= my_start & date <= my_end)
# 		my_females <- dplyr::filter(females_l, grp == my_grp & date >= my_start & date <= my_end)
# 		my_interactions <- interactions_l %>%
# 			dplyr::filter((actor_grp == my_grp | actee_grp == my_grp) & date >= my_start & date <= my_end)
#
# 		# This line allows for the focal animal only to be represented as a non-adult
# 		my_interactions <- my_interactions %>%
# 			dplyr::filter((is_actor_adult & is_actee_adult) |
# 											(is_actor_adult & actee == my_sname) |
# 											(is_actee_adult & actor == my_sname))
#
# 		# Find all distinct members WITHIN GROUP between start and end dates
# 		# For each animal in each group durign this time, calculate:
# 		# number of days present, first date, last date
# 		my_subset <- my_members %>%
# 			dplyr::rename(sname_sex = sex) %>%
# 			dplyr::inner_join(select(df, -sname, -sex), by = c("grp")) %>%
# 			dplyr::group_by(sname, grp, sname_sex) %>%
# 			dplyr::summarise(days_present = n(),
# 											 start = min(date),
# 											 end = max(date))
#
# 	} else {
# 		my_members <- dplyr::filter(members_l, date >= my_start & date <= my_end)
#
# 		# This block filters out all interactions with juveniles if focal is adult
# 		if(df$sex_class %in% c("AF", "AM")) {
# 			my_members <- my_members %>%
# 				dplyr::filter(sex_class %in% c("AF", "AM"))
# 		}
#
# 		# This line allows for the focal animal only to be represented as a non-adult
# 		# my_members <- dplyr::filter(my_members,
# 		#                             (sex == "F" & (date >= matured | sname == my_sname)) |
# 		#                             (sex == "M" & (date >= ranked | sname == my_sname)))
#
# 		my_focals <- dplyr::filter(focals_l, date >= my_start & date <= my_end)
# 		my_females <- dplyr::filter(females_l, date >= my_start & date <= my_end)
#
# 		my_interactions <- dplyr::filter(interactions_l, date >= my_start & date <= my_end)
#
# 		0# This block filters out all interactions with juveniles if focal is adult
# 		if(df$sex_class %in% c("AF", "AM")) {
# 			my_interactions <- my_interactions %>%
# 				dplyr::filter((is_actor_adult & is_actee_adult) )
# 		} else {
# 			# This block filters out all interactions with only adults if focal is juvenile
# 			my_interactions <- my_interactions %>%
# 				dplyr::filter(!(is_actor_adult & is_actee_adult))
# 		}
#
# 		# This line allows for the focal animal only to be represented as a non-adult
# 		# my_interactions <- my_interactions %>%
# 		#   dplyr::filter((is_actor_adult & is_actee_adult) |
# 		#                   (is_actor_adult & actee == my_sname) |
# 		#                   (is_actee_adult & actor == my_sname))
#
# 		if (any(map_int(list(my_members, my_females, my_focals, my_interactions),
# 										nrow) == 0)) {
# 			return(tibble::as_tibble(NULL))
# 		}
#
# 		# Find all distinct members IN POPULATION between start and end dates
# 		# For each animal in each group durign this time, calculate:
# 		# number of days present, first date, last date
# 		my_subset <- my_members %>%
# 			dplyr::rename(sname_sex = sex) %>%
# 			dplyr::rename(sname_sex_class = sex_class) %>%
# 			dplyr::group_by(sname, grp, sname_sex, sname_sex_class) %>%
# 			dplyr::summarise(days_present = n(),
# 											 start = min(date),
# 											 end = max(date),
# 											 .groups = 'drop')
# 	}
#
# 	# For each of these records, find all possible dyad partners
# 	# Store as new list column and unnest to expand
#
#
# 	## some fuctions
# 	get_resident_dates <- function(focal_sname, focal_grp, focal_sex_class) {
# 		my_members %>%
# 			dplyr::filter(sname == focal_sname
# 										& grp == focal_grp
# 										& sex_class == focal_sex_class) %>%
# 			select(date) %>%
# 			pull()
# 	}
#
# 	get_resident_overlap <- function(list1, list2) intersect(list1, list2)
#
# 	get_actor_interactions <- function(focal_sname, focal_grp, focal_sex_class) {
# 		my_interactions %>%
# 			dplyr::filter(actor == focal_sname
# 										& actor_grp == focal_grp
# 										& actor_sex_class == focal_sex_class)
# 	}
#
# 	get_actee_interactions <- function(focal_sname, focal_grp, focal_sex_class) {
# 		my_interactions %>%
# 			dplyr::filter(actee == focal_sname
# 										& actee_grp == focal_grp
# 										& actee_sex_class == focal_sex_class)
# 	}
#
# 	# Return total count of focals during co-residence dates
# 	get_focal_counts <- function(coresidence_dates, focal_grp) {
#
# 		res <- my_focals %>%
# 			dplyr::filter(date %in% coresidence_dates & grp == focal_grp)
#
# 		return(sum(res$sum))
# 	}
#
# 	# Return average number of females present in grp during co-residence dates
# 	get_female_counts <- function(coresidence_dates, focal_grp) {
#
# 		res <- my_females %>%
# 			dplyr::filter(date %in% coresidence_dates & grp == focal_grp)
#
# 		return(mean(res$nr_females))
# 	}
#
# 	# Return grooming by actor to actee during co-residence dates
# 	get_interactions <- function(my_actor, my_actee,
# 															 focal_grp, partner_grp,
# 															 focal_sex_class, partner_sex_class) {
# 		res <- my_interactions %>%
# 			dplyr::filter(actor == my_actor
# 										& actor_grp == focal_grp
# 										& actor_sex_class == focal_sex_class
# 										& actee == my_actee
# 										& actee_grp == partner_grp
# 										& actee_sex_class == partner_sex_class)
#
# 		return(nrow(res))
# 	}
#
#
#
# 	## my subset stuff
# 	my_subset <- my_subset %>%
# 		mutate(resident_dates = pmap(.l = list(sname, grp, sname_sex_class), .f = get_resident_dates)) %>%
# 		mutate(actor_interactions = pmap(.l = list(sname, grp, sname_sex_class), .f = get_actor_interactions)) %>%
# 		mutate(actee_interactions = pmap(.l = list(sname, grp, sname_sex_class), .f = get_actee_interactions)) %>%
# 		unite("temp", sname, sname_sex_class, remove = FALSE)
#
# 	dyads <- my_subset %>%
# 		dplyr::mutate(partner_temp = list(my_subset$temp)) %>%
# 		unnest(cols = partner_temp) %>%
# 		separate(partner_temp, c("partner", "partner_sex_class"), remove = FALSE) %>%
# 		filter(sname != partner)
#
# 	# Get sex and grp of partner
# 	# Remove dyads not in same groups
# 	dyads <- dyads %>%
# 		dplyr::inner_join(dplyr::select(my_subset,
# 																		partner = sname,
# 																		partner_sex = sname_sex,
# 																		partner_grp = grp,
# 																		parnter_resident_dates = resident_dates ),
# 											by = "partner") %>%
# 		dplyr::filter(grp == partner_grp)
#
# 	# Note that the step above created duplicated dyads
# 	# e.g., sname A and partner B, sname B and partner A
# 	# If DSI is symmetric, these can be removed
# 	# That's true here, so remove duplicate dyads
# 	my_subset <- dyads %>%
# 		dplyr::rowwise() %>%
# 		dplyr::mutate(tmp = paste(grp, sort(c(temp, partner_temp)), collapse = '')) %>%
# 		dplyr::distinct(tmp, .keep_all = TRUE) %>%
# 		dplyr::select(-tmp, -temp, -partner_temp) %>%
# 		dplyr::ungroup()
#
# 	# Remove male-male dyads for grooming
# 	if (interactions_l$act[[1]] == "G") {
# 		my_subset <- my_subset %>%
# 			dplyr::filter(!(sname_sex_class == "AM" & partner_sex_class == "AM"))
# 	}
#
# 	# Remove all dyads where not at least one of the two is the correct sex_class
# 	my_subset <- my_subset %>%
# 		dplyr::filter(sname_sex_class == df$sex_class | partner_sex_class == df$sex_class)
#
# 	# Remove all dyads with sub adult males
# 	my_subset <- my_subset %>%
# 		dplyr::filter(!(sname_sex_class == "SM" | partner_sex_class == "SM"))
#
# 	my_subset <- my_subset %>%
# 		mutate(resident_overlap = map2(.x = resident_dates,
# 																	 .y = parnter_resident_dates,
# 																	 .f = get_resident_overlap)) %>%
# 		dplyr::mutate(coresidence_days = length(resident_overlap)) %>%
# 		dplyr::ungroup() %>%
# 		dplyr::filter(coresidence_days >= min_cores_days)
#
# 	## Focal counts
# 	# Get total count of focals during each dyad's co-residence dates
# 	my_subset <- my_subset %>%
# 		dplyr::mutate(n_focals = purrr::pmap_dbl(.l = list(resident_overlap, grp),
# 																						 .f = get_focal_counts)) %>%
# 		dplyr::filter(n_focals > 0)
#
# 	## Female counts
# 	# Get average number of females in group during the dyad's co-residence dates
# 	my_subset <- my_subset %>%
# 		dplyr::mutate(n_females = purrr::pmap_dbl(.l = list(resident_overlap, grp),
# 																							.f = get_female_counts)) %>%
# 		dplyr::filter(n_females > 0)
#
# 	# Remove coresidence_dates to make object simpler
# 	my_subset <- dplyr::select(my_subset, -resident_overlap)
#
# 	# Filter and calculate variables
# 	my_subset <- my_subset %>%
# 		dplyr::mutate(OE = (n_focals / n_females) / coresidence_days,
# 									log2OE = log2(OE)) %>%
# 		dplyr::filter(!is.na(OE))
#
# 	# Exit if no data meet the criteria, exit without proceeding further
# 	# Return NULL
# 	if (nrow(my_subset) == 0) {
# 		return(tibble::as_tibble(NULL))
# 	}
#
# 	## Interactions between dyad during the period of interest
# 	# Since grooming dates are 1st of month for many years, it does not work
# 	# to restrict this to co-resident dates only. Instead, use start and end.
# 	# Co-residence is (usually) implied by the grooming interaction.
# 	# i_given is behavior given by sname to partner
# 	# i_received is behavior received by sname from partner
#
# 	get_i_given <- function(df, partner, partner_grp, partner_sex_class) {
# 		df %>%
# 			filter(actee == partner
# 						 & actee_grp == partner_grp
# 						 & actee_sex_class == partner_sex_class) %>%
# 			nrow()
# 	}
#
# 	get_i_received <- function(df, partner, partner_grp, partner_sex_class) {
# 		df %>%
# 			filter(actor == partner
# 						 & actor_grp == partner_grp
# 						 & actor_sex_class == partner_sex_class) %>%
# 			nrow()
# 	}
#
# 	# my_subset <- my_subset %>%
# 	# 	dplyr::mutate(i_given = purrr::pmap_dbl(list(sname, partner,
# 	# 																							 grp, partner_grp,
# 	# 																							 sname_sex_class, partner_sex_class),
# 	# 																					get_interactions),
# 	# 								i_received = purrr::pmap_dbl(list(partner, sname,
# 	# 																									partner_grp, grp,
# 	# 																									partner_sex_class, sname_sex_class),
# 	# 																						 get_interactions),
# 	# 								i_total = i_given + i_received)
#
# 	my_subset <- my_subset %>%
# 		mutate(i_given = pmap_dbl(.l = list(actor_interactions, partner, partner_grp, partner_sex_class),
# 															.f = get_i_given),
# 					 i_received = pmap_dbl(.l = list(actee_interactions, partner, partner_grp, partner_sex_class),
# 					 											.f = get_i_received),
# 					 i_total = i_given + i_received)
#
# 	if (!directional) {
# 		# Calculate variables
# 		my_subset <- my_subset %>%
# 			dplyr::mutate(i_adj = i_total / coresidence_days,
# 										log2_i_adj = log2(i_adj))
#
# 		# Classify dyads by dyad type ("F-F", "F-M", or "M-M"), and nest by dyad type
# 		# Since this is not directional, "F-M" and "M-F" are combined into one category: F-M
# 		my_subset <- my_subset %>%
# 			dplyr::rowwise() %>%
# 			dplyr::mutate(dyad = paste(sort(c(sname, partner)), collapse = '-'),
# 										dyad_type = paste(sort(c(sname_sex_class, partner_sex_class)), collapse = '-')) %>%
# 			dplyr::ungroup() %>%
# 			dplyr::group_by(dyad_type) %>%
# 			tidyr::nest()
#
# 		# Fit regression separately for the two dyad types and get residuals
# 		my_subset <- my_subset %>%
# 			dplyr::mutate(data = purrr::map(data, fit_dyadic_regression)) %>%
# 			tidyr::unnest()
#
# 		# Reorganize columns
# 		my_subset <- my_subset %>%
# 			dplyr::select(sname, sname_sex, sname_sex_class, partner, partner_sex, partner_sex_class, everything()) %>%
# 			dplyr::arrange(sname, sname_sex_class, partner, partner_sex_class)
# 	} else {
# 		# Behaviors for which direction needs to be preserved, e.g., agonism
# 		message("This has not yet been coded for the juvenile version")
#
# 		# Duplicate each row to have one for sname->partner and one for partner->sname:
# 		sub_d1 <- my_subset %>%
# 			dplyr::rename_at(vars(contains("partner")), funs(sub('partner', 'anim2', .))) %>%
# 			dplyr::rename("anim1" = "sname", "anim1_grp" = "grp", "anim1_sex" = "sname_sex") %>%
# 			dplyr::select(contains("anim1"), contains("anim2"), everything())
#
# 		# Reverse the relevant columns
# 		sub_d2 <- my_subset %>%
# 			dplyr::rename_at(vars(contains("partner")), funs(sub('partner', 'anim1', .))) %>%
# 			dplyr::rename("anim2" = "sname", "anim2_grp" = "grp", "anim2_sex" = "sname_sex",
# 										"i_given" = "i_received", "i_received" = "i_given") %>%
# 			dplyr::select(contains("anim1"), contains("anim2"), everything())
#
# 		# Combine the two sets
# 		my_subset <- dplyr::bind_rows(sub_d1, sub_d2)
#
# 		# Calculate variables
# 		my_subset <- my_subset %>%
# 			dplyr::mutate(i_adj = i_given / coresidence_days,
# 										log2_i_adj = log2(i_adj))
#
#
# 		# Directional dyad type: F-F, F-M, M-M, or M-F
# 		# Since this is directional, "F-M" and "M-F" are distinguished
# 		# Nest by dyad type
# 		my_subset <- my_subset %>%
# 			dplyr::rowwise() %>%
# 			dplyr::mutate(dyad = paste(anim1, anim2, sep = '-'),
# 										dyad_type = paste(anim1_sex, anim2_sex, sep = '-')) %>%
# 			dplyr::ungroup() %>%
# 			dplyr::group_by(dyad_type) %>%
# 			tidyr::nest()
#
# 		# Fit regression separately for the different dyad types and get residuals
# 		my_subset <- my_subset %>%
# 			dplyr::mutate(data = purrr::map(data, fit_dyadic_regression)) %>%
# 			tidyr::unnest()
#
# 		# Reorganize columns
# 		my_subset <- my_subset %>%
# 			dplyr::select(sname = anim1, grp = anim1_grp, sname_sex = anim1_sex,
# 										partner = anim2, partner_grp = anim2_grp,
# 										partner_sex = anim2_sex, everything()) %>%
# 			dplyr::arrange(sname, partner)
# 	}
# 	return(my_subset)
# }
