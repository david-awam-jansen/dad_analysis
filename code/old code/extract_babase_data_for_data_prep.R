library(lubridate)
library(RPostgreSQL)
library(tidyverse)

# Babase data
babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "jansen",
	dbname = "babase",
	password = "Bab00n3455")

## Lets first get the parasite data

## all the samples
master_parasite_l <- tbl(babase, dbplyr::in_schema("gi_parasites", "master_parasite")) %>% 	collect()

## 'all' the samples that have been counted
parasite_total_counts_l <- tbl(babase, dbplyr::in_schema("gi_parasites", "parasite_total_counts")) %>%  collect()

parasite_data <- master_parasite_l %>%
	left_join(parasite_total_counts_l)

parasite_data %>%
	summarise(number_of_samples = n(),
						number_counted = sum(!is.na(count_date)),
						percentage_counted = number_counted/number_of_samples)

## Who are these samples from
biograph_l <- tbl(babase, 'biograph') %>%  collect()

xdata <- parasite_data %>%
	select(sname, sample_date, count_date) %>%
	inner_join(select(biograph_l, sname, sex, birth, statdate)) %>%
	mutate(age_at_sampling = as.numeric((sample_date - birth)/365.25))

## We can make a quick graph
ggplot(data = xdata, aes(age_at_sampling)) +
	 geom_histogram(binwidth = 1)

## critical dates
## Discuss what dates are needed\, but I could think off
## sexual maturation date (F & M)
## age of first kid (F)
## age at with male became ranked (M)
## AGe at with male dispersed (M)

maturedates_l = tbl(babase, "maturedates") %>%
	#inner_join(tbl(babase, 'mstatuses')) %>%
	collect()

rankdates_l = tbl(babase, "rankdates") %>%
	#inner_join(tbl(babase, 'mstatuses')) %>%
		collect()

dispersedates_l <- tbl(babase, "dispersedates") %>%
	#inner_join(select(tbl(babase, 'confidences'), dispconfidence = confidence, everything())) %>%
	collect()


xdata <- xdata %>%
	left_join(maturedates_l, by = "sname") %>%
	left_join(rankdates_l, by = "sname") %>%
	left_join(dispersedates_l, by = "sname")

xdata %>%
	mutate(class = case_when(sex == 'F' & is.na(matured) ~ "JF",
													 sex == 'F' & sample_date < matured ~ "JF",
													 sex == 'F' & sample_date >= matured ~ "AF",
													 sex == 'M' & is.na(matured) ~ "JM",
													 sex == 'M' & sample_date < matured ~ "JM",
													 sex == 'M' & sample_date >= matured ~ "MM")) %>%
	mutate(has_been_counted = if_else(!is.na(count_date) == TRUE, "was counted", "TBD")) %>%
	group_by(sex, has_been_counted, class) %>%
	summarise(counts = n()) %>%
	pivot_wider(names_from = class, values_from = counts)


xdata %>%
	filter(sex == 'M') %>%
	mutate(class = case_when(sex == 'M' & is.na(matured) ~ "JM",
													 sex == 'M' & sample_date < matured ~ "JM",
													 sex == 'M' & sample_date >= matured & is.na(ranked) ~ "SB",
													 sex == 'M' & sample_date >= matured & sample_date < ranked ~ "SB",
				                   sex == 'M' & sample_date >= matured & sample_date >= ranked ~ "AM")) %>%
	mutate(has_been_counted = if_else(!is.na(count_date) == TRUE, "was counted", "TBD")) %>%
	group_by(sex, has_been_counted, class) %>%
	summarise(counts = n()) %>%
	pivot_wider(names_from = class, values_from = counts)

xdata %>%
	filter(sex == 'M') %>%
	mutate(class = case_when(sex == 'M' & is.na(matured) ~ "JM",
													 sex == 'M' & sample_date < matured ~ "JM",
													 sex == 'M' & sample_date >= matured & is.na(ranked) ~ "SB",
													 sex == 'M' & sample_date >= matured & sample_date < ranked ~ "SB",
													 sex == 'M' & sample_date >= matured & sample_date >= ranked ~ "AM")) %>%
	mutate(has_been_counted = if_else(!is.na(count_date) == TRUE, "was counted", "TBD")) %>%
	group_by(sex, has_been_counted, class) %>%
	summarise(counts = n()) %>%
	pivot_wider(names_from = class, values_from = counts)



xdata %>%
	mutate(sample_age = as.numeric(sample_date - birth)/365.25,
				 mature_age = as.numeric(matured - birth)/365.25,
				 ranked_age = as.numeric(ranked - birth)/365.25,
				 disperse_age = as.numeric(dispersed - birth)/365.25) %>%
	group_by(sname) %>%
  mutate(first_sample = min(sample_age)) %>%
	ungroup() %>%
	mutate(sname = forcats::fct_reorder(sname, desc(mature_age))) %>%
	ggplot(aes(y = sname)) +
	  geom_point(aes(x=sample_age)) +
	  geom_point(aes(x=mature_age), shape = 77, color = "red", size = 4) +
	  geom_point(aes(x=ranked_age), shape = 82, color = "red", size = 4) +
	  geom_point(aes(x=disperse_age), shape = 68, color = "red", size = 4) +
		geom_line(aes(x=mature_age, group =1 , color = "matured"), orientation = "y", size = 4) +
	  geom_line(aes(x=ranked_age, group =1 , color = "ranked"), orientation = "y", size = 4) +
	  geom_line(aes(x=disperse_age, group =1 , color = "dispersed"), orientation = "y", size = 4) +
		facet_wrap(. ~ sex, scales = "free_y") +
	  cowplot::theme_cowplot() +
	  theme(legend.position = "bottom")



xdata %>%
	mutate(sample_age = as.numeric(sample_date - birth)/365.25,
				 mature_age = as.numeric(matured - birth)/365.25,
				 ranked_age = as.numeric(ranked - birth)/365.25,
				 disperse_age = as.numeric(dispersed - birth)/365.25) %>%
	group_by(sname) %>%
	select(sname, mature_age) %>%
	distinct() %>%
	ungroup() %>%
	mutate(sname = forcats::fct_reorder(sname, desc(mature_age))) %>%
	ggplot(aes(x=mature_age, y = sname)) +








excluded_members <- selected_behave_gaps %>%
	filter(exclude_group_size == 'YES') %>%
	select(bgid) %>%
	pull()

excluded_ranks <- selected_behave_gaps %>%
	filter(exclude_rank == 'YES') %>%
	select(bgid) %>%
	pull()


current_date <- toupper(format(Sys.Date(), '%d%b%y'))

actor_actees_l <- tbl(babase, "actor_actees") %>%
	filter(act == 'G ') %>%
	inner_join(select(tbl(babase, "parents"), actor = kid, actor_dad = dad)) %>%
	inner_join(select(tbl(babase, "parents"), actee = kid, actee_dad = dad)) %>%
	inner_join(select(tbl(babase, "biograph"),
										actor = sname, actor_birth = birth, actor_sex = sex)) %>%
	inner_join(select(tbl(babase, "biograph"),
										actee = sname, actee_birth = birth, actee_sex = sex)) %>%
	left_join(select(tbl(babase, "rankdates"), actor = sname, actor_ranked = ranked)) %>%
	left_join(select(tbl(babase, "rankdates"), actee = sname, actee_ranked = ranked)) %>%
	mutate(actor_is_ranked = date >= actor_ranked,
				 actee_is_ranked = date >= actee_ranked) %>%
	collect() %>%
	mutate(actor_age = as.numeric((date - actor_birth)/365.25),
				 actee_age = as.numeric((date - actee_birth)/365.25))

## babase dataset
behave_gaps <- tbl(babase, "behave_gaps")
behave_gaps_l <- behave_gaps %>%  collect()

biograph_l <- tbl(babase, "biograph") %>%  collect()  ## biographic data

d_days_l <- tbl(babase, 'members') %>%
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

groups_history <- tbl(babase, "groups_history") %>%
	select(grp = gid, permanent, impermanent, cease_to_exist, last_reg_census) %>%
	group_by(grp) %>%
	mutate(last_date = pmin(impermanent, cease_to_exist, last_reg_census, "2020-07-01")) %>%
	ungroup()

groups_history_l <- groups_history %>%  collect()

groups_l <- tbl(babase, 'groups') %>%  collect()

hybridgene_scores_l <- tbl(babase, "hybridgene_scores") %>%
	filter(hgaid == 2) %>%
	collect()

parents_l <- tbl(babase, "parents") %>%  collect()

potential_dads_l <- tbl(babase, "potential_dads") %>%  collect()## Get babase data

rankdates_l <- tbl(babase, 'rankdates') %>%  collect()

### members and rank datasets controlled for behave gaps
members_partner_overlap_l <- tbl(babase, 'members') %>%
	collect() %>%
	arrange(sname, date)

members <- tbl(babase, 'members') %>%
	filter(grp < 3) %>%
	select(membid, grp, sname, date, grpofresidency)

members <- members %>%
	left_join(groups_history) %>%
	filter(date >= permanent & date <= last_date)

members_bg <- members %>%
	dplyr::left_join(tbl(babase, "behave_gaps"), by = "grp") %>%
	dplyr::filter(date >= gap_start & date <= gap_end) %>%
	select(membid, bgid) %>%
	filter(bgid %in% excluded_members)

members <- members %>%
	anti_join(members_bg, by = 'membid')

members_l <- members %>%
	collect()

ranks <- tbl(babase, 'proportional_ranks') %>%
	filter(grp < 3)

ranks <- ranks %>%
	inner_join(groups_history) %>%
	filter(rnkdate >= permanent & rnkdate <= last_date)

ranks_bg <- ranks %>%
	dplyr::left_join(tbl(babase, "behave_gaps"), by = "grp") %>%
	dplyr::filter(rnkdate >= gap_start & rnkdate <= gap_end) %>%
	select(rnkid, bgid) %>%
	filter(bgid %in% excluded_ranks)

ranks <- ranks %>%
	anti_join(ranks_bg, by = 'rnkid')

## To get it use it in analysis and/or save it you have to 'collect; the data
ranks_l <- ranks %>%
	collect()



## data management of babase datasets
members_AM_l <- members_l %>%  ## Only ranked males
	inner_join(rankdates_l) %>%
	filter(date > ranked) %>%
	select(sname, grp, date, ranked)

members_juveniles_l <- members_l %>%
	inner_join(select(biograph_l, sname, birth), by = 'sname') %>%
	mutate(age = as.numeric(date- birth)/365.25) %>%
	filter(age < 4) %>%
	inner_join(select(parents_l, sname = kid, dad))

ls(pattern = "_l")

save(list = ls(pattern = "_l"),
		 file = paste0("./data/data_set_for_dad_analysis_", current_date, ".RData"))


