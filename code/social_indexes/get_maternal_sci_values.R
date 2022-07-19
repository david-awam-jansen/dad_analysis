## Setup R
##
## This overview is based on an example that was created by Fernando Campos.
## I edited to make it work for the cluster and maternal SCI.
##
## This is already done on the cluster
## Note that tro access the cluster you need an acccount with CRC and have a
## VPN connection to Notre Dame (vpn.nd.edu).
##
## I normally use a combination of the cluster and desktop to get
## the social index calculations done. Especially the DSI parts take a lot of time.
## If you only need SCI values it is likely feasible to do it on just a desktop.

# Install devtools if not already present
if (!("devtools" %in% installed.packages()[,"Package"]))install.packages("devtools")

# Install newest version of ramboseli
devtools::install_github(force = TRUE, repo = "amboseli/ramboseli")

# Install packages that are needed
# Note that RPostgreSQL does not work on cluster

Sys.setenv(TZ = 'UTC')
list.of.packages <- list("foreach", "doSNOW", "parallel", "tidyverse",
												 "lubridate", "dbplyr", "purrrlyr", "RPostgreSQL",
												 "zoo", "ramboseli")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages, require, character.only = T)

## The next part won't work on the cluster.
## It has to be done on a desktop as you need to create a vpn connection to Duke.
## I normally download all the data I need and then upload it to cluster
## See below for details

babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "jansen",
	dbname = "babase",
	password = rstudioapi::askForPassword("Database password"))

# Get local copy of biograph table
biograph_l <- collect(tbl(babase, "biograph"))

# Make a members subset that excludes behavioral observation gaps
members_l <- subset_members(babase)

# Subset other data sets used for sociality indices
focals_l <- subset_focals(babase, members_l)
females_l <- subset_females(members_l)

# Grooming
grooming_l <- subset_interactions(babase, members_l, my_acts = c("G"))

# Make an individual-year-of-life data set for adults
# Normally SCI ect are calculated per individual for every year of its life (from birthday to birthday) after it becomes and adult
# We will filter some of these out later
iyol_adults <- make_iyol(babase, members_l, focals_l, grooming_l)

## Now we will create something that will be needed for the maternal SCI
## Maternal SCI is a bit different then the normal approach as we won't go
## from birthday to birthday of the focal (mother).
## Rather it is based on the birthday of the offsrping.
## Not all steps are really used, but they are needed to link with the iyol structure

kid_data <-  biograph_l %>%
	mutate(mom = stringr::str_sub(pid, start = 1, end =3)) %>%  ## extract sname of mom from pid
	## Next I filter out some kids. You could alter this to fit your resrictions
	filter(!is.na(mom)   ## no unknown moms
				 & sex != 'U'   ## Only kids with know sex
				 & statdate > birth)  ## ONly kids that lived at least 1 year)

## There is one twin in the dataset
## This will create duplicate issues.
## I will throw one of the twin out
kid_data <- kid_data %>%
	filter(sname != 'PIL')

maternal_step1 <- kid_data %>%
	select(kid = sname, sname = mom, birth) %>% ## mom needs to become the focal = sname
	mutate(kid_birth = birth) %>%  ## We want to keep kid_birth for a later stage
	group_by(sname, kid) %>%
	## We nee to calculate social indexes per year so we create two rows
	## 1st from birth of kid to its first birthday; 2nd from 1st birthday to 2nd birthday
	mutate(start = list(seq(from = kid_birth, to = kid_birth + years(2) - days(1), by = 'years'))) %>%
	unnest(cols = c(start)) %>%
	ungroup()

## Combinations of sname (mom), start and grp need to be unique.
## There are a few small issues
## These are caused by siblings being born on the birthday of a previous sibling
maternal_step1 %>%
	group_by(sname, start) %>%
	arrange(sname, start) %>%
	filter(n() != 1) %>%
	mutate(temp = 1:n())

## Lets filter 1 of the duplicates out

maternal_step2 <- maternal_step1 %>%
	group_by(sname, start) %>%
	mutate(temp = 1:n()) %>%
	filter(temp != 2) %>%
	select(-temp) %>%
	mutate(end = start + years(1) - days(1))  %>%
	left_join(select(tbl(babase, 'maturedates'), sname, matured) %>%  collect(), by = 'sname') %>%
	left_join(select(tbl(babase, 'rankdates'), sname, ranked) %>%  collect(), by = 'sname') %>%
	left_join(select(tbl(babase, 'biograph'), sname, statdate, sex) %>%  collect(), by = 'sname') %>%
	mutate(first_start_date = if_else(sex == 'F', matured, ranked)) %>%
	mutate(midpoint = start + floor((end - start)/2),
				 age_start_yrs = as.numeric(start - birth)/365.25,
				 age_class = floor(plyr::round_any(age_start_yrs, 0.005)) + 1) %>%
	arrange(sname, kid, birth, start) %>%
	select(-matured, -ranked)

## The next step is to check in which group the mother was and for how many days.

maternal_iyol <- maternal_step2 %>%
	inner_join(maternal_step2 %>%
						 	left_join(maternal_step2 %>%
						 							left_join(select(members_l, sname, date, grp), by = c("sname")) %>%
						 							filter(date >= start & date <= end) %>%
						 							distinct(sname, start, grp)) %>%
						 	inner_join(select(members_l, sname, grp, date), by = c("sname", "grp")) %>%
						 	filter(date >= start & date <= end) %>%
						 	group_by(sname, grp, start, end) %>%
						 	summarise(days_present = n(), .groups = 'drop') %>%
						 	arrange(sname, grp, start, end),
						 by = c("sname", "start", "end")) %>%
	select(-kid, -kid_birth) %>%
	ungroup()

## Next we combine the normal set od individual years we want to calulate social
## indexes for with the dataset for the maternal SCI.
## We want to keep the regular dataset as including more individuals will provide
## better estimates of the universal slopes.
## (Ask me if you need a more detailed explanation)

iyol <- iyol_adults %>%  ## These are only adults
	mutate(temp = "adult") %>%
	bind_rows(maternal_iyol %>%  mutate(temp = "maternal"))  %>%
	ungroup()

## I have ran into issues with duplicated rows when e.g. a mom has birthday on same day as kid
## So lets find those rows and exclude them
## What is really essential is that the sname, grp, start and end are unique
iyol <- iyol %>%
	group_by(sname, grp, start, end) %>%
	arrange(sname, grp, start, end) %>%
	mutate(cases = n()) %>%
	filter(!(cases == 2 & temp == "maternal")) %>%
	select(-temp, -cases) %>%
	ungroup()

## Restrict to groups where the animal was present for at least 60 days
iyol_sub <- iyol %>%
	filter(days_present >= 60) %>%
	ungroup()

latest_version_date <- read_lines("./data/latest_version_date.txt")

## Lets save all this data
save(biograph_l,
		 members_l,
		 focals_l,
		 females_l,
		 grooming_l,
		 kid_data,
		 maternal_iyol,
		 iyol_adults,
		 iyol,
		 iyol_sub,
		 file = paste0("./data/dataset_for_social_index_calculations_including_maternal_", latest_version_date,
		 							".RData"))

load("./data/dataset_for_social_index_calculations_including_maternal_09JUN22.RData")



## The calculation of sci doesn't take very long, but DSI takes a long time.
## Its recomanded to calculate DSI on cluster

names(iyol_adults)
names(iyol)
# calculate-sci -----------------------------------------------------------

iyol_sub <- iyol_sub %>%  ungroup()

# Calculate grooming social connectedness index
sci_values <- sci(my_iyol = iyol_sub,
									members_l, focals_l, females_l, interactions_l = grooming_l,
									min_res_days = 1, parallel = TRUE, ncores = 4, directional = FALSE)

saveRDS(sci_values, paste0("./data/sci_including_maternal", latest_version_date, ".RDS"))

## The sci_values contains SCI values for all adult females and males as well
## as values that can be used for the maternal_sci calculation

maternal_sci_values <- kid_data %>%
	select(kid = sname, sname = mom, kid_birth = birth) %>%
	left_join(select(maternal_iyol, sname, grp, start, kid_birth = birth),
						by = c("sname", "kid_birth")) %>%
	left_join(sci_values, by = c("sname", "grp", "start")) %>%
	select(kid, kid_birth, sname, start,  grp, days_present, contains('SCI')) %>%
	pivot_longer(names_to = "index", values_to = "value", SCI_F_Dir:SCI_M) %>%
	group_by(kid, sname, kid_birth, start, index) %>%
	filter(!is.na(value) & days_present > 59) %>%
	## Currently there are values per time period per group
	## We need values per year, we get this by using a weighted mean
	summarise(value = weighted.mean(value,days_present, na.rm = TRUE),
						.groups = 'drop') %>%
	mutate(year = round(as.numeric((start - kid_birth)/365.25)) + 1) %>%
	group_by(kid, sname) %>%
	mutate(nr_years = n_distinct(year),
				 years = case_when(n_distinct(year) == 2 ~ 'both',
				 									n_distinct(year) == 1 & year == 1 ~ "1st only",
				 									n_distinct(year) == 1 & year == 2 ~ "2nd only",
				 									n_distinct(year) == 1 & is.na(year) ~ "no data")) %>%
	#filter(!(years %in% c('both', '1st only', 'no data'))) %>%
	group_by(kid, sname, years, index) %>%
	summarise(value = mean(value, na.rm = TRUE),
						.groups = 'drop') %>%
	pivot_wider(names_from = index, values_from = value)

maternal_sci_values <- write_csv(maternal_sci_values, paste0("./data/maternal_sci_values_", latest_version_date, ".csv"))



# calculate-dsi -----------------------------------------------------------

# Calculate population-level dyadic grooming index
# Warning: takes a really long time! Use the cluster
# I added some script and job scirpt to your folder on the cluster

dsi_pop <- dyadic_index(iyol_sub, biograph_l, members_l, focals_l, females_l,
												grooming_l, min_cores_days = 1, within_grp = FALSE,
												parallel = TRUE, directional = FALSE)

saveRDS(dsi_pop, paste0("data/dsi-pop_", Sys.Date(), ".RDS"))

# Summarize DSI variables for top partners in each year of life
dsi_pop_summary <- dyadic_index_summary(dsi_pop)

saveRDS(dsi_pop_summary, paste0("data/dsi-pop_summary_", Sys.Date(), ".RDS"))
