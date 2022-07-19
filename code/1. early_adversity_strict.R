library(tidyverse)
library(RPostgreSQL)
## setting up data connection
dbuser <- readline(prompt="Enter babase username: ")
dbpass <- readline(prompt="Enter babase password: ")
dbname <- "babase";
dbhost <- "localhost";
dbport <- "22222"
drv <- dbDriver("PostgreSQL")

##open a connection to the database
babase <- dbConnect(drv, host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)

selected_behave_gaps <- read_csv('./data/selected_behave_gaps.csv')

emp_snames <- tbl(con, dbplyr::in_schema("genetic_sequencing", "emp_tissue")) %>%
  filter(!is.na(sname)) %>%
  select(sname) %>%
  distinct() %>%
  pull()


biographic_data <- tbl(con, "biograph")

xdata <- biographic_data %>%
  filter(sname %in% emp_snames) %>%
  mutate(ea_check = (birth > '1971-07-01'     ## exclude animals born before the onset of monitoring
         & birth > '1976-08-15'    ## exclude animals born before rainfall estimates
         & matgrp < 3                ## exclude lodge group
         & matgrp != 1.300           ## exclude proton's group because the data are so thin
         & bstatus == 0             ## exclude animals with bad age estimates
         & statdate > birth               ## exclude animals who lived for < 1 day
         & !is.na(pid)                   ## exclude individuals without known moms
         & sex != 'U'
         #& birth <= '2019-12-01' - years(4)
         # & birth <= '2016-01-01'
         )
  ) %>%
  mutate(mom := str_sub(pid, 1, 3)) %>%
  mutate(rnkdate = rnkdate(birth)) %>%
  select(ea_check, sname, sex, birth, rnkdate, matgrp, pid, mom,
         bstatus, status, dcause, statdate)

xdata_l <- xdata %>%  collect()

## members data
members_filtered <- tbl(con, 'members') %>%
  inner_join(select(xdata, grp = matgrp, date = birth) %>%  distinct()) %>%
  inner_join(select(tbl(con, "groups_history"), grp = gid,
                    permanent, impermanent, cease_to_exist, last_reg_census, study_grp)) %>%
  filter(!is.na(study_grp) & grp != 9) %>%
  group_by(grp) %>%
  mutate(last_date = pmin(impermanent, cease_to_exist,
                          last_reg_census, "2020-01-01")) %>%
  ungroup() %>%
  select(-impermanent, cease_to_exist, last_reg_census) %>%
  filter(date >= permanent & date <= last_date) %>%
  left_join(tbl(con, "maturedates"), by = "sname") %>%
  filter(date > matured) %>%
  select(membid, sname, grp, date)

# Find behavior gaps that overlap records in members_keep
excluded_density <- selected_behave_gaps %>%
  filter(exclude_density) %>%
  select(bgid) %>%
  pull()

bg <- members_filtered %>%
  dplyr::left_join(tbl(con, "behave_gaps"), by = "grp") %>%
  dplyr::filter(date >= gap_start & date <= gap_end) %>%
  select(membid, bgid, gap_start, gap_end, notes) %>%
  filter(bgid %in% excluded_density)

members_filtered <- members_filtered %>%
  anti_join(bg, by = 'membid')

## i. Density
density_raw <- members_filtered %>%
    group_by(grp, date) %>%
    summarise(density = n()) %>%
  collect()

xdata_l <- xdata_l %>%
  left_join(select(density_raw, matgrp = grp, birth = date, density))

# xdata_l %>%
#   mutate(has_density = is.na(density)) %>%
#   janitor::tabyl(has_density)

## 2. maternal rank
excluded_rank <- selected_behave_gaps %>%
  filter(exclude_rank) %>%
  select(bgid) %>%
  pull()

maternal_rank_raw <- tbl(con, "ranks") %>%
  filter(rnktype == 'ADF') %>%
  inner_join(select(xdata, grp = matgrp, rnkdate, birth, sname = mom)) %>%
  distinct() %>%
  inner_join(select(tbl(con, "groups_history"), grp = gid,
                    permanent, impermanent, cease_to_exist, last_reg_census)) %>%
  group_by(grp) %>%
  mutate(last_date = pmin(impermanent, cease_to_exist,
                          last_reg_census, "2020-01-01")) %>%
  ungroup() %>%
  filter(birth >= permanent & birth <= last_date) %>%
  select(rnkid, mom = sname, birth, rnkdate, grp, rnktype, rank)

rank_bg <- maternal_rank_raw %>%
  dplyr::left_join(tbl(con, "behave_gaps"), by = "grp") %>%
  dplyr::filter(birth >= gap_start & birth <= gap_end) %>%
  select(rnkid, bgid) %>%
  filter(bgid %in% excluded_rank)

maternal_rank_raw <- maternal_rank_raw %>%
  anti_join(rank_bg, by = 'rnkid') %>%
  collect()

xdata_l <- xdata_l %>%
  left_join(maternal_rank_raw %>%
               select(mom, birth, maternal_rank = rank),
            by = c("mom", "birth"))

xdata_l %>%
  mutate(has_rank = is.na(maternal_rank)) %>%
  janitor::tabyl(sex, has_rank) %>%
  janitor::adorn_totals(where = c("row", "col"))


## iii. Rainfall and drought
rainfall <- tbl(con, "wreadings") %>%
  left_join(tbl(con, "raingauges")) %>%
  select(wrdaytime, rain) %>%
  collect() %>%
  #mutate(wrdaytime = as.character(wrdaytime)) %>%
  separate(wrdaytime, into = c("date", NA), sep = ' ')

library(lubridate)
get_rainfall <- function(focal, start) {
  start <- ymd(start)
  x <- rainfall %>%
    filter(date >= start &
            date <= start + years(1) - days(1)) %>%
    select(rain) %>%
    pull() %>%
    sum(na.rm = TRUE)

  print(paste(focal, start))

  return(x)
}

xdata_l <- xdata_l %>%
  group_by(birth) %>%
  mutate(rainfall = map2_dbl(.x = sname, .y = birth, .f = get_rainfall),
         drought = rainfall <= 200)

xdata_l %>%
  mutate(has_rainfall = !is.na(rainfall)) %>%
  janitor::tabyl(sex, drought) %>%
  janitor::adorn_totals(where = c("row", "col"))

## v. sibling
sibling_data <- tbl(con, 'parents') %>%
  select(sname = kid, mom) %>%
  inner_join(select(tbl(con, 'biograph'), sname, birth)) %>%
  left_join(select(tbl(con, 'parents'), sibling_sname = kid, mom), by = "mom") %>%
  filter(sname != sibling_sname) %>%
  inner_join(select(tbl(con, 'biograph'), sibling_sname = sname, sibling_birth = birth, statdate)) %>%
  collect() %>%
  filter(sibling_birth > birth &
         sibling_birth < birth + months(18))  %>%
  group_by(sname) %>%
  #filter(n() >1) %>%
  summarise(sibling_sname = paste0(sibling_sname, collapse = " "),
            sibling_birth = min(sibling_birth))


xdata_l <- xdata_l %>%
  left_join(sibling_data, by = c('sname')) %>%
  mutate(sibling = !is.na(sibling_sname)) %>%
  mutate(age_at_birth_sibling = as.numeric((ymd(sibling_birth) - ymd(birth))/365.25))

## vi. maternal loss
maternal_loss_data <- tbl(con, 'biograph') %>%
  select(sname, birth) %>%
  left_join(select(tbl(con, 'parents'), sname = kid, mom)) %>%
  left_join(select(tbl(con, 'biograph'), mom = sname, statdate)) %>%
  mutate(maternal_loss_age = as.numeric((statdate - birth)/365.25),
         maternal_loss = maternal_loss_age < 4) %>%
  select(sname, contains("maternal")) %>%
  collect()

xdata_l <- xdata_l %>%
  inner_join(maternal_loss_data)

## iV. maternal SCI to be added
xdata_ea_maternal_SCI_focals <- xdata_l %>%  select(sname, mom, birth)
## uncomment if new set is needed
# write.csv(x = xdata_ea_maternal_SCI_focals, file = paste("./data/", "xdata_for_maternal_SCI_", today(), ".csv", sep=""))

#maternal_iyol <- xdata_ea_maternal_SCI_focals %>%
# xdata_ea_maternal_SCI_focals %>%
#   select(kid = sname, sname = mom, kid_birth = birth, birth) %>%
#   mutate(birth = kid_birth) %>%
#   group_by(sname, kid) %>%
#   mutate(start = list(seq(from = kid_birth, to = kid_birth + years(2) - days(1), by = 'years'))) %>%   unnest(cols = c(start)) %>%
#   ungroup() %>%
#   mutate(end = start + years(1) - days(1)) %>%
#   left_join(select(tbl(babase, 'maturedates'), sname, matured) %>%  collect(), by = 'sname') %>%
#   left_join(select(tbl(babase, 'rankdates'), sname, ranked) %>%  collect(), by = 'sname') %>%
#   left_join(select(tbl(babase, 'biograph'), sname, statdate, sex) %>%  collect(), by = 'sname') %>%
#   mutate(first_start_date = if_else(sex == 'F', matured, ranked)) %>%
#   mutate(midpoint = start + floor((end - start)/2),
#          age_start_yrs = as.numeric(start - birth)/365.25,
#          age_class = floor(plyr::round_any(age_start_yrs, 0.005)) + 1,
#          sex_class = case_when(start < matured | is.na(matured) ~ "JUV",
#                                sex == "F" & start >=  matured ~ "AF",
#                                sex == "M" & start >= ranked ~  "AM",
#                                sex == "M" & start >= matured &
#                                  (start < ranked | is.na(ranked)) ~ "SM")

# sci <- read_rds('./data/sci_with_universal_values_20191204.RDS')
# maternal_iyol <- read_rds("maternal_iyol.RDS")

maternal_sci_values <- xdata_ea_maternal_SCI_focals %>%
  select(kid = sname, sname = mom, kid_birth = birth) %>%
  left_join(select(maternal_iyol, sname, grp, start, kid_birth = birth)) %>%
  left_join(sci) %>%
  select(kid, kid_birth, sname, start,  grp, days_present, SCI_F,SCI_M) %>%
  pivot_longer(names_to = "index", values_to = "value", SCI_F:SCI_M) %>%
  group_by(kid, sname, kid_birth, start, index) %>%
  filter(!is.na(value) & days_present > 59) %>%
  summarise(value = weighted.mean(value,days_present, na.rm = TRUE)) %>%
  mutate(year = round(as.numeric((start - kid_birth)/365.25)) + 1) %>%
  group_by(kid, sname) %>%
  mutate(nr_years = n_distinct(year),
         years = case_when(n_distinct(year) == 2 ~ 'both',
                           n_distinct(year) == 1 & year == 1 ~ "1st only",
                           n_distinct(year) == 1 & year == 2 ~ "2nd only",
                           n_distinct(year) == 1 & is.na(year) ~ "no data")) %>%
  #filter(!(years %in% c('both', '1st only', 'no data'))) %>%
  group_by(kid, sname, years, index) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  pivot_wider(names_from = index, values_from = value)

xdata_l <- xdata_l %>%
  left_join(maternal_sci_values %>%
              filter(years %in% c("both", "1st only")) %>%
              ungroup() %>%
              select(mom = sname, sname = kid,
                     years,
                     maternal_SCI_F = SCI_F),
            by = c("mom", "sname"))

mauna_data <- read_rds("./data/metadataForDavid_DSIF_ea.rds")
mauna_snames <- mauna_data %>%
  select(sname) %>%
  distinct() %>%
  pull()

ea_data_for_mauna <- xdata_l %>%
  filter(sname %in% mauna_snames) %>%
  select(sname, sex, birth, matgrp, mom, bstatus, status, statdate,
         density, maternal_loss, maternal_rank, maternal_SCI_F, sibling,
         group_size = density, drought,
         ea_check) %>%
  group_by(sname) %>%
  mutate(included_cases = complete.cases(
    maternal_loss, maternal_SCI_F,
    maternal_rank, sibling, group_size, drought))

write.csv(x = ea_data_for_mauna,
          file = paste("./data/", "ea_dataset_for_mauna", ".csv", sep=""))

mauna_iyol <- mauna_data %>%
  filter(sex == "F") %>%
  select(sname, collection_date) %>%
  distinct() %>%
  as_tibble() %>%
  mutate(start = collection_date - years(1) + days(1),
         end = collection_date) %>%
  left_join(select(tbl(babase, 'maturedates'), sname, matured) %>%  collect(), by = 'sname') %>%
  left_join(select(tbl(babase, 'rankdates'), sname, ranked) %>%  collect(), by = 'sname') %>%
  left_join(select(tbl(babase, 'biograph'), sname, birth, statdate, sex) %>%  collect(), by = 'sname') %>%
  mutate(first_start_date = if_else(sex == 'F', matured, ranked)) %>%
  mutate(midpoint = start + floor((end - start)/2),
         age_start_yrs = as.numeric(start - birth)/365.25,
         age_class = floor(plyr::round_any(age_start_yrs, 0.005)) + 1)

mauna_iyol

mauna_iyol <-mauna_iyol %>%
  inner_join(mauna_iyol %>%
               left_join(mauna_iyol %>%
                           left_join(select(members_l, sname, date, grp), by = c("sname")) %>%
                           filter(date >= start & date <= end) %>%
                           distinct(sname, start, grp)) %>%
               inner_join(dplyr::select(members_l, sname, grp, date), by = c("sname", "grp")) %>%
               filter(date >= start & date <= end) %>%
               group_by(sname, grp, start, end) %>%
               summarise(days_present = n()) %>%
               arrange(sname, grp, start, end),
             by = c("sname", "start", "end")) %>%
  mutate(birth_dates = start) %>%
  select(sname, grp, start, end, days_present, sex, birth, first_start_date,
         ranked, matured, statdate, midpoint, birth_dates, midpoint,
         age_start_yrs,  age_class)

mauna_iyol_sub <- mauna_iyol %>%
  filter(days_present > 60) %>%
  select(-matured, - ranked)

write_rds(mauna_iyol_sub, "./data/mauna_iyol_sub.rds")

list.of.packages <- list("foreach", "doSNOW", "parallel", "tidyverse",
                         "lubridate", "dbplyr", "purrrlyr", "RPostgreSQL",
                         "zoo", "ramboseli")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages, require, character.only = T)


load("./data/social_index_dataset_20200312.RData")
mauna_iyol_sub <- read_rds("./data/mauna_iyol_sub.rds")

dsi_pop <- dyadic_index(mauna_iyol_sub,
                        biograph_l, members_l, focals_l, females_l,
                        grooming_l, min_cores_days = 1, within_grp = FALSE,
                        parallel = TRUE, directional = FALSE, ncores = 7)


save(dsi_pop, file = "dsi_pop_for_Mauna.Rdata")



