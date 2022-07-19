library(tidyverse)

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_dad_analysis_", latest_version_date, ".RData"))

xdata_l <- biograph_l %>%
  filter(birth > '1971-07-01'     ## exclude animals born before the onset of monitoring
         & birth > '1976-08-15'   ## exclude animals born before rainfall estimates
         & matgrp < 3             ## exclude lodge group
         & matgrp != 1.300        ## exclude proton's group because the data are so thin
         & bstatus == 0           ## exclude animals with bad age estimates
         & statdate > birth       ## exclude animals who lived for < 1 day
         & sex != 'U'             ## exclude individuals that died before they could be sexed
  ) %>%
  filter(birth + years(4) < '2020-05-31') %>%
  inner_join(select(parents_l, sname = kid, mom, dad)) %>%
  mutate(birth_rnkdate = floor_date(birth, 'month'))

## i. Density
## Number of adults in group on date of birth
## Only use if group was a was a study group (This is stricter the Tung et al analysis)
## Let for instance look at an example of filtering out members data based on
## groups history and behave gaps
members_l <- members_l %>%
  inner_join(select(biograph_l, sname, sex)) %>%
  left_join(select(maturedates_l, sname, matured))

adult_group_sizes_l <- members_l %>%
  filter(date > matured) %>%
  select(grp, date, sname) %>%
  inner_join(xdata_l %>%
               select(grp = matgrp, date = birth) %>%
               distinct(),
             by = c("grp", "date")) %>%
  group_by(grp, date) %>%
  summarise(nr_adults = n()) %>%
  arrange(desc(date))

xdata_l <- xdata_l %>%
  left_join(select(adult_group_sizes_l, matgrp = grp, birth = date, density = nr_adults))

## ii. maternal rank
maternal_ranks_l <- ranks_l %>%
  filter(rnktype == 'ADF') %>%
  select(sname, grp, rnkdate, rank) %>%
  inner_join(select(xdata_l, sname = mom, kid = sname, grp = matgrp, rnkdate = birth_rnkdate)) %>%
  select(mom = sname, sname = kid, maternal_rank = rank)

xdata_l <- xdata_l %>%
    left_join(select(maternal_ranks_l, mom, sname, maternal_rank))
xdata_l %>% filter(is.na(maternal_rank))

## iii. Rainfall and drought
get_rainfall <- function(focal_birth) {
  rainfall_l %>%
    filter(date >= focal_birth,
           date <= focal_birth + years(1)) %>%
    select(rain) %>%
    pull() %>%
    sum(na.rm = TRUE)
}

xdata_l <- xdata_l %>%
  mutate(rainfall = map_dbl(.x = birth, .f = get_rainfall)) %>%
  mutate(drought = rainfall <= 200)

## iv. sibling
get_siblings <- function(focal_birth, focal_mom) {
  xx <- pregdata_l %>%
    filter(mom == focal_mom) %>%
    filter(birth > ymd(focal_birth),
           birth <= ymd(focal_birth) + months(18)) %>%
    select(sibling = kid, sibling_birth = birth, sibling_dad = dad)

    if(nrow(xx) == 0) {
      tibble(sibling = NA, sibling_birth = NA, sibling_dad = NA, has_sibling = 0)
    } else(
      xx %>%
      mutate(has_sibling = n())
    )

}

xdata_l <- xdata_l %>%
	 mutate(sibling = map2(.x = birth, .y = mom, .f = get_siblings)) %>%
  unnest(cols = c(sibling), keep_empty = TRUE)

xdata_l <- xdata_l %>%
  mutate(age_at_birth_sibling = (as.numeric((ymd(sibling_birth) - ymd(birth))/365.25)))

## v. maternal loss
xdata_l <- xdata_l %>%
  inner_join(maternal_loss_l)

## vi. maternal SCI to be added
maternal_sci_values <- read_csv(paste0("./data/maternal_sci_values_31MAY22.csv"))

xdata_l <- xdata_l %>%
	left_join(select(maternal_sci_values, mom = sname, sname = kid, years, contains("SCI")))

xdata_l <- xdata_l %>%
	mutate(sex= factor(sex, labels = c("Females", "Males")),
				 maternal_loss = factor(maternal_loss,
				 											 levels = c(FALSE, TRUE),
				 											 labels = c("Mom alive", "Mom dies")),
				 maternal_SCI_binary = factor(maternal_SCI_binary,
				 														 levels = c(FALSE, TRUE),
				 														 labels = c("Socially connected", "Socially isolated")),
				 maternal_rank_binary = factor(maternal_rank_binary,
				 															levels = c(FALSE, TRUE),
				 															labels = c("Normal rank", "Low rank")),
				 sibling = factor(sibling,
				 								 levels = c(FALSE, TRUE),
				 								 labels = c("Absent", "Com sib. present")),
				 drought = factor(drought,
				 								 levels = c(FALSE, TRUE),
				 								 labels = c("Wet", "Dry")),
				 density_binary = factor(density_binary,
				 												levels = c(FALSE, TRUE),
				 												labels = c("Normal group", "Large group")),
				 first_born = parity == 1,
				 first_born = factor(first_born,
				 										levels = c(FALSE, TRUE),
				 										labels = c("Experienced mother", "First born")),
				 habitat_quality = case_when(
				 	matgrp == 1 & year(ymd(birth)) <= 1987 ~  "low",
				 	matgrp == 2 & year(ymd(birth)) <= 1991 ~ "low",
				 	TRUE ~ "high"),
				 habitat_quality = factor(habitat_quality))

write.csv(x = xdata_l, file = paste("./data/", "ea_dataset_less_restricted_",latest_version_date, ".csv", sep=""))


# custom_order <- c("SCI_F", "SCI_M", "SCI_J", "DSI_F", "DSI_F_mom_excluded", "DSI_M", "SUB_J")
desc_missing <- NULL

numbers <- xdata_l %>%
	dplyr::group_by(sex) %>%
	summarise(number = n())

xdata_l %>%
	select(sex,dad, density, maternal_rank, maternal_loss, SCI_F, has_sibling, drought) %>%
	dplyr::group_by(sex) %>%
	naniar::miss_var_summary() %>%
	#mutate((n_miss/pct_miss) * 100)
	dplyr::mutate(variable = factor(variable,
																	levels = (sort(unique(variable),
																								 decreasing = TRUE))))  %>%
	# forcats::fct_relevel(rev(desc_missing), after = Inf) %>%
	# forcats::fct_relevel(rev(custom_order), after = Inf)) %>%
	left_join(numbers, by = c("sex")) %>%
	mutate(number = paste(n_miss, " out of ", number, sep = "\n")) %>%
	ggplot(aes(sex, variable, fill = pct_miss)) +
	geom_tile(color = "black", size = 2)  +
	geom_text(aes(label = number), color = "black", size = 3)  +
	facet_wrap(~sex) +
	scale_fill_gradient2(low = "white", high = "firebrick",
											 limits = c(0, 100), name = "% Miss") +
	theme_minimal() + theme(legend.position = "bottom") +
	ggtitle("Missing data overview for included juveniles")



##After adding maternal_SCI there are still 7 moms for which maternal rank is unknow.
missing_ranks <- xdata %>%  filter(is.na(maternal_rank) & !is.na(maternal_SCI_F)) %>%
  select(mom = mom.x, birth, statdate, sex, matgrp)

missing_ranks %>%
  mutate(statage = as.numeric((statdate - birth)/365.25)) %>%
  mutate(age_check = statage > 4) %>%
  janitor::tabyl(sex, age_check)

# ## get all rannk data for those females
# extra_ranks <-data.table(dbGetQuery(paste("SELECT *
#                                           FROM ranks r
#                                           JOIN groups_history gh ON gh.gid = r.grp
#                                           WHERE rnkdate > gh.permanent  AND rnkdate < LEAST(gh.impermanent, gh.cease_to_exist, gh.last_reg_census, '2020-01-01') AND
#                                           rnktype = 'ADF' AND
#                                                 sname IN ('",paste(missing_ranks$mom, collapse = "', '") ,"');", sep=""),
#                                     con = con))
#
# extra_ranks[, rnkdate := ymd(rnkdate) ]
#
#
# missing_ranks <- missing_ranks %>%
#   mutate(sname = mom)
#
# missing_ranks <- missing_ranks %>%
#   left_join(missing_ranks %>%
#                inner_join(extra_ranks, by = "sname") %>%
#                filter(rnkdate > birth) %>%
#                group_by(sname, grp) %>%
#                summarise(nearest_rnkdate_ahead = min(rnkdate))) %>%
# left_join(missing_ranks %>%
#   left_join(missing_ranks %>%
#               inner_join(extra_ranks, by = "sname") %>%
#               filter(rnkdate < birth) %>%
#               group_by(sname, grp) %>%
#               summarise(nearest_rnkdate_before = max(rnkdate))),
#   by = c("sname", "birth", "mom", "matgrp", "grp")) %>%
#   mutate(birth_rnkdate = floor_date(birth, unit = "month")) %>%
#   mutate(rnkdate = case_when(
#     abs(birth_rnkdate - nearest_rnkdate_before) > abs(birth_rnkdate - nearest_rnkdate_ahead) ~
#       nearest_rnkdate_ahead,
#     abs(birth_rnkdate - nearest_rnkdate_before) < abs(birth_rnkdate - nearest_rnkdate_ahead) ~
#       nearest_rnkdate_before,
#     TRUE ~ nearest_rnkdate_ahead)) %>%
#   left_join(select(extra_ranks, mom = sname, grp, rnkdate, rank),
#             by = c("mom", "grp", "rnkdate"))
#
# missing_ranks %>%
#   mutate(test = ceiling(as.numeric(lubridate::as.difftime(rnkdate - birth)/30)))
#
# ## merge this with the data
# xdata <- xdata %>%
#   left_join(missing_ranks %>%
#               select(mom, birth, nearest_rankdate = rnkdate, nearest_rank = rank)
#   )
#
# xdata <- xdata %>%
#   mutate(maternal_rank =
#            case_when(is.na(maternal_rank) ~ nearest_rank,
#                      TRUE ~ maternal_rank))
# xdata <- xdata %>%
#   select(-nearest_rank, -nearest_rankdate)


