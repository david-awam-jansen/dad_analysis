library(broom)
library(tidyverse)
library(ramboseli)

# Open the 'traditional dsi' --------------------------------------------------
#  I ran into memmory issues when openimng all at onces
#  I will open and reopen files when needed
dsi <- readRDS('../ramboseli_juveniles_version2/own_code/dsi_pop_1_days_2019-06-10.RDS')
dsi <- dsi %>%
  mutate(nr_dyads = map_dbl(.x = subset, length)) %>%
  filter(nr_dyads > 0)


dsi_summary <- readRDS('../ramboseli_juveniles_version2/own_code/dsi_pop_summary_days_1_2019-06-10.RDS')
biograph_l <- readRDS('../ramboseli_juveniles_version2/data/biograph_l_2019-03-04.RDS')

moms <- biograph_l %>%
  select(sname, pid)  %>%
  mutate(mom = str_sub(pid,1,3))

# Functions -------------------------------------------------------------------
## explicit NA to deal with
explicit_na <- function(x) {
dims <- length(dim(x))
if (dims == 0L && length(x) == 0) {
  x <- ifelse(is.list(x) && !is.data.frame(x), list(NA_integer_), NA_integer_)
} else if (dims == 2L && nrow(x) == 0) {
  x[TRUE, ] <- NA_integer_
}
x
}

zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

## Adding focal data to the subset dataset in DSI
add_dsi_focal <- function(df) {
  df %>%
    select(focal_sex = sex, focal_sname = sname, focal_grp = grp,
           focal_start = start, subset) %>%
    group_by(focal_sex, focal_sname, focal_grp, focal_start) %>%
    unnest() %>%
    filter(!(focal_sex == "M" & dyad_type == 'F-F')) %>%
    mutate(focal = ((focal_sname == sname | focal_sname == partner) &
                      (focal_grp == grp | focal_grp == partner_grp))) %>%
    nest(-focal_sex, -focal_sname, -focal_grp, -focal_start, .key = subset) %>%
    rename(sname = focal_sname, grp = focal_grp, start = focal_start) %>%
    select(-focal_sex)
}

juv_subset <- function(df){
  df %>%
    select(focal_sex = sex, focal_sname = sname, focal_grp = grp,
           focal_start = start, subset) %>%
    group_by(focal_sex, focal_sname, focal_grp, focal_start) %>%
    unnest(subset) %>%
    unite("sname_key", sname, grp, remove = FALSE) %>%
    unite("partner_key", partner, partner_grp, remove = FALSE) %>%
    left_join(select(moms, sname = sname, sname_mom = mom), by="sname")  %>%
    left_join(select(moms, partner = sname, partner_mom= mom), by="partner") %>%
    mutate(maternal_grooms = if_else(sname == partner_mom | partner == sname_mom,
                                     true = "DSI_maternal",
                                     false = "DSI_maternal_excluded")) %>%
    nest(-focal_sex, -focal_sname, -focal_grp, -focal_start, .key = subset) %>%
    rename(sname = focal_sname, grp = focal_grp, start = focal_start) %>%
    select(-focal_sex)
}

#juv_subset <- purrr::safely(juv_subset, otherwise = tibble(no_focal_data = TRUE))

## running linear models for dsi
Slope_Model_LM_dsi <- function(df) {
  summary(lm(formula = log2_i_adj ~ log2OE, data = df))
}

## get the model results
get_models_dsi <- function(df) {
  df %>%
    ungroup() %>%
    filter(i_adj > 0) %>%
    select(dyad_type, contains("log2")) %>%
    group_by(dyad_type) %>%
    nest() %>%
    mutate(models = map(.x = data, .f = Slope_Model_LM_dsi)) %>%
    select(-data)
}

## filter on focals
focal_filter <- function(df) {
  df %>%
    filter(focal == TRUE)
}

## Doing the regression for DSI
get_dsi_regression <- function(df) {
  if(all(is.na(df$i_adj))) {
    NULL
  } else {
    df %>%
      filter(df$i_adj > 0) %>%
      lm(formula = log2_i_adj ~ log2OE)
  }
}

## Getting the residual values
get_dsi_residuals <- function (df) {
  if (all(df$i_adj == 0) | all(is.na(df$i_adj))) {
    return(tbl_df(NULL))
  }
  df <- df %>%
    mutate(log2_i_adj = if_else(i_adj > 0, log2_i_adj, NA_real_))

  df$res_value <-
    as.numeric(
      residuals(
        lm(data = df, log2_i_adj ~ log2OE, na.action=na.exclude)
      )
    )

  return(df)
}

universal_dsi_zscored <- function(df, focal_sex) {

  sv <- dsi_universal_slopes %>% filter(sex == focal_sex)
  names(sv) <- c("sex", "dyad_type", "B0", "B1")
  df <- df %>%inner_join(sv, by = "dyad_type")

  df <- df %>%
    filter(i_adj > 0) %>%
    group_by(dyad_type) %>%
    mutate(res_value =
             log2_i_adj - (B0 + (B1 * log2OE))) %>%
    ungroup()

  df <- df %>%
    group_by(dyad_type) %>%
    mutate_at(vars(starts_with('res')),.funs = list(zscore = ~zscore(.)))

  df <- df %>%
    filter(focal == TRUE) %>%
    select(dyad_type, res_i_adj, contains("zscore"))
}

## trick to deal with empty models (e.g. no DSI values)
possible_get_models_dsi <- possibly(get_models_dsi, otherwise = tibble(model_failed = TRUE))

safe_filter <- possibly(focal_filter, otherwise = tibble(no_focal_data = TRUE))

tidy2 <- purrr::possibly(tidy, "model_absent" == TRUE)

#universal_dsi_slope_values <- purrr::possibly(universal_dsi_slope_values, "model_absent" == TRUE)

# DSI manipulations ------------------------------------------------------------
dsi <- dsi %>%
  select(-subset) %>%
  left_join(add_dsi_focal(dsi), by = c("sname", "grp", "start"))

dsi <- dsi %>%
  select(-subset) %>%
  left_join(juv_subset(dsi), by = c("sname", "grp", "start"))

write_rds(dsi, "./data/dsi_with_additional_values.RDS.gz", compress = "gz")

dsi_only_focal <- dsi %>%
  mutate(focal_data = map(subset, safe_filter))  %>%
  select(-di)

write_rds(dsi_only_focal, "./data/dsi_only_focal.RDS.gz", compress = "gz")
rm(ls=dsi)

## Getting universal slopes values for DSI -------------------------------------
dsi_universal_slopes <- dsi_only_focal %>%
  unnest(focal_data) %>%
  group_by(sex, dyad_type) %>%
  nest(.key = focal_data)  %>%
  mutate(regression = purrr::map(focal_data,get_dsi_regression)) %>%
  mutate(model_results = regression %>%  map(.,.f = tidy2)) %>%
  unnest(model_results) %>%
  select(sex, dyad_type, term, estimate) %>%
  spread(key = term, value = estimate)  %>%
  set_names("sex", "dyad_type", "OE_dsi_intercept", "OE_dsi_slope")

## Running universal models for focals -----------------------------------------
dsi_focal_universal <- dsi_only_focal %>%
  unnest(focal_data) %>%
  group_by(sex, dyad_type) %>%
  nest(.key = focal_data)  %>%
  mutate(regression = purrr::map(focal_data,get_dsi_regression))  %>%
  mutate(res_value = purrr::map(focal_data,get_dsi_residuals)) %>%
  select(sex, dyad_type, res_value) %>%
  unnest(res_value) %>%
  group_by(sex, sex_class, sname, grp, start, dyad_type) %>%
  arrange(-res_value) %>%
  slice(1:3) %>%
  filter(res_i_adj > -9999) %>%
  summarise(DSI = mean(res_value)) %>%
  dplyr::mutate(DSI_type = case_when(
    sex_class == "AF" & dyad_type == "AF-AF" ~ "DSI_F_us",
    sex_class == "AF" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_M_us",
    sex_class == "AM" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_F_us",
    sex_class == "JUV" & dyad_type %in% c("AF-JUV", "JUV-AF") ~ "DSI_F_us",
    sex_class == "JUV" & dyad_type %in% c("AM-JUV", "JUV-AM" )~ "DSI_M_us",
    sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "DSI_J_us")) %>%
  select(-dyad_type) %>%
  spread(key = DSI_type, value = DSI)

rm(ls=dsi_only_focal)
dsi <- read_rds("./data/dsi_with_additional_values.RDS.gz")
dsi_summary <- readRDS('../ramboseli_juveniles_version2/own_code/dsi_pop_summary_days_1_2019-06-10.RDS')

## getting universal values that are zcored ------------------------------------
dsi_with_universal_zcores <- dsi %>%
  inner_join(select(dsi_summary, sname, grp, start, age_class),
             by = c("sname", "grp", "start", "age_class")) %>%
  mutate(zscored_dsi =
           pmap(.l = list(subset,sex), .f = universal_dsi_zscored)) %>%
  unnest(zscored_dsi) %>%
  mutate(res_value = res_value_zscore) %>%
  group_by(sex, sex_class, sname, grp, start, dyad_type) %>%
  arrange(-res_value) %>%
  slice(1:3) %>%
  filter(res_i_adj > -9999) %>%
  summarise(DSI = mean(res_value)) %>%
  dplyr::mutate(DSI_type = case_when(
    sex_class == "AF" & dyad_type == "AF-AF" ~ "DSI_F_zscore",
    sex_class == "AF" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_M_zscore",
    sex_class == "AM" & dyad_type %in% c("AF-AM", "AM-AF") ~ "DSI_F_zscore",
    sex_class == "JUV" & dyad_type %in% c("AF-JUV", "JUV-AF") ~ "DSI_F_zscore",
    sex_class == "JUV" & dyad_type %in% c("AM-JUV", "JUV-AM" )~ "DSI_M_zscore",
    sex_class == "JUV" & dyad_type == "JUV-JUV" ~ "DSI_J_zscore")) %>%
  select(-dyad_type) %>%
  spread(key = DSI_type, value = DSI)

dsi_summary <- dsi_summary %>%
  left_join(dsi_focal_universal) %>%
  left_join(dsi_with_universal_zcores)

write_rds(dsi_summary, "./data/dsi_summary_with_universal.RDS")

# Adding maternal grooms to data --------------------------------------------
## extra functions for maternal dsi

devtools::install_github("david-awam-jansen/ramboseli")

get_maternal_DSI <- function(df, sex_class, known_mom, focal_key) {
  if(sex_class == "JUV" & known_mom) {
    df %>%
      filter(!is.na(maternal_grooms))  %>%
      group_by(maternal_grooms) %>%
      nest() %>%
      mutate(regression = purrr::map(data, ramboseli::fit_dyadic_regression)) %>%
      unnest(regression) %>%
      filter(sname_key == focal_key | partner_key == focal_key) %>%
      group_by(maternal_grooms) %>%
      select(maternal_grooms, res_i_adj) %>%
      dplyr::filter(res_i_adj > -9999) %>%
      top_n(n = 3, res_i_adj)  %>%
      summarise(value = mean(res_i_adj, na.rm=TRUE)) %>%
      ungroup()  %>%
      spread(key = maternal_grooms, value = value) %>%
      explicit_na()
  } else {
    tibble(DSI_maternal = NA)
  }
}

## Check this part
get_maternal_DSI_universal <- function(df, sex_class, known_mom, focal_key) {
  if (sex_class == "JUV" & known_mom) {
    msv <- maternal_universal_slopes
    names(msv) <- c("maternal_grooms", "B0", "B1")
    df <- df %>% inner_join(msv, by = "maternal_grooms")

    df %>%
      filter(!is.na(maternal_grooms))  %>%
      group_by(maternal_grooms) %>%
      nest() %>%
      mutate(regression = purrr::map(data, ramboseli::fit_dyadic_regression)) %>%
      unnest(regression) %>%
      filter(sname_key == focal_key |
               partner_key == focal_key) %>%
      group_by(maternal_grooms) %>%
      select(maternal_grooms, res_i_adj) %>%
      dplyr::filter(res_i_adj > -9999) %>%
      top_n(n = 3, res_i_adj)  %>%
      summarise(value = mean(res_i_adj, na.rm = TRUE)) %>%
      ungroup()  %>%
      spread(key = maternal_grooms, value = value) %>%
      explicit_na()
  } else {
    tibble(DSI_maternal = NA)
  }
}

get_universal_maternal_zscores <- function(df, sex_class, known_mom) {

  if(sex_class == "JUV" & known_mom) {
    msv <- maternal_universal_slopes
    names(msv) <- c("maternal_grooms", "B0", "B1")
    df <- df %>%inner_join(msv, by = "maternal_grooms")

    df <- df %>%
      filter(i_adj > 0) %>%
      group_by(maternal_grooms) %>%
      mutate(res_value =
               log2_i_adj - (B0 + (B1 * log2OE))) %>%
      ungroup()

    df <- df %>%
      group_by(maternal_grooms) %>%
      mutate_at(vars(starts_with('res')),.funs = list(zscore = ~zscore(.)))

    df <- df %>%
      filter(focal == TRUE) %>%
      select(maternal_grooms, res_i_adj, contains("zscore"))
  } else {
    tibble(DSI_maternal_zscore = NA)
  }
}

dsi_summary <- read_rds("../R01_grant/data/dsi_summary_with_universal.RDS")

biograph_l <- readRDS('../ramboseli_juveniles_version2/data/biograph_l_2019-03-04.RDS')

moms <- biograph_l %>%
  select(sname, pid)  %>%
  mutate(mom = str_sub(pid,1,3))

dsi_summary <- dsi_summary %>%
  left_join(moms, by = "sname") %>%
  mutate(known_mom = !is.na(mom))

dsi <- read_rds("../R01_grant/data/dsi_with_additional_values.RDS.gz")

#get_maternal_DSI <- purrr::safely(get_maternal_DSI, otherwise = tibble(test = NA))

# Data manipulation
dsi_maternal <- dsi_summary   %>%
  filter(!is.na(DSI_F)) %>%
  select(sname, grp, start, sex, sex_class, known_mom) %>%
  left_join(select(dsi, sname, grp, start, sex_class, subset)) %>%
  unite("focal_key", sname, grp, remove = FALSE)  %>%
  #mutate(subset = map(.x = subset, .f = juv_subset))  %>%
  mutate(maternal_DSI_values =
           pmap(.l = list(subset, sex_class, known_mom, focal_key),
                .f = get_maternal_DSI)) %>%
  left_join(select(moms, sname = sname, mom = mom), by="sname")  %>%
  explicit_na() %>%
  unnest(maternal_DSI_values) %>%
  select(-subset) %>%
  explicit_na()


rm(ls=dsi_summary)
rm(ls=dsi)
rm(ls=dsi_focal)
dsi_only_focal <- read_rds("../R01_grant/data/dsi_only_focal.RDS.gz")


#  Getting universal values for maternal versions -----------------------------
maternal_universal_slopes <- dsi_only_focal %>%
  #mutate(subset = map(.x = subset, .f = juv_subset))  %>%
  filter(sex_class == "JUV") %>%
  unnest(subset) %>%
  filter(!is.na(maternal_grooms))  %>%
  group_by(maternal_grooms) %>%
  nest(.key = focal_data)  %>%
  mutate(regression = purrr::map(focal_data,get_dsi_regression)) %>%
  mutate(model_results = regression %>%  map(.,.f = tidy2)) %>%
  unnest(model_results)  %>%
  select(maternal_grooms, term, estimate)  %>%
  spread(key = term, value = estimate)  %>%
  set_names( "maternal_grooms", "OE_dsi_intercept", "OE_dsi_slope")

write.csv(maternal_universal_slopes, "../R01_grant/data/maternal_universal_slopes_values.cvs")
rm(ls=dsi_only_focal)


maternal_universal_slopes <- read.csv("../R01_grant/data/maternal_universal_slopes_values.cvs")
dsi_only_focal <- read_rds("../R01_grant/data/dsi_only_focal.RDS.gz")

## Running universal models for focals -----------------------------------------
dsi_focal_juvs <- dsi_only_focal %>%
  filter(sex_class == "JUV")

write_rds(dsi_focal_juvs, "../R01_grant/data/dsi_focal_juvs.RDS")

rm(ls=dsi_only_focal)

dsi_maternal_universal <- dsi_focal_juvs %>%
  filter(sex_class == "JUV") %>%
  unnest(subset) %>%
  filter(!is.na(maternal_grooms))  %>%
  group_by(maternal_grooms) %>%
  nest(.key = focal_data)   %>%
  mutate(regression = purrr::map(focal_data,get_dsi_regression))  %>%
  mutate(res_value = purrr::map(focal_data,get_dsi_residuals)) %>%
  unnest(res_value)

dsi_maternal_universal <- dsi_maternal_universal %>%
  group_by(sex, sex_class, sname, grp, start, maternal_grooms) %>%
  arrange(-res_value) %>%
  slice(1:3) %>%
  filter(res_i_adj > -9999) %>%
  summarise(DSI = mean(res_value)) %>%
  dplyr::mutate(DSI_type = case_when(
    maternal_grooms == "DSI_maternal" ~ "DSI_maternal_us",
    maternal_grooms == "DSI_maternal_excluded" ~ "DSI_maternal_excluded_us",
    maternal_grooms == "DSI_none" ~ NA_character_)) %>%
  select(-maternal_grooms)  %>%
  spread(key = DSI_type, value = DSI) %>%
  explicit_na()




dsi_summary <- read_rds("../R01_grant/data/dsi_summary_with_universal.RDS")
dsi_summary <- dsi_summary %>%
  left_join(moms, by = "sname") %>%
  mutate(known_mom = !is.na(mom))
dsi <- read_rds("../R01_grant/data/dsi_with_additional_values.RDS.gz")


## getting universal values that are zcored for maternal indexes
dsi_maternal_zcores <- dsi %>%
  inner_join(select(dsi_summary, sname, grp, start, age_class, known_mom),
             by = c("sname", "grp", "start", "age_class")) %>%
  mutate(zscored_materal =
           pmap(.l = list(subset,sex_class, known_mom), .f = get_universal_maternal_zscores)) %>%
  unnest(zscored_materal) %>%
  mutate(res_value = res_value_zscore) %>%
  mutate(maternal_grooms =
           if_else(is.na(maternal_grooms), "DSI_none", maternal_grooms)) %>%
  group_by(sex, sex_class, sname, grp, start, maternal_grooms) %>%
  arrange(-res_value) %>%
  slice(1:3) %>%
  filter(res_i_adj > -9999) %>%
  summarise(DSI = mean(res_value)) %>%
  dplyr::mutate(DSI_type = case_when(
    maternal_grooms == "DSI_maternal" ~ "DSI_maternal_zscore",
    maternal_grooms == "DSI_maternal_excluded" ~ "DSI_maternal_excluded_zscore",
    maternal_grooms == "DSI_none" ~ NA_character_)) %>%
  select(-maternal_grooms)  %>%
  spread(key = DSI_type, value = DSI) %>%
  explicit_na()


dsi_summary <- dsi_summary %>%
  left_join(dsi_maternal) %>%
  left_join(dsi_maternal_universal) %>%
  left_join(dsi_maternal_zcores)

write_rds(dsi_summary, "../R01_grant/data/dsi_summary_with_new_values_maternal.RDS")

