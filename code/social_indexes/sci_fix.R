library(broom)
library(tidyverse)

# Open the 'traditional dsi' --------------------------------------------------
#  I ran into memmory issues when openimng all at onces
#  I will open and reopen files when needed
#sci <- readRDS('../ramboseli_juveniles_version2/own_code/sci_60_days_2019-06-25.RDS')
sci <- readRDS("./data/David_sci_60_days_2019-07-22.RDS")

sci <-
  sci %>%
#  filter(sex_class == "AF") %>%
  inner_join(select(xdata, sname, age_class),
             by = c("sname", "age_class") )

library(cowplot)
sci %>%
  select(sex_class, SCI_F, SCI_F_Dir, SCI_F_Rec, SCI_F_tot) %>%
  pivot_longer(names_to = 'index', values_to = 'value',
               cols = SCI_F_Dir:SCI_F_tot) %>%
  ggplot(aes(x=SCI_F, y = value)) +
    geom_point() +
    facet_grid(index ~sex_class) +
    geom_abline(slope=1, intercept = 0, color = 'red') +
    cowplot::theme_half_open()

## Adding focal data to the subset dataset in DSI
add_sci_focal <- function(df) {
  df %>%
    select(focal_sex = sex, focal_sname = sname, focal_grp = grp,
           focal_start = start, subset) %>%
    group_by(focal_sex, focal_sname, focal_grp, focal_start) %>%
    unnest() %>%
    mutate(focal = ((focal_sname == sname) &
                      (focal_grp == grp))) %>%
    nest(-focal_sex, -focal_sname, -focal_grp, -focal_start, .key = subset) %>%
    rename(sname = focal_sname, grp = focal_grp, start = focal_start) %>%
    select(-focal_sex)
}

sci <- sci %>%
  select(-subset) %>%
  left_join(add_sci_focal(sci), by = c("sname", "grp", "start"))

get_sci_regression <- function(df) {
  df %>%
    lm(formula = value ~ log2OE)
  }

sci_universal_slopes <- sci %>%
  select(focal_sex = sex, focal_sname = sname, focal_grp = grp,
         focal_start = start, subset) %>%
  unnest() %>%
  filter(focal == TRUE) %>%
  select(sex_class, contains('log')) %>%
  gather(key = 'index', value = 'value',
               log2ItoF_daily:log2ItFme_daily) %>%
  filter(!is.na(value)) %>%
  group_by(sex_class, index) %>%
  nest() %>%
  mutate(regression = purrr::map(data,get_sci_regression)) %>%
  mutate(model_results = regression %>%  map(.,.f = tidy)) %>%
  unnest(model_results) %>%
  select(sex_class, index, term, estimate) %>%
  spread(key = term, value = estimate)  %>%
  set_names("sex_class", "index", "B0", "B1")

sci_universal_slopes

get_sci_residuals <- function (df) {
  if(sum(!is.na(df$value)==0)) {
    df$res_value <- NA_real_

  } else {
    df$res_value <- as.numeric(
      residuals(
        lm(data = df, value ~ log2OE, na.action=na.exclude)
        )
      )
  }
    return(df)
}

zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

sci_focal_universal <- sci %>%
  filter(!is.na(SCI_F)) %>%
  select(start, subset) %>%
  unnest() %>%
  filter(focal == TRUE) %>%
  select(sex_class, sname, grp, start, contains('log2')) %>%
  gather(key = 'index', value = 'value',
               log2ItoF_daily:log2ItFme_daily) %>%
  # pivot_longer(names_to = 'index', values_to = 'value',
  #              log2ItoF_daily:log2ItFme_daily) %>%  filter(!is.na(value)) %>%
  group_by(sex_class, index) %>%
  nest() %>%
  mutate(res_value = purrr::map(data,get_sci_residuals)) %>%
  unnest(res_value) %>%
  select(sname, grp, start, sex_class, index, res_value)  %>%
  mutate(index = stringr::str_replace(index, pattern = "log2", replacement = "res")) %>%
  #pivot_wider(names_from = index, values_from = res_value) %>%
  spread(key = index, value = res_value) %>%
  mutate(SCI_F_us = (resItoF_daily + resIfromF_daily)/2,
         SCI_F_tot_us = resItF_daily,
         SCI_Fme_us = (resItoFme_daily + resIfromFme_daily)/2,
         SCI_Fme_tot_us = (resItFme_daily)/2,
         SCI_M_us = (resItoM_daily + resIfromM_daily)/2,
         SCI_M_tot_us = resItM_daily,
         SCI_J_us = (resItoJ_daily + resIfromJ_daily)/2,
         SCI_J_tot_us = resItJ_daily
         ) %>%
  select(sex_class, sname, grp, start, contains('SCI_'))

sci_with_universal_zcores <- sci %>%
    select(focal_sex = sex, focal_sname = sname, focal_grp = grp,
           focal_start = start, subset) %>%
    unnest() %>%
    select(focal_sex, focal_sname, focal_grp, focal_start,
           sex_class, sname, grp, start, focal, contains("log2")) %>%
    gather(key = 'index', value = 'value', log2ItoF_daily:log2ItFme_daily) %>%
    #pivot_longer(names_to = 'index', values_to = 'value', log2ItoF_daily:log2ItFme_daily) %>%
    filter(!is.na(value)) %>%
    inner_join(sci_universal_slopes, by=c('sex_class', 'index')) %>%
    mutate(res_value = value - (B0 + (B1 * log2OE))) %>%
    group_by(focal_sex, focal_sname, focal_grp, focal_start, index) %>%
    select(focal_sex, focal_sname, focal_grp, focal_start,
           sex_class, sname, grp, start, focal, index, contains('res'))  %>%
    mutate_at(vars(starts_with('res')),.funs = list(zscore = ~zscore(.))) %>%
    ungroup() %>%
    filter(focal == TRUE) %>%
    mutate(index = stringr::str_replace(index, pattern = "log2", replacement = "res")) %>%
    select(-res_value) %>%
    #pivot_wider(names_from = index, values_from = zscore) %>%
    spread(key = index, value =  zscore) %>%
    mutate(SCI_F_zscore = (resItoF_daily + resIfromF_daily)/2,
           SCI_F_tot_zscore = resItF_daily,
           SCI_Fme_zscore = (resItoFme_daily + resIfromFme_daily)/2,
           SCI_Fme_tot_zscore = resItFme_daily,
           SCI_M_zscore = (resItoM_daily + resIfromM_daily)/2,
           SCI_M_tot_zscore = resItM_daily,
           SCI_J_zscore = (resItoJ_daily + resIfromJ_daily)/2,
           SCI_J_tot_zscore = resItJ_daily) %>%
    select(sex_class, sname, grp, start=focal_start, contains('SCI_'))

sci2 <- sci %>%
  left_join(sci_focal_universal) %>%
  left_join(sci_with_universal_zcores)

write_rds(sci2, "./data/sci_with_new_values_zscores.RDS")
write_rds(sci2, "./data/sci_with_new_values_zscores_22JUL19.RDS")
