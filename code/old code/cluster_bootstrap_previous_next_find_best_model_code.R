# libraries
library(broom.mixed)
library(doSNOW)
library(foreach)
library(lmtest)
library(lmerTest)
library(parallel)
library(tidyverse)

setwd('/afs/crc.nd.edu/group/ArchieLab/David')

who_grooms_bootstrap <- read_csv("who_grooms_bootstrap.csv")

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

run_model <- function(model, variable) {
	print(paste0("looking at ", variable))
	update(model, paste0(' ~ . - ', variable))
}

get_aic <- function(model) {
	model %>%
		glance() %>%
		select(AIC) %>%
		pull()
}

get_lmtest <- function(full_model, new_model) {
	df <- lmtest::lrtest(full_model, new_model)
  df$`Pr(>Chisq)`[2]

  }

run_bootstrap_optimize <- function() {

	temp_data <- who_grooms_bootstrap %>%
		select(focal, focal_birth, focal_grp, mom, dad, AMales, start, end, nr_days
					 ## ages
					 , kid_age, mom_age, AMales_age
					 ## early adversity
					 , maternal_loss, maternal_rank, maternal_SCI_F, density, sibling, drought, cumulative_adversity
					 , observer_effort

					 ## male details
					 , is_dad, mean_ordrank, mean_proprank
					 , contains('estrous')
					 , daily_d_days_rate, offspring_years, nr_pdads
					 , previous_kid_with_mom, next_kid_with_mom ## see bootstrap

					 ## for subset
					 , hybrid_score, anubis_admix

					 ## responses
					 , does_groom,  zscore

					 ## create random
					 , previous_random := !!sym(paste0("previous_", ii))
					 , next_random := !!sym(paste0("next_", ii))) %>%
		mutate(bootstrap_previous = if_else(previous_kid_with_mom == .5, previous_random, previous_kid_with_mom),
					 bootstrap_next = if_else(next_kid_with_mom == .5, next_random, next_kid_with_mom)) %>%
		select(-previous_random, -next_random) %>%
		filter(!is.na(observer_effort))

	print(ii)

  ## Full model
  full_model <- glmer(data = temp_data %>%  sample_n(1000)
														, formula = does_groom ~ 1
											      + observer_effort
														## ages
														+ kid_age+ mom_age+ poly(AMales_age, 2)
														## early adversity
														+ cumulative_adversity
														## male details
														+ is_dad + mean_ordrank
														+ estrous_me +
														+	daily_d_days_rate + offspring_years + nr_pdads
														+ bootstrap_previous + bootstrap_next
														+ (1|focal) +  (1|AMales)
														, weights = log(nr_days)
														, family = binomial
														, control=glmerControl(optimizer = "bobyqa",
																									 optCtrl=list(maxfun=2e5))
														, na.action = 'na.fail')

  not_optimal <- TRUE
  model_to_improve <- full_model

  while(not_optimal == TRUE) {

  	last_model <- model_to_improve
  	print(paste0("There are ", length(names(fixef(last_model, add.dropped=FALSE))[-1]), " variables in the nmodel"))
  	model_checks <- tibble(variables = names(fixef(last_model, add.dropped=FALSE))[-1]) %>%
  		mutate(variables = case_when(str_detect(variables, "poly") == TRUE ~ "poly(AMales_age, 2)",
  																 str_detect(variables, "TRUE") ~ str_remove(variables, "TRUE"),
  																 TRUE ~ variables)) %>%
  		distinct() %>%
  		bind_cols(tibble(full_model = list(last_model))) %>%
  		tibble(new_model = pmap(.l = list(model = full_model, variable = variables),
  														.f = run_model))

  	to_be_removed <- model_checks %>%
  		mutate(new_AIC = map_dbl(new_model, get_aic),
  					 full_AIC = map_dbl(full_model, get_aic),
  					 dAIC = full_AIC  - new_AIC,
  					 lmtest = map2_dbl(full_model, new_model, .f = get_lmtest)) %>%
  		arrange(-lmtest) %>%
  		mutate(l_order = seq(1:n())) %>%
  		filter(lmtest > 0.001) %>%
  		slice(1)

  	if(length(to_be_removed$variables) == 1) {
  		not_optimal = TRUE
  		print(paste(to_be_removed$variables,  "is removed"))
  		model_to_improve <- to_be_removed$new_model[[1]]
  		} else {
  			not_optimal = FALSE
  			best_model <- last_model
  		}


}

  xx <- tibble(variables = names(fixef(full_model, add.dropped=FALSE))[-1],
  						 included = 1,
  						 model = "full") %>%
  	pivot_wider(names_from = "variables", values_from = "included") %>%
  	bind_rows(tibble(variables = names(fixef(best_model, add.dropped=FALSE))[-1],
  									 included = 1,
  									 model = as.character(ii)) %>%
  							pivot_wider(names_from = "variables", values_from = "included"))

  if(ii == 1) {
  	write_csv(xx, "./results/best_model_variables.csv", append = TRUE)
  } else {
  	write_csv(xx[2,], "./results/best_model_variables.csv", append = TRUE)
  }
}

nr_draws = 1000

ncores <- 32

cl <- makeCluster(ncores,  outfile="")
registerDoSNOW(cl)

foreach(ii = 1:nr_draws,
				.packages = c('Matrix', 'lmerTest', 'tidyverse')) %dopar% {
					run_bootstrap_optimize()
				}
stopCluster(cl)
