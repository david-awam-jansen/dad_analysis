# libraries
library(doSNOW)
library(foreach)
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

run_bootstrap_model <- function() {

	temp_data <- who_grooms_bootstrap %>%
		select(focal, focal_birth, focal_grp, mom, dad, AMales, start, end, nr_days
					 ## ages
					 , kid_age, mom_age, AMales_age
					 ## early adversity
					 , maternal_loss, maternal_rank, maternal_SCI_F, density, sibling, drought, cumulative_adversity

					 ## male details
					 , is_dad, mean_ordrank, mean_proprank
					 , estrous_presence, daily_d_days_rate, offspring_years, nr_pdads
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
		select(-previous_random, -next_random)

	print(ii)

################################################################################
################################################################################
## Full model dyadic strength
	temp_data_groom_strength <- temp_data %>%
		filter(does_groom == TRUE)

	full_model_groom_strength <- lmer(data = temp_data_groom_strength
																		, formula = zscore ~ 1
																		## ages
																		+ kid_age+ mom_age+ AMales_age
																		## early adversity
																		+ cumulative_adversity
																		## male details
																		+ is_dad + mean_ordrank
																		+ estrous_presence + daily_d_days_rate + offspring_years + nr_pdads
																		+ bootstrap_previous + bootstrap_next
																			+ (1|focal) +  (1|AMales)
																		, weights = log(nr_days)
																		, control=lmerControl(optimizer = "bobyqa",
																													optCtrl=list(maxfun=2e5)),
																		, na.action = 'na.fail')

	full_model_groom_strength_results <- full_model_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_groom_strength_results
							, file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_groom_strength_results_raw_v2.text'
							,  append = TRUE, col.names = FALSE)

	## full_strength_hybrid
	temp_data_groom_strength_hybrid <- temp_data %>%
		filter(does_groom == TRUE) %>%
		filter(!is.na(hybrid_score))

	full_model_hybrid_groom_strength <- update(object = full_model_groom_strength,
															 														. ~ . + hybrid_score,
															 														data = temp_data_groom_strength_hybrid)

	full_model_hybrid_groom_strength_results <- full_model_hybrid_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_hybrid_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_hybrid_groom_strength_results
							, file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_hybrid_groom_strength_results_raw.text'
							,  append = TRUE, col.names = FALSE)

	## full_strength_admix
	temp_data_groom_strength_admix <- temp_data %>%
		filter(does_groom == TRUE) %>%
		filter(!is.na(anubis_admix))

	full_model_admix_groom_strength <- update(object = full_model_groom_strength,
																						 . ~ . + anubis_admix,
																						 data = temp_data_groom_strength_admix)

	full_model_admix_groom_strength_results <- full_model_admix_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_admix_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_admix_groom_strength_results
							, file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_admix_groom_strength_results_raw.text'
							,  append = TRUE, col.names = FALSE)

	## Full model with maternal loss
	full_model_maternal_loss_groom_strength <- update(object = full_model_groom_strength,
															. ~ . - cumulative_adversity + maternal_loss,
															data = temp_data_groom_strength)


	full_model_maternal_loss_groom_strength_results <- full_model_maternal_loss_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_maternal_loss_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_maternal_loss_groom_strength_results,
							file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_maternal_loss_groom_strength_results_raw.text'
							,  append = TRUE, col.names = FALSE)

## no dad	strength
	temp_data_groom_strength_no_dad <- temp_data %>%
		filter(does_groom == TRUE) %>%
		filter(is_dad == FALSE)

	full_model_no_dad_groom_strength <- update(object = full_model_groom_strength,
																						. ~ . -is_dad,
																						data = temp_data_groom_strength_no_dad)

	full_model_no_dad_groom_strength_results <- full_model_no_dad_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_no_dad_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_no_dad_groom_strength_results
							, file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_no_dad_groom_strength_results_raw.text'
							,  append = TRUE, col.names = FALSE)

	## only dad	strength
	temp_data_groom_strength_only_dad <- temp_data %>%
		filter(does_groom == TRUE) %>%
		filter(is_dad == TRUE)

	full_model_only_dad_groom_strength <- update(object = full_model_groom_strength,
																							. ~ . -is_dad - estrous_presence,
																							data = temp_data_groom_strength_only_dad)

	full_model_only_dad_groom_strength_results <- full_model_only_dad_groom_strength %>%
		broom.mixed::tidy() %>%
		left_join(vif.mer(full_model_only_dad_groom_strength), by = 'term') %>%
		mutate(run = paste0("model_", ii))

	write.table(full_model_only_dad_groom_strength_results
							, file = '/afs/crc.nd.edu/group/ArchieLab/David/results/full_model_only_dad_groom_strength_results_raw.text'
							,  append = TRUE, col.names = FALSE)
}

nr_draws = 4

ncores <- 4

cl <- makeCluster(ncores,  outfile="")
registerDoSNOW(cl)

foreach(ii = 1:nr_draws,
				.packages = c('Matrix', 'lmerTest', 'tidyverse')) %dopar% {
					run_bootstrap_model()
				}
stopCluster(cl)
