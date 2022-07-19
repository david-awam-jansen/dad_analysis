source("./code/social_indexes/biographical-data.R")
latest_version_date <- read_lines("./data/latest_version_date.txt")

babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "jansen",
	dbname = "babase",
	password = "Bab00n3455")

# Get local copy of biograph table
biograph_l <- collect(tbl(babase, "biograph"))

# Make a members subset that excludes behavioral observation gaps
members_l <- subset_members(babase, .adults_only = FALSE)

# Subset other data sets used for sociality indices
focals_l <- subset_focals(babase, members_l)
females_l <- subset_females(members_l)

# Grooming
grooming_l <- subset_interactions(babase, members_l, my_acts = c("G"), .adults_only = FALSE)

# Make an individual-year-of-life data set for adults
iyol <- make_iyol(babase, members_l, focals_l, grooming_l, .adults_only = FALSE)

iyol_sub <- iyol %>%
	filter(days_present >= 60)

current_date <- toupper(format(today(), '%d%b%y'))
local_files <- ls(pattern = "_l")

save(list = c(local_files, "iyol", "iyol_sub"),
		 file = paste0("./data/data_set_for_social_index_analysis_", current_date, ".RData"))
