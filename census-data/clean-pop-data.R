
# Load the population data ---------------

pop <- readxl::read_excel(here::here("census-data/co-est2023-pop.xlsx"),
                          skip = 4,
                          col_names = F)

names(pop) <- c("area",
                "base_april_2020",
                "est_july_2020",
                "est_july_2021",
                "est_july_2022",
                "est_july_2023")


# split the string by comma ----------
pop$county_clean <- stringr::str_replace_all(pop$area, "\\.", "")

pop <- pop[-1,]
pop <- pop[-c(3145:3150),]

spt <- stringr::str_split(pop$county_clean,
                          "\\,")

counties <- do.call(rbind,(lapply(spt,function(x){
 county <- x[[1]]
 state <- x[[2]]
 out <- data.frame(county = county,
                   state = state)
 out
})))

pop$county_name <- counties$county
pop$state_name <- counties$state
pop$state_name <- stringr::str_replace_all(pop$state_name,
                                           "^ ","")

# map on abbreviations ------

sts <- read.csv(here::here("usa_adj_mat_files/states.csv"))

pop$abbrev <-  plyr::mapvalues(pop$state_name,
                               from = sts$State,
                               to = sts$Abbreviation)


pop$full_county_name <- paste(pop$county_name,
                              pop$abbrev)
# join fips ids -----------

fips <- read.csv(here::here("usa_adj_mat_files/fips.csv"))

pop$fips_num <- plyr::mapvalues(pop$full_county_name,
                                from = fips$full_name,
                                to = fips$fips_num)


# Note currently some of the fips numbers do not match,
# I think this is mostly the parishes, so this may need to be revisited
# once we move to states with parishes!