# load all data ---------------
adj_mat <- read.csv(here::here("usa_adj_mat_files/countyadj.csv"))

fips_dt <- readxl::read_excel(here::here("usa_adj_mat_files/US_FIPS_Codes.xls"),
                              skip = 1)

sts <- read.csv(here::here("usa_adj_mat_files/states.csv"))


# manipulate the fips dt to match on adj mat -------
fips_dt$fips_num <- paste0(fips_dt$`FIPS State`,
                           fips_dt$`FIPS County`)

fips_dt$abbrev <- plyr::mapvalues(fips_dt$State,
                                  from = sts$State,
                                  to = sts$Abbreviation)


fips_dt$full_name <- paste0(fips_dt$`County Name`,
                            " County ",
                            fips_dt$abbrev)


write.csv(fips_dt,here::here("usa_adj_mat_files/fips.csv"))

# join on adj_mat --------



