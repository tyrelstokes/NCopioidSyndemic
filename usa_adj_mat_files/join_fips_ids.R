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

adj_names <- tolower(adj_mat$X)
fips_names <- tolower(fips_dt$full_name)



adj_names <- stringr::str_replace_all(adj_names," ","")
adj_names <- stringr::str_replace_all(adj_names,"\\.","")
adj_names <- stringr::str_replace_all(adj_names,"county","")
adj_names <- stringr::str_replace_all(adj_names,"parish","")

fips_names <- stringr::str_replace_all(fips_names," ","")
fips_names <- stringr::str_replace_all(fips_names,"\\.","")
fips_names <- stringr::str_replace_all(fips_names,"county","")
fips_names <- stringr::str_replace_all(fips_names,"parish","")



matched_counties <- fips_names[which(fips_names%in% adj_names)]
unmatched_counties <- fips_names[-which(fips_names%in% adj_names)]

ff <- data.frame(adj_names)

ff2 <- data.frame(unmatched_counties)

mps <- plyr::mapvalues(adj_names,
                       from = fips_names,
                       to = fips_dt$fips_num)


adj_mat$county_id <- mps


length(fips_dt$full_name[which(fips_dt$full_name %in% adj_mat$X)])

write.csv(adj_mat,here::here("usa_adj_mat_files/aug_adj_mat.csv"))
