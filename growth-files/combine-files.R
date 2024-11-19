# combine files -----


lh <- read.csv(here::here("growth-files/logistic_hier.csv"))
hh <- read.csv(here::here("growth-files/hurdle_hier.csv"))
eh <- read.csv(here::here("growth-files/exp_hier.csv"))
lh$type = "Hierarchical"
hh$type <- "Hierarchical"
eh$type <- "Hierarchical"

ls <- read.csv(here::here("growth-files/logistic_single.csv"))
hs <- read.csv(here::here("growth-files/hurdle_single.csv"))

ls$type <- "Single"
hs$type <- "Single"



total <- rbind(lh,hh,eh,ls,hs)






pp <- ggplot2::ggplot(total,
                      ggplot2::aes(x = state_id,
                                   y = mean)) +
      ggplot2::geom_point(ggplot2::aes(color = model),alpha = 0.7)+
    ggplot2::xlab("State ID")+
    ggplot2::ylab("P(Delta > 0)")+
  hrbrthemes::theme_ipsum() + 
  ggplot2::facet_wrap(~type)+
  ggplot2::geom_hline(yintercept = 0.5)

pp
