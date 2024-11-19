# formula functions ----------------------------


## Add fixed effects ---------------------------

add_fes <- function(init_form = NULL,
                    new_fes){
  
  if(is.null(init_form)==FALSE){
    out <-  paste(init_form,paste(new_fes,collapse = "+"),collapse = "+")
  }else{
    out <-   new_fes
  }
  out
}


## Add single random effect -----------------

add_re <- function(init_form,
                   prior = "list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))",
                   constr = T,
                   name_re,
                   model = "iid",
                   values_re = NULL){
  
  
  if(is.null(values_re)==F){
    mx <- max(values_re,na.rm =T)
  }else{
    mx <- 2
  }
  
  if(is.null(prior) == T){
    prior = "list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))"
  }
  
  if(mx >1){
    re_form <- paste0("f(",
                     name_re,
                     ", model =\"",
                      model,
                     "\",hyper =",
                     prior,
                     ",constr=",
                     constr,
                     ")")
    
    if(is.null(init_form)==FALSE){
      out <-  paste(c(init_form,re_form),collapse = "+")
    }else{
      out <-   re_form
    } 
    
  }else{
    out <- init_form
  }
  out
  
}


## Add many res --------------------------------------

add_res_many <- function(init_form,
                         re_names,
                         prior_list = NULL,
                         constr_vec = NULL,
                         model_vec = NULL,
                         values_list = NULL){
  
  
  P <- length(re_names)
  
  if(is.null(prior_list) == T){
    prior_list = vector("list",length = P)
  }
  
  if(is.null(constr_vec) == T){
    constr_vec = rep(NULL,P)
  }
  
  if(is.null(model_vec) == T){
    model_vec = rep(NULL,P)
  }
  
  if(is.null(values_list) == T){
    values_list = vector("list",length = P)
  }
  
  
  for(i in 1:P){
    
    init_form <-  add_re(init_form = init_form,
                         prior = prior_list[[i]],
                         constr = constr_vec[i],
                         name_re = re_names[i],
                         model = model_vec[i],
                         values_re = values_list[[i]])
    
  }
  
  init_form  
}


## Add single spatial -------------------------------------
spat_single <- function(form = NULL,
                        id_base,
                        group_base = 'm',
                        spatial_model = 'bym',
                        grouped_type = 'ar',
                        ar_order = 5,
                        scale_model = T,
                        adjust_con = T,
                        copy = F,
                        i,
                        copy_id = NULL){
  
  if(copy == F){
    inter<- paste0("f(",
                   id_base,
                   i,
                   ",E",
                   i,
                   ",model='",
                   spatial_model,
                   "',graph=as.matrix(adj.mat),group =",
                   group_base,
                   i,
                   ",scale.model=",
                   scale_model,
                   ",adjust.for.con.comp = ",
                   adjust_con,
                   ",control.group = list(model = '",
                   grouped_type,
                   "', order = ",
                   ar_order,
                   ",hyper = hyper.prec,scale.model = TRUE),initial = 3)")
    
  }else{
    
    inter <- paste0("f(",
                    id_base,
                    i,
                    ",E",
                    i,
                    ",copy = '",
                    copy_id,
                    "',group = ",
                    group_base,
                    i,
                    ",hyper = list(beta = beta.prior))")
    
    
  }
  
  
  
  inter <- stringr::str_remove_all(inter,"\\n")  
  
  if(is.null(form)==FALSE){ 
    out <- paste0(c(form,
                    inter),
                  collapse = "+")
    
  }else{
    out <- inter
  }
  
  out
  
}


# Add all spatial ----------------


add_spatial <- function(init_form,
                        id_base,
                        spatial_model = 'bym',
                        group_base = 'm',
                        n_groups,
                        grouped_type = 'ar',
                        ar_order = 5,
                        scale_model = T,
                        adjust_con = T){
  
  
  start <- spat_single(form = init_form,
                       id_base = id_base,
                       spatial_model = spatial_model,
                       group_base = group_base,
                       grouped_type = grouped_type,
                       scale_model = scale_model,
                       adjust_con = adjust_con,
                       copy = F,
                       i = 1,
                       ar_order = ar_order)
  
  
  for(i in 2:n_groups){
    start <- spat_single(form = start,
                         id_base = id_base,
                         spatial_model = spatial_model,
                         group_base = group_base,
                         grouped_type = grouped_type,
                         scale_model = scale_model,
                         adjust_con = adjust_con,
                         copy = T,
                         i = i,
                         copy_id = paste0(id_base,1),
                         ar_order = ar_order)
  }
  
  start
  
}

## Add spatial copies only -----------------------------

add_spatial_copies <- function(init_form,
                               id_base,
                               group_base,
                               copy = T,
                               which_groups,
                               copy_id = "id1"){
  
  n_groups <- length(which_groups)
  
  start <- spat_single(form = init_form,
                       id_base = id_base,
                       spatial_model = T,
                       group_base = group_base,
                       grouped_type = 'ar2',
                       scale_model = scale_model,
                       adjust_con = adjust_con,
                       copy = T,
                       i = which_groups[1],
                       copy_id = copy_id)
  
  if(length(which_groups) >1){
    for(i in 2:n_groups){
      start <- spat_single(form = start,
                           id_base = id_base,
                           spatial_model = spatial_model,
                           group_base = group_base,
                           grouped_type = 'ar2',
                           scale_model = T,
                           adjust_con = T,
                           copy = T,
                           i = which_groups[i],
                           copy_id = copy_id)
    }
  }
  
  start 
}


# add ar effect ------------

add_ar <- function(init_form,
                   prior = NULL,
                   id_base = "m",
                   which_ind,
                   ar_order = 5,
                   copy =F,
                   copied_model = "m1"){
  
  if(is.null(prior)){
    prior <- "hyper.prec"
  }

  if(copy == FALSE){
  out <- paste0(init_form,
                  "+f(",
                  id_base,
                  which_ind,
                  ",model = 'ar', order = ",
                  ar_order,
                  ", adjust.for.con.comp = TRUE, hyper =",
                  prior,
                  ")")
   
  }else{
   out <- paste0(init_form,
                   "+f(",
                   id_base,
                   which_ind,
                   ",copy ='",
                   copied_model,
                   "',hyper = list(beta = beta.prior))")
    

  }
  
  out
  
}


add_ar_many <- function(init_form,
                        prior = NULL,
                        id_base = "m",
                        which_inds,
                        ar_order = 5,
                        copy_vec,
                        copied_model = "m1"){
  
 nl <- length(which_inds) 
  
for(i in 1:nl){
  if(i ==1){
    cur_form <- init_form
  }
  
  cur_form <- add_ar(init_form = cur_form,
                prior = prior,
                id_base = id_base,
                which_ind = which_inds[i],
                ar_order = ar_order,
                copy =copy_vec[i],
                copied_model = copied_model)
}
  
  
cur_form  
  
}




# Full joint formula ----------------------------------

joint_formula <- function(outcome = "y",
                          rm_int = T,
                          fixed_effects,
                          re_effects,
                          full_spatial =T,
                          lag_spatial_hosp = T,
                          lag_spatial_deaths = T,
                          n_outcomes_cases,
                          n_outcomes_hosp,
                          n_outcomes_deaths,
                          ar_order = 5
){
  
  # Start formulat
  
  init_form <- paste0(outcome, " ~ ")
  
  if(rm_int ==T){
    init_form <-  paste0(init_form,"-1 +")
  }
  
  init_form <- add_fes(init_form = init_form,
                       new_fes = fixed_effects)
  
  
  n_re <- length(re_effects)
  init_form <- add_res_many(init_form = init_form,
                            re_names = re_effects,
                            prior_list = NULL,
                            constr_vec = rep(T,n_re),
                            model_vec = rep("iid",n_re),
                            values_list = NULL)
  
  
  if(full_spatial == T){
  n_groups <- sum(n_outcomes_cases,n_outcomes_hosp,n_outcomes_deaths)
  
  init_form <- add_spatial(init_form = init_form,
                                       id_base = "id",
                                       spatial_model = 'bym2',
                                       group_base = 'm',
                                       n_groups = n_groups,
                                       grouped_type = 'ar',
                                       ar_order = ar_order,
                                       scale_model = T,
                                       adjust_con = T)
  
  

  }
  
  if(lag_spatial_hosp == T){
    
    which_hosp <- c((n_outcomes_cases +1):(n_outcomes_cases +n_outcomes_hosp))
    
    
    
    init_form <- add_spatial_copies(init_form = init_form,
                                    id_base = "idc",
                                    group_base = "mlag",
                                    copy = T,
                                    which_groups = which_hosp,
                                    copy_id = "id1")
  }
  
  if( lag_spatial_deaths == T){
    
    which_deaths <- c((n_outcomes_cases +n_outcomes_hosp+1):(n_outcomes_cases +n_outcomes_hosp+n_outcomes_deaths))
    
    init_form <- add_spatial_copies(init_form = init_form,
                                    id_base = "idc",
                                    group_base = "mlag",
                                    copy = T,
                                    which_groups = which_deaths,
                                    copy_id = "id1")
    
    
  }
  
  init_form
  
}


formula_county <- function(outcome = "y",
                           rm_int = T,
                           fixed_effects,
                           re_effects,
                           lag_spatial_hosp = T,
                           lag_spatial_deaths = T,
                           n_outcomes_cases,
                           n_outcomes_hosp,
                           n_outcomes_deaths){
  
  
  
  init_form <- paste0(outcome, " ~ ")
  
  if(rm_int ==T){
    init_form <-  paste0(init_form,"-1 +")
  }
  
  init_form <- add_fes(init_form = init_form,
                       new_fes = fixed_effects)
  
  
  n_re <- length(re_effects)
  init_form <- add_res_many(init_form = init_form,
                            re_names = re_effects,
                            prior_list = NULL,
                            constr_vec = rep(T,n_re),
                            model_vec = rep("iid",n_re),
                            values_list = NULL)
  
  
  nl <- sum(n_outcomes_cases,n_outcomes_hosp,n_outcomes_deaths)
  
  init_form <- add_ar_many(init_form = init_form,
                           prior = NULL,
                           id_base = "m",
                           which_inds = c(1:nl),
                           ar_order = 5,
                           copy_vec = c(F,rep(T,(nl-1))),
                           copied_model = "m1")
  
  if(lag_spatial_hosp == T){
    
    which_hosp <- c((n_outcomes_cases +1):(n_outcomes_cases +n_outcomes_hosp))
    
    
    
    init_form <- add_ar_many(init_form = init_form,
                             prior = NULL, 
                             id_base = "mlag",
                             which_inds = which_hosp,
                             copy_vec = rep(T,length(which_hosp)),
                             copied_model = "m1")
  }

  
  if( lag_spatial_deaths == T){
    
    which_deaths <- c((n_outcomes_cases +n_outcomes_hosp+1):(n_outcomes_cases +n_outcomes_hosp+n_outcomes_deaths))
    
    init_form <- add_ar_many(init_form = init_form,
                             prior = NULL, 
                             id_base = "mlag",
                             which_inds = which_deaths,
                             copy_vec = rep(T,length(which_deaths)),
                             copied_model = "m1")
    
    
  }  
  
  
init_form  
  
}
