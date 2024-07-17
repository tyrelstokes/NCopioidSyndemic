# Find the appropriate indexes (two vectors) ------------

index_finder <- function(n,
                         np){
  
 beg <- 1
 end <- n
 max_int <- n*np

 beg_vec <- vector(length = np)
 end_vec <- vector(length = np)
 
 beg_vec[1] <- 1
 end_vec[1] <- n
 j <- 2
 while(end < max_int){
   beg <- end + 1
   end <- end + n
   beg_vec[j] <- beg
   end_vec[j] <- end
   j <- j+1
 }
 out <- list(beg_vec = beg_vec,
             end_vec = end_vec)
 
 out
}


# Create INLA vectors with appropriate NAs -----------

inla_vectors <- function(n,
                        np,
                        x_list,
                        base_name = "y"){
  indices <- index_finder(n = n,
                          np = np)
  beg_vec <- indices$beg_vec
  end_vec <- indices$end_vec
  
  
  out <- vector("list", length = np)
  
  full_vec <- c(1:max(end_vec))
  for(i in 1:np){
    inter <- vector(length = n*np)
    
    beg <- beg_vec[i]
    end <- end_vec[i]
    
    
    
    inter[beg:end] <- x_list[[i]]
    
    na_ind <- full_vec[-(beg:end)]
    
    inter[na_ind] <- rep(NA,n)
    
  out[[i]] <- inter
  }
  
  names(out) <- paste0(base_name,"_",c(1:np))
 out 
}  


# Inla matrix function --------------

inla_matrix <- function(n,
                        np,
                        y_list){
  
  y <- matrix(ncol = np, nrow = n*np)
  
  indices <- index_finder(n = n,
                          np = np)
  beg_vec <- indices$beg_vec
  end_vec <- indices$end_vec
  
  full_vec <- c(1:max(end_vec))
  
  for(i in 1:np){
    inter <- vector(length = n*np)
    
    beg <- beg_vec[i]
    end <- end_vec[i]
    
    
    
    inter[beg:end] <- y_list[[i]]
    
    na_ind <- full_vec[-(beg:end)]
    
    inter[na_ind] <- rep(NA,n)
    
    y[,i] <- inter
  }
  
  y
}
#############

# inla formula function -------------------

inla_formula <- function(np,
                         spatial = T,
                         spatial_model = "'bym2'",
                         grouped_effect = T,
                         grouped_type = "'ar2'",
                         temporal_model = "ar",
                         order = 5,
                         mar = T,
                         vac = T,
                         late = T){
 
if(np == 3){ 
form <-  "y ~ -1+ int1 + int2 +int3+"

if(mar == T){
form <- paste0(form,
               "mar1 + mar2 + mar3+")

}

if(vac == T){
  
  form <- paste0(form,
                 "vac1 + vac2 + vac3+") 
}
 
if(late == T){
  
  form <- paste0(form,
                 "late1 + late2 + late3+") 
}


if(spatial == T){
  if(grouped_effect == T){
  
  form <- paste0(form,
                 "f(id1,E1,model=",spatial_model,
      ",graph=as.matrix(adj.mat),group = m1,
      scale.model=TRUE,
      adjust.for.con.comp = TRUE,",
      "control.group = list(model = ",
      grouped_type,
      ",hyper = hyper.prec,
                           scale.model = TRUE),
      initial = 3)+
    f(id2,
      E2,
      copy = 'id1',
      group = m2,
      hyper = list(beta = beta.prior))+
    f(id3,
      E3,
      copy = 'id1',
      group = m3,
      hyper = list(beta = beta.prior)+", sep = "")
  
  form <- stringr::str_remove_all(form,"\\n")

  
  }else{
    
    form <- paste0(form,
                   "f(id1,E1,model=",spatial_model,
                   ",graph=as.matrix(adj.mat),group = m1,
      scale.model=TRUE,
      adjust.for.con.comp = TRUE,
      initial = 3)+
    f(id2,
      E2,
      copy = 'id1',
      group = m2,
      hyper = list(beta = beta.prior))+
    f(id3,
      E3,
      copy = 'id1',
      group = m3,
      hyper = list(beta = beta.prior)+", sep = "")
    
    form <- stringr::str_remove_all(form,"\\n") 
  }
  
  
}


if(temporal_model == "ar"){
  
 form <- paste0(form, "f(m1,
    model= 'ar',
    order = ",
    order,",adjust.for.con.comp = TRUE)+
    f(m2,
      copy = 'm1',
      hyper = list(beta = beta.prior))+
    f(m3,
      copy = 'm1',
      hyper = list(beta = beta.prior))+")
  
 form <- stringr::str_remove_all(form,"\\n") 
  
}



  
}else{
 
  form <-  "y ~ -1+ int1 + int2+"
  
  if(mar == T){
    form <- paste0(form,
                   "mar1 + mar2+")
    
  }
  
  if(vac == T){
    
    form <- paste0(form,
                   "vac1 + vac2+") 
  }
  
  if(late == T){
    
    form <- paste0(form,
                   "late1 + late2+") 
  }
  
  
  
  if(spatial == T){
    if(grouped_effect == T){
      
      form <- paste0(form,
                     "f(id1,E1,model=",spatial_model,
                     ",graph=as.matrix(adj.mat),group = m1,
      scale.model=TRUE,
      adjust.for.con.comp = TRUE,",
                     "control.group = list(model = ",
                     grouped_type,
                     ",hyper = hyper.prec,
                           scale.model = TRUE),
      initial = 3)+
    f(id2,
      E2,
      copy = 'id1',
      group = m2,
      hyper = list(beta = beta.prior))+", sep = "")
      
      form <- stringr::str_remove_all(form,"\\n")
      
      
    }else{
      
      form <- paste0(form,
                     "f(id1,E1,model=",spatial_model,
                     ",graph=as.matrix(adj.mat),group = m1,
      scale.model=TRUE,
      adjust.for.con.comp = TRUE,
      initial = 3)+
    f(id2,
      E2,
      copy = 'id1',
      group = m2,
      hyper = list(beta = beta.prior))+", sep = "")
      
      form <- stringr::str_remove_all(form,"\\n") 
    }
    
    
  }
  
  
  if(temporal_model == "ar"){
    
    form <- paste0(form, "f(m1,
    model= 'ar',
    order = ",
                   order,",adjust.for.con.comp = TRUE)+
    f(m2,
      copy = 'm1',
      hyper = list(beta = beta.prior))+")
    
    form <- stringr::str_remove_all(form,"\\n") 
    
  }
  
  
}
  
  lchar_plus <- (stringr::str_sub(form,-1) =="+")
  if(lchar_plus == TRUE){
    form <- stringr::str_sub(form,end = -2)
  }
  
  form
}  
  
  

                    



  
  
  
  
  
  
  