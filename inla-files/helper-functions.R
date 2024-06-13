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
    
  
  
  
  

                    



  
  
  
  
  
  
  