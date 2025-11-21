library(patchwork)
library(ggplot2)

valid_heterochronous_diagonals <- function(i, n, previous) {
  if (i == 1 || previous == 1) {
    # Return the tuple with the value (2, 2)
    values <- c(2, 2)
  } else if (previous >= (2 * n - i - 1)) {
    # Return the tuple with the value (previous - 1, previous - 1)
    values <- c(previous - 1, previous - 1)
  } else {
    # Return the tuple with the values (previous - 1, previous + 1)
    values <- c(previous - 1, previous + 1)
  }
  return(values)
}

valid_non_diagonal_values <- function(i, j, left, left_up, up, prev_diag = NULL, isochronous = FALSE) {
  if (j == 1) {
    # Inequalities for first column.
    lower <- max(0, up - 1)
    upper <- up
  } else {
    # Generic inequalities.
    lower <- max(left, up - 1, left + up - left_up - 1)
    upper <- min(up, left + up - left_up)
  }
  
  if (!isochronous && up == prev_diag) {
    forced_value <- prev_diag - 1
    if (lower <= forced_value && forced_value <= upper) {
      lower <- forced_value
      upper <- forced_value
    } else {
      stop(sprintf("Problem with entry (%d,%d) with forced value %d; bounds %d and %d; neighboring entries %d, %d, %d; and previous diagonal value %d appearing.",
                   i, j, forced_value, lower, upper, left, left_up, up, prev_diag))
    }
  }
  
  if (lower >= 0 && lower <= upper) {
    return(c(lower, upper))
  } else {
    stop(sprintf("Problem with entry (%d,%d) with bounds %d and %d, and neighboring entries %d, %d, %d.",
                 i, j, lower, upper, left, left_up, up))
  }
}

# beta2_prob <- function(n_tips,alpha,beta){
#   
#   m <- 2 * n_tips - 2
#   
#   probs = matrix(1,nrow = m, ncol=m)
#   for (i in 1:(m-1)) {
#     probs[i, i] <- rbeta(1,alpha,beta)
#     for (j in (i+1):m) {
#       print(paste0("i:", i, "j:" ,j))
#       prev = probs[j-1,i]
#       probs[j, i] <- prev + (1-prev)* rbeta(1,alpha,beta)
#     }
#   }
#   
#   return(probs)
# }

# beta2_prob(4,1.5,1.5)
# this is wrong because the prob of min is always growing, if we choose max, then we want the 
# prob of max to grow

# instead, we will do it in the sampling function
# first fix prob on the diagonal, then sample the diag (based on sampling/bifurcating, 
# we know if we need max or min)
# then we fill in each row (move one by one......?????????

#maybe we need to do the D-matrix and then transfer into the F-matrix => this feels better

sample_tree_beta2 <- function(n_tips, alpha,beta) {
  m <- 2 * n_tips - 2
  
  # Create F-matrix & Fill constant entries
  F_mat <- matrix(0, nrow = m, ncol = m)
  F_mat[1, 1] <- 2
  F_mat[2, 1] <- 1
  F_mat[m, m] <- 1
  
  # Create P-matrix
  prob_mat <- matrix(Inf, nrow=m, ncol=m )
  
  # Fill the diagonal, which forces values on the subdiagonal.
  for (i in 2:(m - 1)) {
    previous <- F_mat[i - 1, i - 1]
    diagonals <- valid_heterochronous_diagonals(i, n_tips, previous)
    lower <- diagonals[1]
    upper <- diagonals[2]
    
    prob <- rbeta(1,alpha, beta)
    prob_mat[i,i] = ifelse(lower==upper, Inf, prob)
    
    F_mat[i, i] <- ifelse(runif(1) <= prob, lower, upper)
    F_mat[i + 1, i] <- F_mat[i, i] - 1
    
    # D_mat[i,i] <- (F_mat[i,i] - previous)+1
  }
  
  
  
  # Fill entries off the diagonal and subdiagonal
  for (i in 3:m) {
    for (j in 1:(i - 2)) {
      left <- if (j == 1) 0 else F_mat[i, j - 1]
      left_up <- if (j == 1) 0 else F_mat[i - 1, j - 1]
      up <- F_mat[i - 1, j]
      prev_diag <- F_mat[i - 1, i - 1]
      
      non_diag <- valid_non_diagonal_values(i, j, left, left_up, up, prev_diag, isochronous = FALSE)
      lower <- non_diag[1]
      upper <- non_diag[2]
      
      current_prob <- ifelse(upper == lower, Inf, rbeta(1,alpha,beta))
      up_prob <- prob_mat[i-1,j]
      
      prob_mat[i,j] = ifelse(up_prob == Inf, current_prob, up_prob *current_prob)
      prob = prob_mat[i,j]
      
      F_mat[i, j] <- ifelse(runif(1) <= prob, lower, upper)
    }
  }
  
  return(F_mat)
}




# sample_example <- function(n_samples) {
#   n_tips <- 5
#   m <- 2 * n_tips - 2
#   prob_matrix <- matrix(0, nrow = m, ncol = m)
#   for (i in 1:m) {
#     for (j in 1:i) {
#       prob_matrix[i, j] <- runif(1)
#     }
#   }
#   
#   for (sample in 1:n_samples) {
#     mat <- sample_heterochronous_f_matrix_from_probs(n_tips, prob_matrix)
#     cat(sprintf("Sampled F-matrix:\n %s\n", mat))
#   }
# }




beta2_sample <- function(sample_size, tree_size,alpha,beta, seed=781){
  
  set.seed(781)
  
  m <- 2 * tree_size - 2
  
  internal_length_vec = rep(0,sample_size)
  cherry_count_vec = c(0,sample_size)
  total_length_vec = c(0,sample_size)
  
  for (i in 1:sample_size){
    
    current_tree = sample_tree_beta2(tree_size,alpha,beta)
    r = 2*tree_size-2
    internal_tree_length = 0 
    total_cherry = 0
    
    D = t(diff(t(cbind(rep(0,r),current_tree))))
    D [upper.tri (D )]=0 
    #padded_D = rbind(D,rep(0,r))
    E = diff(rbind(D,rep(0,r)))*-1
    E[upper.tri(E)]=0
    
    # count number of cherries
    padded_E = cbind(E,rep(0,r))
    for ( j in which(diag(D)==2)){
      #branches = E [, which( c(0, E[,j])==1) ] 
      if (sum( padded_E [, which( c(0, E[,j])==1) ] !=0)==0 ){
        total_cherry=total_cherry+1 
      }
    }
    #print(total_cherry)
    cherry_count_vec[i] = total_cherry
    
    # internal length 
    for (k in which(diag(D)[2:r]==2)+1){
      branch_length =k - which( (D[k-1,] - D[k,])==1)[1]
      internal_tree_length= internal_tree_length + branch_length 
    }
    internal_length_vec[i] = internal_tree_length
    
    total_length_vec[i] = sum(diag(current_tree))
    
    
  }# end of iterating sample size
  
  
  print(paste0("avg number of cherries: ", sum(cherry_count_vec)/sample_size))
  print(paste0("avg internal length: ", sum(internal_length_vec) / sample_size))
  print(paste0("avg total length: ", sum(total_length_vec) / sample_size))
  
  
  return( list(
    int_length = internal_length_vec,
    cherry = cherry_count_vec,
    total_length = total_length_vec, 
    avg_int_length = sum(internal_length_vec) / sample_size, 
    avg_total_length = sum(total_length_vec) / sample_size,
    avg_cherry = sum(cherry_count_vec)/sample_size
  ))
  
}


