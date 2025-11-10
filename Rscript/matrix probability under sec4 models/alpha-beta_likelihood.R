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








ab_likelihood <- function(mat,alpha,beta){
  
  total_prob = 1
  m = dim(mat)[1]
  n_tips = (m+2)/2
  
  for (i in 2:(m - 1)) {
    F_ii <- mat[i , i ]
    
    previous = mat[i-1,i-1]
    diagonals <- valid_heterochronous_diagonals(i, n_tips, previous)
    f_l<- diagonals[1]
    f_u <- diagonals[2]
    prob_ii = ifelse(f_l == f_u , 1, 
                     beta(alpha + (f_u - F_ii)/2, beta+(F_ii-f_l)/2 ) / beta(alpha,beta))
    total_prob = total_prob * prob_ii
  }
  
  # Fill entries off the diagonal and subdiagonal
  for (i in 3:m) {
    for (j in 1:(i - 2)) {
      F_ij = mat[i,j] 
      
      left <- if (j == 1) 0 else mat[i, j - 1]
      left_up <- if (j == 1) 0 else mat[i - 1, j - 1]
      up <- mat[i - 1, j]
      prev_diag <- mat[i - 1, i - 1]
      
      non_diag <- valid_non_diagonal_values(i, j, left, left_up, up, prev_diag, isochronous = FALSE)
      f_l <- non_diag[1]
      f_u <- non_diag[2]
      prob_ij = ifelse(f_l == f_u , 1, 
                       beta(alpha + f_u - F_ij, beta+F_ij-f_l) / beta(alpha,beta))
      total_prob = total_prob * prob_ij
      
    }
  }
  
  return (total_prob)
}


# test= matrix(data=c(2,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1), nrow=4, byrow=TRUE)
# ab_likelihood(test,alpha=10,beta=1)
# beta(10+1,1)/beta(10,1)
# 
# test2 = matrix(data=c(2,0,0,0,1,3,0,0,1,2,2,0,1,1,1,1), nrow=4, byrow=TRUE)
# ab_likelihood(test2,alpha=10,beta=1)
# (beta(10,2) * beta(10,2)*beta(10,2) ) / (beta(10,1) )^3




