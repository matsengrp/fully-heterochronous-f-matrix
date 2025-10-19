valid_heterochronous_diagonals <- function(i, n, previous) {
  if (i == 1 || previous == 1) {
    # Return the tuple with the value (2, 2)
    values <- c(2, 2)
  } else if (previous == (2 * n - i - 1)) {
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

sample_heterochronous_f_matrix_from_probs <- function(n_tips, prob_matrix,beta=1) {
  # Matrix side length.
  m <- 2 * n_tips - 2
  
  # Fill constant entries.
  mat <- matrix(0, nrow = m, ncol = m)
  mat[1, 1] <- 2
  mat[2, 1] <- 1
  mat[m, m] <- 1
  
  # Fill the diagonal, which forces values on the subdiagonal.
  for (i in 2:(m - 1)) {
    previous <- mat[i - 1, i - 1]
    diagonals <- valid_heterochronous_diagonals(i, n_tips, previous)
    lower <- diagonals[1]
    upper <- diagonals[2]
    prob <- prob_matrix[i, i]
    mat[i, i] <- ifelse(runif(1) <= prob, lower, upper)
    mat[i + 1, i] <- mat[i, i] - 1
  }
  
  # Fill entries off the diagonal and subdiagonal
  for (i in 3:m) {
    for (j in 1:(i - 2)) {
      left <- if (j == 1) 0 else mat[i, j - 1]
      left_up <- if (j == 1) 0 else mat[i - 1, j - 1]
      up <- mat[i - 1, j]
      prev_diag <- mat[i - 1, i - 1]
      
      non_diag <- valid_non_diagonal_values(i, j, left, left_up, up, prev_diag, isochronous = FALSE)
      lower <- non_diag[1]
      upper <- non_diag[2]
      prob <- prob_matrix[i, j]
      mat[i, j] <- ifelse(rbeta(1,beta,beta) <= prob, lower, upper)
    }
  }
  
  return(mat)
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




bernoulli_sample <- function(sample_size, tree_size,beta=1, seed=781){
  
  set.seed(781)
  
  m <- 2 * tree_size - 2
  prob_matrix <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m) {
      for (j in 1:i) {
        prob_matrix[i, j] <- rbeta(1,beta,beta)
      }
  }
  
  internal_length_vec = rep(0,sample_size)
  cherry_count_vec = c(0,sample_size)
  total_length_vec = c(0,sample_size)
  
  for (i in 1:sample_size){
    
    current_tree = sample_heterochronous_f_matrix_from_probs(tree_size,prob_matrix)
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

bernoulli_100_5=bernoulli_sample(sample_size = 100,tree_size = 5)
bernoulli_100_20=bernoulli_sample(sample_size = 100,tree_size = 20)
bernoulli_100_50=bernoulli_sample(sample_size = 100,tree_size = 50)

bernoulli_500_5=bernoulli_sample(sample_size = 500,tree_size = 5)
bernoulli_500_20=bernoulli_sample(sample_size = 500,tree_size = 20)
bernoulli_500_50=bernoulli_sample(sample_size = 500,tree_size = 50)

bernoulli_1000_5=bernoulli_sample(sample_size = 1000,tree_size = 5)
bernoulli_1000_20=bernoulli_sample(sample_size = 1000,tree_size = 20)
bernoulli_1000_50=bernoulli_sample(sample_size = 1000,tree_size = 50)

# beta = 10
bernoulli_100_5.10=bernoulli_sample(sample_size = 100,tree_size = 5,beta=10)
bernoulli_100_20.10=bernoulli_sample(sample_size = 100,tree_size = 20,beta=10)
bernoulli_100_50.10=bernoulli_sample(sample_size = 100,tree_size = 50,beta=10)

bernoulli_500_5.10=bernoulli_sample(sample_size = 500,tree_size = 5,beta=10)
bernoulli_500_20.10=bernoulli_sample(sample_size = 500,tree_size = 20,beta=10)
bernoulli_500_50.10=bernoulli_sample(sample_size = 500,tree_size = 50,beta=10)

bernoulli_1000_5.10=bernoulli_sample(sample_size = 1000,tree_size = 5,beta=10)
bernoulli_1000_20.10=bernoulli_sample(sample_size = 1000,tree_size = 20,beta=10)
bernoulli_1000_50.10=bernoulli_sample(sample_size = 1000,tree_size = 50,beta=10)

# beta = 100
bernoulli_100_5.100=bernoulli_sample(sample_size = 100,tree_size = 5,beta=100)
bernoulli_100_20.100=bernoulli_sample(sample_size = 100,tree_size = 20,beta=100)
bernoulli_100_50.100=bernoulli_sample(sample_size = 100,tree_size = 50,beta=100)

bernoulli_500_5.100=bernoulli_sample(sample_size = 500,tree_size = 5,beta=100)
bernoulli_500_20.100=bernoulli_sample(sample_size = 500,tree_size = 20,beta=100)
bernoulli_500_50.100=bernoulli_sample(sample_size = 500,tree_size = 50,beta=100)

bernoulli_1000_5.100=bernoulli_sample(sample_size = 1000,tree_size = 5,beta=100)
bernoulli_1000_20.100=bernoulli_sample(sample_size = 1000,tree_size = 20,beta=100)
bernoulli_1000_50.100=bernoulli_sample(sample_size = 1000,tree_size = 50,beta=100)









# beta = 1
p1 <- ggplot(data.frame(int_length = bernoulli_1000_20$int_length),
              aes(x = int_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "internal tree length")+ 
  xlim(34,245)+
  ylim(0,45)

p2 <- ggplot(data.frame(total_length = bernoulli_1000_20$total_length),
              aes(x = total_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "total tree length", title="Bernoulli splitting (beta =1)",y="")+
  xlim(60,546)+
  ylim(0,40)

p3 <- ggplot(data.frame(cherry = bernoulli_1000_20$cherry),
              aes(x = cherry)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "number of cherries",y="")+
  xlim(0,13)+
  ylim(0,300)

# beta=10
p4 <- ggplot(data.frame(int_length = bernoulli_1000_20.10$int_length),
              aes(x = int_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "internal tree length")+ 
  xlim(34,245)+
  ylim(0,45)

p5 <- ggplot(data.frame(total_length = bernoulli_1000_20.10$total_length),
              aes(x = total_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "total tree length", title="Bernoulli splitting (beta =10)",y="")+
  xlim(60,546)+
  ylim(0,40)

p6 <- ggplot(data.frame(cherry = bernoulli_1000_20.10$cherry),
              aes(x = cherry)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "number of cherries",y="")+
  xlim(0,13)+
  ylim(0,300)


# beta =100
p7 <- ggplot(data.frame(int_length = bernoulli_1000_20.100$int_length),
              aes(x = int_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "internal tree length")+ 
  xlim(34,245)+
  ylim(0,45)

p8 <- ggplot(data.frame(total_length = bernoulli_1000_20.100$total_length),
              aes(x = total_length)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "total tree length", title="Bernoulli splitting (beta =100)",y="")+
  xlim(60,546)+
  ylim(0,40)

p9 <- ggplot(data.frame(cherry = bernoulli_1000_20.100$cherry),
              aes(x = cherry)) +
  geom_bar(color = "black", fill = "white") + theme_bw() +
  labs(x = "number of cherries",y="")+
  xlim(0,13)+
  ylim(0,300)

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,
          nrow = 3, ncol = 3, align = "hv")







