library(patchwork)

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

sample_heterochronous_f_matrix_from_probs <- function(n_tips, prob_matrix,alpha=1,beta=1) {
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
      mat[i, j] <- ifelse(rbeta(1,alpha,beta) <= prob, lower, upper)
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




bernoulli_sample <- function(sample_size, tree_size,alpha=1,beta=1, seed=781){
  
  set.seed(781)
  
  m <- 2 * tree_size - 2
  prob_matrix <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m) {
      for (j in 1:i) {
        prob_matrix[i, j] <- rbeta(1,alpha,beta)
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

# alpha = 10, beta =1
bernoulli_100_5.10_1=bernoulli_sample(sample_size = 100,tree_size = 5,alpha=10)
bernoulli_100_20.10_1=bernoulli_sample(sample_size = 100,tree_size = 20,alpha=10)
bernoulli_100_50.10_1=bernoulli_sample(sample_size = 100,tree_size = 50,alpha=10)

bernoulli_500_5.10_1=bernoulli_sample(sample_size = 500,tree_size = 5,alpha=10)
bernoulli_500_20.10_1=bernoulli_sample(sample_size = 500,tree_size = 20,alpha=10)
bernoulli_500_50.10_1=bernoulli_sample(sample_size = 500,tree_size = 50,alpha=10)

bernoulli_1000_5.10_1=bernoulli_sample(sample_size = 1000,tree_size = 5,alpha=10)
bernoulli_1000_20.10_1=bernoulli_sample(sample_size = 1000,tree_size = 20,alpha=10)
bernoulli_1000_50.10_1=bernoulli_sample(sample_size = 1000,tree_size = 50,alpha=10)

# alpha = 10, beta = 10
bernoulli_100_5.10_10=bernoulli_sample(sample_size = 100,tree_size = 5,beta=10,alpha=10)
bernoulli_100_20.10_10=bernoulli_sample(sample_size = 100,tree_size = 20,beta=10,alpha=10)
bernoulli_100_50.10_10=bernoulli_sample(sample_size = 100,tree_size = 50,beta=10,alpha=10)

bernoulli_500_5.10_10=bernoulli_sample(sample_size = 500,tree_size = 5,beta=10,alpha=10)
bernoulli_500_20.10_10=bernoulli_sample(sample_size = 500,tree_size = 20,beta=10,alpha=10)
bernoulli_500_50.10_10=bernoulli_sample(sample_size = 500,tree_size = 50,beta=10,alpha=10)

bernoulli_1000_5.10_10=bernoulli_sample(sample_size = 1000,tree_size = 5,beta=10,alpha=10)
bernoulli_1000_20.10_10=bernoulli_sample(sample_size = 1000,tree_size = 20,beta=10,alpha=10)
bernoulli_1000_50.10_10=bernoulli_sample(sample_size = 1000,tree_size = 50,beta=10,alpha=10)

# alpha=1, beta=10
bernoulli_100_5.1_10=bernoulli_sample(sample_size = 100,tree_size = 5,beta=10)
bernoulli_100_20.1_10=bernoulli_sample(sample_size = 100,tree_size = 20,beta=10)
bernoulli_100_50.1_10=bernoulli_sample(sample_size = 100,tree_size = 50,beta=10)

bernoulli_500_5.1_10=bernoulli_sample(sample_size = 500,tree_size = 5,beta=10)
bernoulli_500_20.1_10=bernoulli_sample(sample_size = 500,tree_size = 20,beta=10)
bernoulli_500_50.1_10=bernoulli_sample(sample_size = 500,tree_size = 50,beta=10)

bernoulli_1000_5.1_10=bernoulli_sample(sample_size = 1000,tree_size = 5,beta=10)
bernoulli_1000_20.1_10=bernoulli_sample(sample_size = 1000,tree_size = 20,beta=10)
bernoulli_1000_50.1_10=bernoulli_sample(sample_size = 1000,tree_size = 50,beta=10)












# beta = 1
df=data.frame(int_length = bernoulli_1000_20.10_1$int_length)
hist_data <- df %>%
  mutate(bin = cut(int_length, breaks = 6)) %>%
  group_by(bin) %>%
  summarise(
    int_length = mean(int_length, na.rm = TRUE),
    count = n()
  )
p1 <- ggplot(data.frame(int_length = bernoulli_1000_20.10_1$int_length),
              aes(x = int_length)) +
  geom_histogram(bins = 6, color = "gray80", fill = "gray80") +
  geom_line(data=hist_data,
            aes(y = count,x=int_length),           # <- forces 0 outside [min(x), max(x)]
            linewidth = 0.5, color = "blue"
  ) +
  theme_minimal() +
  labs(y="count")+
  theme(axis.title.x = element_blank(), 
        axis.ticks = element_blank() )

df=data.frame(total_length = bernoulli_1000_20.10_1$total_length)
hist_data <- df %>%
  mutate(bin = cut(total_length, breaks = 15)) %>%
  group_by(bin) %>%
  summarise(
    total_length = mean(total_length, na.rm = TRUE),
    count = n()
  )
p2 <- ggplot(data.frame(total_length = bernoulli_1000_20.10_1$total_length),
              aes(x = total_length)) +
  geom_histogram(bins = 15, color = "gray80", fill = "gray80") +
  geom_line(data=hist_data,
    aes(y = count,x=total_length),           # <- forces 0 outside [min(x), max(x)]
    linewidth = 0.5, color = "blue"
  ) + 
  theme_minimal() +
  labs(x = "total tree length", title=expression("Bernoulli splitting (" ~alpha == 10~","~ beta==1 ~ ")"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.ticks = element_blank() )+
  ylim(-10,570)


df=data.frame(total_length = bernoulli_1000_20.10_1$cherry)
hist_data <- df %>%
  mutate(bin = cut(total_length, breaks = 4)) %>%
  group_by(bin) %>%
  summarise(
    total_length = mean(total_length, na.rm = TRUE),
    count = n()
  )
p3 <- ggplot(data.frame(cherry = bernoulli_1000_20.10_1$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray", fill = "gray") +
  geom_line(data=hist_data,
            aes(y = count,x=total_length),           # <- forces 0 outside [min(x), max(x)]
            linewidth = 0.5, color = "blue"
  ) + 
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank() )

# beta=10
p4 <- ggplot(data.frame(int_length = bernoulli_1000_20.10_10$int_length),
              aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "white") +
  geom_line(
    aes(y = after_stat(density) * after_stat(n) * 1.5),
    stat = "density",
    trim = TRUE,              # <- forces 0 outside [min(x), max(x)]
    linewidth = 0.5, color = "blue"
  ) + 
  theme_minimal() +
  labs(x = "internal tree length",y="count")+
  theme(axis.title.x = element_blank(), 
        axis.ticks = element_blank() )+
  xlim(34,330)+
  ylim(0,60)

p5 <- ggplot(data.frame(total_length = bernoulli_1000_20.10_10$total_length),
              aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  geom_line(
    aes(y = after_stat(density) * after_stat(n) * 2.5),
    stat = "density",
    trim = TRUE,              # <- forces 0 outside [min(x), max(x)]
    linewidth = 0.5, color = "blue"
  )+
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "total tree length", title=expression("Bernoulli splitting (" ~alpha == 10~","~ beta==10 ~ ")"))+
  xlim(70,742)+
  ylim(0,65)+
  theme(plot.title = element_text(size = 10),
        axis.ticks = element_blank() )

df=data.frame(total_length = bernoulli_1000_20.10_10$cherry)
hist_data <- df %>%
  mutate(bin = cut(total_length, breaks = 17)) %>%
  group_by(bin) %>%
  summarise(
    total_length = mean(total_length, na.rm = TRUE),
    count = n()
  )

p6 <- ggplot(data.frame(cherry = bernoulli_1000_20.10_10$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  geom_line(data=hist_data,
            aes(y = count,x=total_length),           # <- forces 0 outside [min(x), max(x)]
            linewidth = 0.5, color = "blue"
  ) + 
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+ 
  xlim(0,17)+
  ylim(0,500)


# beta =100
p7 <- ggplot(data.frame(int_length = bernoulli_1000_20.1_10$int_length),
              aes(x = int_length)) +
  geom_bar(color = "gray", fill = "gray") +
  geom_line(
    aes(y = after_stat(density) * after_stat(n) * 1.1),
    stat = "density",
    trim = TRUE,              # <- forces 0 outside [min(x), max(x)]
    linewidth = 0.5, color = "blue"
  )+
  theme_minimal() +
  labs(x = "internal tree length",y="count")+
  xlim(34,330)+
  ylim(0,60)+
  theme(axis.ticks = element_blank() )

p8 <- ggplot(data.frame(total_length = bernoulli_1000_20.1_10$total_length),
              aes(x = total_length)) +
  geom_bar(color = "gray", fill = "gray") +
  geom_line(
    aes(y = after_stat(density) * after_stat(n) * 2.2),
    stat = "density",
    trim = TRUE,              # <- forces 0 outside [min(x), max(x)]
    linewidth = 0.5, color = "blue"
  )+
  theme_minimal() +
  labs(x = "total tree length", title=expression("Bernoulli splitting (" ~alpha == 1~","~ beta==10 ~ ")"))+
  xlim(70,742)+
  ylim(0,65)+
  theme(plot.title = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )

df=data.frame(total_length = bernoulli_1000_20.1_10$cherry)
hist_data <- df %>%
  mutate(bin = cut(total_length, breaks = 17)) %>%
  group_by(bin) %>%
  summarise(
    total_length = mean(total_length, na.rm = TRUE),
    count = n()
  )
p9 <- ggplot(data.frame(cherry = bernoulli_1000_20.1_10$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") + 
  geom_line(data=hist_data,
            aes(y = count,x=total_length),           # <- forces 0 outside [min(x), max(x)]
            linewidth = 0.5, color = "blue"
  ) + 
  theme_minimal() +
  labs(x = "number of cherries")+ 
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,17)+
  ylim(0,500)




(p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9) & theme(plot.margin = margin(2,2,2,2))


