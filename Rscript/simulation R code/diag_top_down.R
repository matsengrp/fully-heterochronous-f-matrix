library(gdata)
library(ggplot2)
library(cowplot)

source("bernoulli_splitting_final.R")
source("sampling_coales_final.R")

catalan<-function(n){
  return(choose(2*n,n)/(n+1))
}
# catalan(6)

####uniformly sample dyck path####
dyck_steps <- function(n) {
  stopifnot(n >= 2)
  
  # N(i+1,j)/N(i,j)
  p_right <- function(i, j) {
    ((i - j + 2.0) / (i - j + 1.0)) * ((n - 1.0 - i) / (2.0 * n - 2.0 - i - j))
  }
  
  len <- 2L * n - 3L                  # (1,0) -> (n-1,n-1)
  steps <- integer(len)
  i <- 1L; j <- 0L
  
  for (t in 1:len) {
    if (i == n - 1L) { steps[t] <- -1L; j <- j + 1L; next }
    if (j == n - 1L) { steps[t] <- +1L; i <- i + 1L; next }
    if (j + 1L > i)  { steps[t] <- +1L; i <- i + 1L; next }
    
    p <- p_right(i, j)
    if (p < 0) p <- 0 else if (p > 1) p <- 1
    if (runif(1) < p) { steps[t] <- +1L; i <- i + 1L }
    else              { steps[t] <- -1L; j <- j + 1L }
  }
  steps
}

#### sample one tree under diagonal top-down model####
rhetdiag <- function(n){
  prob= 1/catalan(n-1)
  
  dyck = dyck_steps(n)
  diag = 2+ c(0, cumsum(dyck))
  D = matrix(0, 2*n-2,2*n-2)
  D_diag = c(2, 2*(dyck ==1 ))
  diag(D)=D_diag
  for (i in 1:(2*n-3)){
    non_zero = D[i,]!=0
    F_ii = sum(D[i,non_zero])
    D[i+1,non_zero] = D[i,non_zero]
    choice = which(D[i,]!=0)
    i_j = resample( x=choice,size = 1, prob = D[i,non_zero]/F_ii)
    prob = prob * D[i,i_j] / F_ii
    #print(paste0(i,': ', prob))
    D[i+1, i_j] = D[i,i_j]-1
  }
  F_mat <- t(apply(D, 1L, cumsum)) * (row(D) >= col(D))
  
  #print(prob)
  return(F_mat)
}

#test
rhetdiag(4)


















#### obtain a sample ####
diag_sample <- function(sample_size, tree_size,seed=781){
  set.seed(781)
  internal_length_vec = rep(0,sample_size)
  cherry_count_vec = c(0,sample_size)
  total_length_vec = c(0,sample_size)
  
  for (i in 1:sample_size){
    
    current_tree = rhetdiag(tree_size)
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



# diag_sample_100_5 = diag_sample(100,5)
# diag_sample_100_20 = diag_sample(100,20)
# diag_sample_100_50 = diag_sample(100,50)
# 
# diag_sample_500_5 = diag_sample(500,5)
# diag_sample_500_20 = diag_sample(500,20)
# diag_sample_500_50 = diag_sample(500,50)

diag_sample_1000_5 = diag_sample(1000,5)
diag_sample_1000_20 = diag_sample(1000,20)
diag_sample_1000_50 = diag_sample(1000,50)






#### generating graphs####
# coales
p21 <- ggplot(data.frame(int_length = coales_sample_1000_20$int_length),
              aes(x = int_length)) +
  geom_bar(color = "lightpink", fill = "lightpink") +
  theme_minimal() +
  labs(y="count")+
  theme(axis.title.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_text(size=10),)+
  xlim(34,330)+
  ylim(0,60)

p22 <- ggplot(data.frame(total_length = coales_sample_1000_20$total_length),
              aes(x = total_length)) +
  geom_bar(color = "lightpink", fill = "lightpink") +
  theme_minimal() +
  labs( title= "Coalescence")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10,face="bold"))+
  xlim(70,742)+
  ylim(0,65)


p23 <- ggplot(data.frame(cherry = coales_sample_1000_20$cherry),
              aes(x = cherry)) +
  geom_bar(color = "lightpink", fill = "lightpink") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,17)+
  ylim(0,500)

# diagonal
p24 <- ggplot(data.frame(int_length = diag_sample_1000_20$int_length),
              aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(y="count")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.ticks = element_blank() )+
  xlim(34,330)+
  ylim(0,60)

p25 <- ggplot(data.frame(total_length = diag_sample_1000_20$total_length),
              aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(title="Diagonal top-down")+
  xlim(70,742)+
  ylim(0,65)+
  theme(plot.title = element_text(size = 10,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())



p26 <- ggplot(data.frame(cherry = diag_sample_1000_20$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  xlim(0,17)+
  ylim(0,500)




(p21|p22|p23)/ (p24|p25|p26) /(p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9) & theme(plot.margin = margin(2,2,2,2))





