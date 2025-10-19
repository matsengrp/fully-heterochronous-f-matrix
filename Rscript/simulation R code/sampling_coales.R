library(philentropy)
library("ape")
library(gdata)
library(gmp)
####(not really used) function for reduced tangent number####
tangent_numbers_std <- function(n) {
  stopifnot(n >= 1L)
  T <- as.bigz(rep(0, n))
  T[1] <- as.bigz(1)
  if (n >= 2) {
    for (k in 2:n) T[k] <- as.bigz(k - 1) * T[k - 1]         # init: (k-1)!
    for (k in 2:n) {
      for (j in k:n) {
        T[j] <- as.bigz(j - k) * T[j - 1] + as.bigz(j - k + 2) * T[j]
      }
    }
  }
  T
}

# Reduced tangent numbers: T_n = T_std_n / 2^(n-1)
reduced_tangent_numbers <- function(n) {
  Tstd <- tangent_numbers_std(n)
  pow2 <- pow.bigz(as.bigz(2), 0:(n-1))
  Tred <- Tstd / pow2  # exact divition
  as.bigz(Tred)
}

####generate one tree under coalescence model####

rhetcoal<-function(n){
  Fmat<-matrix(0,nrow=2*n-2,ncol=2*n-2)
  Fmat[1,1]<-2; Fmat[2*n-2,2*n-2]<-1;Fmat[2*n-3,2*n-3]<-2; Fmat[2*n-2,2*n-3]<-1
  prob<-1
  leaves<-2*n-4 #first two events are necessarily cherries
  vintages<-c(1,2)
  lv<-2
  for (j in 3:(2*n-2)){
    nume1<-leaves*(leaves-1)
    nume2<-nume1+lv*(lv-1)
    if (runif(1)<nume1/nume2){
      prob<-prob*nume1/nume2
      #print("merge two leaves")
      vintages<-c(vintages,j)
      leaves<-leaves-2
      lv<-lv+1
      #this will be a sampling event so diagonal is+1
      diag(Fmat)[2*n-1-j]<-diag(Fmat)[2*n-j]+1
      Fmat[(2*n-j):(2*n-2),2*n-j-1]<-Fmat[(2*n-j):(2*n-2),2*n-j]
    }else{
      #print("merge two vintages")
      who<-sample(vintages,2)
      prob<-prob* 1/ (nume2/2)  ###TYPO HERE!!!
      #print(who)
      lv<-lv-1
      where<-c(which(vintages==who[1]),which(vintages==who[2]))
      vintages<-c(vintages[-where],j)
      diag(Fmat)[2*n-1-j]<-diag(Fmat)[2*n-j]-1
      Fmat[(2*n-j):(2*n-1-max(who)),2*n-j-1]<-Fmat[(2*n-j):(2*n-1-max(who)),2*n-j]-2
      Fmat[(2*n-max(who)):(2*n-1-min(who)),2*n-j-1]<-Fmat[(2*n-max(who)):(2*n-1-min(who)),2*n-j]-1
      if (min(who)>1){
      Fmat[(2*n-min(who)):(2*n-2),2*n-j-1]<-Fmat[(2*n-min(who)):(2*n-2),2*n-j]}
    }
  }
  #print(prob)
  return(Fmat)
}

# Fmat<-rhetcoal(50)
# mytree<-tree_from_hetF(Fmat)
# plot(mytree)





#### obtain a sample ####
coales_sample <- function(sample_size, tree_size,seed=781){
  set.seed(781)
  internal_length_vec = rep(0,sample_size)
  cherry_count_vec = c(0,sample_size)
  total_length_vec = c(0,sample_size)
  
  for (i in 1:sample_size){
    
    current_tree = rhetcoal(tree_size)
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




coales_sample_100_5 = coales_sample(100,5)
coales_sample_100_20 = coales_sample(100,20)
coales_sample_100_50 = coales_sample(100,50)

coales_sample_500_5 = coales_sample(500,5)
coales_sample_500_20 = coales_sample(500,20)
coales_sample_500_50 = coales_sample(500,50)

coales_sample_1000_5 = coales_sample(1000,5)
coales_sample_1000_20 = coales_sample(1000,20)
coales_sample_1000_50 = coales_sample(1000,50)



# set.seed(781)
# sample_size = 100
# tree_size = 5
# test_output = coales_sample(100,n)
# tree_list = test_output$tree
# tr = matrix(tree_list[[24]],nrow = 2*n-2 ,byrow=TRUE)
# plot(tree_from_hetF(tr))
# test_output$cherry;test_output$count;test_output$length
# sum(test_output$cherry * test_output$count)
# 
# data.frame(test_output$count)


