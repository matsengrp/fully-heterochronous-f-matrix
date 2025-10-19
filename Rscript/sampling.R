##Code for simulating fully heterochronous trees with n tips
library("ape")
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
      prob<-prob*(1-nume1/nume2)
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
  print(prob)
  return(Fmat)
}

Fmat<-rhetcoal(6)
mytree<-tree_from_hetF(Fmat)
plot(mytree)
