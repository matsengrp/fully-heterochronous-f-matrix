#### Functions used for isochronous F-matrices (implemented in fmatrix)


#Generates a table with all isochronous trees with n leaves (the format is a bit weird)
construct_trees<-function(n){
  ##only works for n>4
  ##This function generates a form of F matrices for all trees with n tips, It takes long for n~10
  
  k<-(n-2)*(n-1)/2
  k2<-(n-3)*(n-2)/2
  initial<-rep(1,k)
  first.column<-seq(1,k2)*(seq(1,k2)-1)/2 + 1
  first.column<-first.column[first.column<k]
  diagonal<-(seq(1,k2+1)*(seq(1,k2+1)-1)/2)[-1]
  diagonal<-diagonal[diagonal<=k]
  initial[diagonal]<-seq(1,length(first.column))
  initial<-rbind(initial,initial)
  initial[2,2]<-0
  #x_index, column, row
  F.matrixRel<-cbind(first.column,rep(1,length(first.column)),seq(1,length(first.column)))
  j<-2
  F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1)))
  for (j in 3:length(first.column)){ #row in F matrix
    F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1)))
    #column in F matrix
    for (cc in 1:(j-1)){
      xval<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==j,1]
      nrows.to<-nrow(initial)
      for (i in 1:nrows.to){
        if (cc==1){ #it is a first column element
          prevxval<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==(j-1),1]
          cond1<-max(c(0,initial[i,prevxval]-1))
          cond2<-initial[i,prevxval]
          newval<-seq(cond1,cond2)
          initial[i,xval]<-newval[1]
          if (length(newval)>1){
            for (z in 2:length(newval)){
              temp<-initial[i,]
              temp[xval]<-newval[z]
              initial<-rbind(initial,temp)
            }
          }
        }else{
          prevxval<-F.matrixRel[F.matrixRel[,2]==(cc-1) & F.matrixRel[,3]==j,1]
          prevxval2<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==(j-1),1]
          prevxval3<-F.matrixRel[F.matrixRel[,2]==(cc-1) & F.matrixRel[,3]==(j-1),1]
          cond1<-max(c(0, initial[i,prevxval], initial[i,prevxval2]-1, initial[i,prevxval]+initial[i,prevxval2]-initial[i,prevxval3]-1))
          cond2<-min(initial[i,prevxval]+initial[i,prevxval2]-initial[i,prevxval3],initial[i,prevxval2])
          newval<-seq(cond1,cond2)
          initial[i,xval]<-newval[1]
          if (length(newval)>1){
            for (z in 2:length(newval)){
              temp<-initial[i,]
              temp[xval]<-newval[z]
              initial<-rbind(initial,temp)
            }
          }
        }
        
      }
      
    }
  }
  return(list(F.matrixRel=F.matrixRel,res=initial))
}

##Generates a table with all isochronous trees with n leaves (the format is a bit weird)
Fmat_from_construct <- function(n, output, relpos){
  #This function converts the output of the construct_trees function to
  
  # Expected syntax:
  # F.list<-list()
  # for (j in 1:nrow(res$res)){
  #   F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
  # }
  out <- matrix(0,nrow=n-1,ncol=n-1)
  diag(out) <- seq(2,n)
  for (i in 1:length(output)){
    coord <- relpos[relpos[,1]==i,]
    out[(coord[3]+1),coord[2]] <- output[i]
  }
  return(out)
}

tt<-construct_trees(5)
F.list<-list()
for (j in 1:nrow(tt$res)){
  F.list[[j]]<- Fmat_from_construct(5, tt$res[j,], tt$F.matrixRel)
}
length(F.list)


##Computes the Euler zigzag number
zigzag<-function(n){
  factor<-rep(1,n)
  zig<-rep(1,n)
  for (j in 2:n){
    factor[j]<-factor[j-1]*(j-1)
  }
  for ( j in 3:n){
    sum=0
    for (k in 1:(j-1)){
      sum=sum+(factor[j-1]/(factor[j-k]*factor[k]))*zig[k]*zig[j-k]
    }
    zig[j]<-sum/2
  }
  return(zig)
}

##The idea is to figure out the diagonals first 
##June 16, 2025

diagon<-function(n){
  #n is the number of tips
  m<-factorial(2*(n-1))/(factorial(n)*factorial(n-1))
  print(m)
  cc<-2*n-2
  out<-matrix(1,2,2)
  out[,1]<-2
  out[,2]<-c(3,1)
  budget_m<-c(n-1,n-2) #number of -1 needed
  budget_p<-c(n-3,n-2) #number of +1 needed
  k<-3
  while(k<=cc){
    out<-cbind(out,rep(1,nrow(out)))
    trow<-nrow(out)
    for ( j in 1:trow){
      if (((budget_m[j]>0)+(budget_p[j]>0)==2) & (out[j,k-1]>1)){
        out<-rbind(out,out[j,])
        out[j,k]<-out[j,k-1]+1
        out[nrow(out),k]<-out[nrow(out),k-1]-1
        budget_m[nrow(out)]<-budget_m[j]-1
        budget_p[nrow(out)]<-budget_p[j]
        budget_p[j]<-budget_p[j]-1
      }else{
        if (budget_m[j]>0 & out[j,k-1]>1){
          out[j,k]<-out[j,k-1]-1
          budget_m[j]<-budget_m[j]-1
        }else{
          out[j,k]<-out[j,k-1]+1
          budget_p[j]<-budget_p[j]-1
        }
      }
    }
    k<-k+1
  }
return(out)
}

diagonal_list<-diagon(4)

catalan<-function(n){
  return(choose(2*n,n)/(n+1))
}
catalan(3)

construct_trees_het_D<-function(n){
  diaglist<-diagon(n)
  ndiagonal<-rep(0,nrow(diaglist))
  F.list<-list()
  n1<-n
  n<-2*n-1
  list_out<-list()
  for (d in 1:nrow(diaglist)){
    diag<-diaglist[d,]
    ##This function generates all D matrices
    ##We know the value of the diagonal, so we just figure out the k other values, where k is
    k<-(n-2)*(n-1)/2
    k2<-(n-3)*(n-2)/2
    initial<-matrix(rep(1,k),ncol=k,nrow=1) #Each row of initial will have one D-matrix 
    #indices are by row in the triangular array
    first.column<-seq(1,k2)*(seq(1,k2)-1)/2 + 1 #indices for the values of the first column
    first.column<-first.column[first.column<k] #indices for first column
    diagonal<-(seq(1,k2+1)*(seq(1,k2+1)-1)/2)[-1]
    diagonal<-diagonal[diagonal<=k] #indices for lower diagonal
    #initial[diagonal]<-diag[-length(diag)]-1 #setting the values of the diagonal 
    #initial<-rbind(initial,initial) #repeat and set the 2,2 to 0
    #initial[2,2]<-0
    #x_index, column, row
    type=c(2,diff(diag)+1)
    F.matrixRel<-cbind(first.column,rep(1,length(first.column)),seq(1,length(first.column)))
    j<-2
    F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1))) #indices of the first column and diagonal
    xprev<-F.matrixRel[,1][F.matrixRel[,3]==j-1] #first row
    prev_val<-initial[xprev]
    xnow<-F.matrixRel[,1][F.matrixRel[,3]==j]
    if (type[j]==2){
      prev_val<-c(prev_val,2)
      initial<-rbind(initial,initial) 
      initial[1,xnow]<-prev_val-c(1,0)
      initial[2,xnow]<-prev_val-c(0,1)}else{
        prev_val<-c(prev_val,0)
        initial[xnow]<-c(0,0)
      }
    for (j in 3:length(first.column)){ #row in F matrix
      F.matrixRel<-rbind(F.matrixRel,cbind(first.column[j]+seq(1,j-1),2:j,rep(j,j-1)))
      xprev<-F.matrixRel[,1][F.matrixRel[,3]==(j-1)]
      xnow<-F.matrixRel[,1][F.matrixRel[,3]==j]
      #column in F matrix
      #for (cc in 1:(j-1)){
      # xval<-F.matrixRel[F.matrixRel[,2]==cc & F.matrixRel[,3]==j,1] # the index of coordinate in cc row and j column
      nrows.to<-nrow(initial)
      for (i in 1:nrows.to){
        prev_val<-c(initial[i,xprev],type[j])
        lprev<-length(prev_val)
        where.to<-seq(1:lprev)[prev_val>0]
        lnewval<-length(where.to)
        subs<-rep(0,lprev)
        subs[where.to[1]]<-1
        initial[i,xnow]<-prev_val-subs
        x<-2
        if (lnewval>1){
          for (z in 2:lnewval){
            subs<-rep(0,lprev)
            subs[where.to[x]]<-1
            temp<-initial[i,]
            temp[xnow]<-prev_val-subs
            initial<-rbind(initial,temp)
            x<-x+1
          }
        }
        
        
      }
      
      #}
    }
    ##Paste diagonal
   ndiagonal[d]<-nrow(initial)
   list_out[[d]]<-initial
    
  }
  return(list(list_out=list_out,ndiagonal=ndiagonal,diaglist=diaglist,F.matrixRel=F.matrixRel))

}


n<-4
out<-construct_trees_het_D(n)
sum(out$ndiagonal)
zigzag(2*n)[2*n]/2^(n-1)
##Yay, it matches!

diagval<-diagonal_list[[1]]

Fmat_from_construct_het <- function(out){
  n<-max(out$F.matrixRel[,2])
  ncols<-ncol(out$list_out[[1]])
  outM <- matrix(0,nrow=n-1,ncol=n-1)
  F.list<-list()
  for (d in 1:nrow(out$diaglist)){
    diag(outM)<-out$diaglist[d,]
    for (i in 1:out$ndiagonal[d]){
      for ( j in 1:ncols){
        coord <- out$F.matrixRel[out$F.matrixRel[,1]==j,]
        outM[(coord[3]+1),coord[2]] <- output[i]
      }
    }
  }
  
  #This function converts the output of the construct_trees function to
  
  # Expected syntax:
  # F.list<-list()
  # for (j in 1:nrow(res$res)){
  #   F.list[[j]]<- Fmat_from_construct(n, res$res[j,], res$F.matrixRel)
  # }
 
  diag(out) <- diagval
  for (i in 1:length(output)){
    coord <- relpos[relpos[,1]==i,]
    out[(coord[3]+1),coord[2]] <- output[i]
  }
  return(out)
}

Fmat_from_construct_het(7, diagval, res$res[j,], res$F.matrixRel)
#This should give me 20 F-matrices, see for example line 12 is invalid


F.list2<-list()
for (j in 1:nrow(res$res)){
  F.list2[[j]]<- Fmat_from_construct_het(7, diagval, res$res[j,], res$F.matrixRel)
}
length(F.list2)

