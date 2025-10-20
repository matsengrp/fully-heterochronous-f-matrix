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
  #we will generate all possible diagonals (one per row)
  cc<-2*n-2. #total number of entries in the diagonal
  out<-matrix(1,nrow=2,ncol=cc) #the first 2 obvious diagonals (one per row)
  out[,1]<-2 #always start with 2
  out[,2]<-c(3,1) #it can then be 3 or 1
  budget_m<-c(n-1,n-2) #number of -1 needed (sampling)
  budget_p<-c(n-3,n-2) #number of +1 needed. (coalescence)
  k<-3
  while(k<=cc){
    #out<-cbind(out,rep(1,nrow(out)))
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
        if (budget_m[j]>0 & out[j,k-1]>1){ #only sampling events left
          out[j,k]<-out[j,k-1]-1
          budget_m[j]<-budget_m[j]-1
        }else{
          out[j,k]<-out[j,k-1]+1 #only coalescent events left
          budget_p[j]<-budget_p[j]-1
        }
      }
    }
    k<-k+1
  }
  return(out)
}

diagon_fast <- function(n) {
  cc <- 2 * n - 2
  max_paths <- choose(2 * n - 2, n - 1)  # upper bound
  out <- matrix(NA, nrow = max_paths, ncol = cc)
  budget_m <- rep(NA, max_paths)
  budget_p <- rep(NA, max_paths)
  
  # Initial state
  out[1, 1:2] <- c(2, 3)
  out[2, 1:2] <- c(2, 1)
  budget_m[1:2] <- c(n - 1, n - 2)
  budget_p[1:2] <- c(n - 3, n - 2)
  row_count <- 2
  
  for (k in 3:cc) {
    for (j in 1:row_count) {
      current_val <- out[j, k - 1]
      if (is.na(current_val)) next
      
      if ((budget_m[j] > 0) & (budget_p[j] > 0) & current_val > 1) {
        # +1 branch
        row_count <- row_count + 1
        out[row_count, 1:(k - 1)] <- out[j, 1:(k - 1)]
        out[row_count, k] <- current_val - 1
        budget_m[row_count] <- budget_m[j] - 1
        budget_p[row_count] <- budget_p[j]
        
        # original row gets +1
        out[j, k] <- current_val + 1
        budget_p[j] <- budget_p[j] - 1
      } else if (budget_m[j] > 0 & current_val > 1) {
        out[j, k] <- current_val - 1
        budget_m[j] <- budget_m[j] - 1
      } else {
        out[j, k] <- current_val + 1
        budget_p[j] <- budget_p[j] - 1
      }
    }
  }
  
  return(out[1:row_count, ])
}


diagonal_list<-diagon(7)
diagonal_list2<-diagon_fast(7)
#the cardinality of the diagonal space is is Catalan n-1
catalan<-function(n){
  return(choose(2*n,n)/(n+1))
}
catalan(6)


#we now construct all F-matrices
construct_trees_het_D<-function(n){
  diaglist<-diagon_fast(n)
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
    initial<-matrix(rep(1,k),ncol=k,nrow=1) #Each row of initial will correspond to one D-matrix 
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
      initial<-rbind(initial,initial) ;initial
      initial[1,xnow]<-prev_val-c(1,0)
      initial[2,xnow]<-prev_val-c(0,1)
    }else{
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
  n<-ncol(out$diaglist)
  ncols<-ncol(out$list_out[[1]])
  outM <- matrix(0,nrow=n,ncol=n)
  F.list<-list()
  x<-1
  for (d in 1:nrow(out$diaglist)){
    diag(outM)<-out$diaglist[d,]
    for (i in 1:out$ndiagonal[d]){
      for ( j in 1:ncols){
        coord <- out$F.matrixRel[out$F.matrixRel[,1]==j,]
        if (coord[2]==1){ #if first column
          outM[(coord[3]+1),coord[2]] <- out$list_out[[d]][i,coord[1]]
        }else{
          outM[(coord[3]+1),coord[2]] <- outM[(coord[3]+1),coord[2]-1]+out$list_out[[d]][i,coord[1]]
        }
      }
      F.list[[x]]<-outM
      x<-x+1
    }
  }
  return(F.list)
}

F.list<-Fmat_from_construct_het(out)

tree_from_hetF <- function(Fmat){
  #assume time between events have length 1
  n <- (ncol(Fmat)+2)/2 #max diagonal value gives the number of tips
  #if(is.null(coal.times)) coal.times <- seq(1, n-1)
  #first column parent, second column offspring
  edge=matrix(rep(0,4*n-4),ncol=2)
  edge[,2]= (1:(2*n-2))
  #vintages = c()
  #times=c(rep(0,n),coal.times)
  N<-ncol(Fmat)
  F_difference = c(0,diff(diag(Fmat)))
  for (j in 1:N){
    new_node = 2*n-j
    setleaves<-seq(2*n-j-1,2)[diff(Fmat[j:N,j])==-1]
    for (i in setleaves){
      if (edge[which(edge[,2]==i),1]==0){ #this is eq. to looking at D matr. first -1
        edge[which(edge[,2]==i),1]<-new_node 
      }
    }
    
  }
  #who does not have a parent?
  pa<-max(seq(2*n-1,2)[Fmat[N,]==1])
  edge[which(edge[,2]==1),1]<-pa
  edgelength<-edge[,1]-edge[,2]
  new.edge<-edge
  pa<-sort(unique(edge[,1]))
  ch<-seq(1,2*n-1)
  for (i in pa){
    ch<-ch[-which(ch==i)]
  }
  for (j in 1:(n-1)){
    new.edge[,1][edge[,1]==pa[j]]<-2*n-j
    new.edge[,2][edge[,2]==pa[j]]<-2*n-j
  }
  for (j in 1:n){
    new.edge[,2][edge[,2]==ch[j]]<-j
  }
  final_tree=rcoal(n)
  final_tree$edge=new.edge
  final_tree$edge.length=edgelength
  #final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  #final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
  
  #tree_from_F(rbind(rep(0,n),cbind(rep(0,n-1),Fmat)), coal.times)
}

ex3=F.list[[1]]
ex3=tree_from_hetF(F.list[[1]])
plot(ex3)

## Basic Tree Plot
ex3$tip.label<-c("6","5","4","3")
ex3$node.label<-c(0,1,2)
p3 <- ggtree(ex3, right = TRUE) +
  geom_tiplab(align = FALSE, linesize = 0.5, size = 5) +
  layout_dendrogram() +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5, vjust = -0.5, size = 5)
#+ 
# geom_hline(data = internal_nodes,
#                                                                                                  aes(yintercept = y),
#                                                                                                  linetype = "dashed",
#                                                                                                  color = "grey50")

plot(p3)


p <- ggtree(ex3, right = TRUE) + layout_dendrogram()
tree_data <- p$data
all_times <- sort(unique(tree_data$x))
n<-ex3$Nnode+1
# Create subscript labels as plotmath expressions
labels <- paste0("u[", seq(0,2*n-2), "]")

# Positioning
y_max <- max(tree_data$y) + 0.5

label_df <- data.frame(
  x = all_times,
  y = y_max,
  label = labels
)

# Final plot with subscripted labels
p3 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  geom_vline(xintercept = all_times, linetype = "dashed", color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), parse = TRUE, hjust = -0.1, size = 5)

print(p3)





y_min <- min(tree_data$y)

# Create segment data frame for vertical dashed lines cut at y_max (label height)
segment_df <- data.frame(
  x = all_times,
  y_start = y_min,
  y_end = y_max - 0.1  # slightly below label so line stops just before label
)

prob_df <-data.frame(
  x= all_times, 
  y= rep(0.5,7),
  label = c("","2/2","1/3","2/4","1/3","2/2","")
)

# Update plot:
p3 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  # Replace geom_vline with geom_segment for cut lines
  geom_segment(data = segment_df,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               linetype = "dashed",
               color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), 
          parse = TRUE, hjust = -0.1, size = 5)+
  geom_text(data=prob_df, aes(x = x, y = y, label = label),
            parse=TRUE, size=5,color="gray60")

print(p3)


