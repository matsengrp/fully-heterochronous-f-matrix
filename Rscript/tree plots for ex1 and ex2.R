library(ape)
library(dendextend)
library(ggtree)
library(ggplot2)

#### Example 1: heterochronous case ####
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

ex1=matrix(data=c(2,0,0,0,1,3,0,0,1,2,2,0,0,1,1,1),nrow=4,byrow=TRUE)
ex1=tree_from_hetF(ex1)
plot(ex1)

## Basic Tree Plot
ex1$tip.label<-c("4","3","2")
ex1$node.label<-c(0,1)
p1 <- ggtree(ex1, right = TRUE) +
  geom_tiplab(align = FALSE, linesize = 0.5, size = 5) +
  layout_dendrogram() +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5, vjust = -0.5, size = 5)
#+ 
# geom_hline(data = internal_nodes,
#                                                                                                  aes(yintercept = y),
#                                                                                                  linetype = "dashed",
#                                                                                                  color = "grey50")

plot(p1)


p <- ggtree(ex1, right = TRUE) + layout_dendrogram()
tree_data <- p$data
all_times <- sort(unique(tree_data$x))
n<-ex1$Nnode+1
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
p1 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  geom_vline(xintercept = all_times, linetype = "dashed", color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), parse = TRUE, hjust = -0.1, size = 5)

print(p1)





y_min <- min(tree_data$y)

# Create segment data frame for vertical dashed lines cut at y_max (label height)
segment_df <- data.frame(
  x = all_times,
  y_start = y_min,
  y_end = y_max - 0.1  # slightly below label so line stops just before label
)

# Update plot:
p1 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  # Replace geom_vline with geom_segment for cut lines
  geom_segment(data = segment_df,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               linetype = "dashed",
               color = "gray60") 
  geom_text(data = label_df, aes(x = x, y = y, label = label), 
             parse = TRUE, hjust = -0.1, size = 5)

print(p1)





#### Example 2: isochronous tree ####
library(ape)
library(ggtree)

ex2=rcoal(5,tip.label=NULL)

ex2 <- structure(
  list(
    edge = matrix(c(
      8,1,8,2,7,3,9,4,9,5,7,8,6,7,6,9
    ), ncol = 2, byrow = TRUE),
    edge.length = c(2,2,3,1,1,1,1,3),
    tip.label   = c("t5","t2","t4","t3","t1"),
    Nnode       = 4      
  ),
  class = "phylo"
)

plot(ex2)

ex2$tip.label<- c(" ", " "," "," "," ")
ex2$node.label<-c(0,1,2,3)
p2  <- ggtree(ex2, right = TRUE) +
  geom_tiplab(align = FALSE, linesize = 0.5, size = 5) +
  layout_dendrogram() +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5, vjust = -0.5, size = 5)
plot(p2 )


p <- ggtree(ex2, right = TRUE) + layout_dendrogram()
tree_data <- p$data
all_times <- sort(unique(tree_data$x))
n<-ex2$Nnode+1
# Create subscript labels as plotmath expressions
labels <- paste0("u[", seq(0,n-1), "]")

# Positioning
y_max <- max(tree_data$y) + 0.5

label_df <- data.frame(
  x = all_times,
  y = rep(y_max,5),
  label = labels
)

# Final plot with subscripted labels
p2  <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  geom_vline(xintercept = all_times, linetype = "dashed", color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), parse = TRUE, hjust = -0.1, size = 5)

print(p2 )





y_min <- min(tree_data$y)

# Create segment data frame for vertical dashed lines cut at y_max (label height)
segment_df <- data.frame(
  x = all_times,
  y_start = y_min,
  y_end = y_max - 0.1  # slightly below label so line stops just before label
)

# Update plot:
p2  <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  # Replace geom_vline with geom_segment for cut lines
  geom_segment(data = segment_df,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               linetype = "dashed",
               color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), 
            parse = TRUE, hjust = -0.1, size = 5)

print(p2 )






print(p1)
print(p2)
