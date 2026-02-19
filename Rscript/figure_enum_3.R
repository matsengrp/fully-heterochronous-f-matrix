library(ape)
library(dendextend)
library(ggtree)
library(ggplot2)

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




ex0=matrix(data=c(2,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1),nrow=4,byrow=TRUE)

paper_tree_plot <- function(fmat, tip_label , node_label ){

ex0=tree_from_hetF(fmat)
plot(ex0)

## Basic Tree Plot
ex0$tip.label<-tip_label
ex0$node.label<-node_label
p0 <- ggtree(ex0, right = TRUE) +
  geom_tiplab(align = FALSE, linesize = 0.5, size = 5) +
  layout_dendrogram() +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5, vjust = -0.5, size = 5)
#+ 
# geom_hline(data = internal_nodes,
#                                                                                                  aes(yintercept = y),
#                                                                                                  linetype = "dashed",
#                                                                                                  color = "grey50")

plot(p0)


p <- ggtree(ex0, right = TRUE) + layout_dendrogram()
tree_data <- p$data
all_times <- sort(unique(tree_data$x))
n<-ex0$Nnode+1
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
p0 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  geom_vline(xintercept = all_times, linetype = "dashed", color = "gray60") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), parse = TRUE, hjust = -0.1, size = 5)

print(p0)





y_min <- min(tree_data$y)

# Create segment data frame for vertical dashed lines cut at y_max (label height)
segment_df <- data.frame(
  x = all_times,
  y_start = y_min,
  y_end = y_max - 0.1  # slightly below label so line stops just before label
)

# Update plot:
p0 <- p +
  geom_tiplab(angle=0,size = 5,hjust=-0.5) +
  geom_text2(aes(subset = !isTip, label = label), hjust=-0.5,vjust = -0.5, size = 5) +
  # Replace geom_vline with geom_segment for cut lines
  geom_segment(data = segment_df,
               aes(x = x, xend = x, y = y_start, yend = y_end),
               linetype = "dashed",
               color = "gray60") 
geom_text(data = label_df, aes(x = x, y = y, label = label), 
          parse = TRUE, hjust = -0.1, size = 5)

return(p0)
}

ex0=matrix(data=c(2,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1),nrow=4,byrow=TRUE)
tip_label<-c("4","3","1") #from bottom up
node_label<-c(0,2) #from root down
p0= paper_tree_plot(fmat=ex0, tip_label,node_label)

ex1=matrix(data=c(2,0,0,0,1,3,0,0,0,2,2,0,0,1,1,1),nrow=4,byrow=TRUE)
tip_label<-c("4","3","2") #from bottom up
node_label<-c(0,1) #from root down
p1= paper_tree_plot(fmat=ex1, tip_label,node_label)

ex2=matrix(data=c(2,0,0,0,1,3,0,0,1,2,2,0,0,1,1,1),nrow=4,byrow=TRUE)
tip_label<-c("4","3","2") #from bottom up
node_label<-c(0,1) #from root down
p2= paper_tree_plot(fmat=ex2, tip_label,node_label)

ex3=matrix(data=c(2,0,0,0,1,3,0,0,1,2,2,0,1,1,1,1),nrow=4,byrow=TRUE)
tip_label<-c("4","3","2") #from bottom up
node_label<-c(0,1) #from root down
p3= paper_tree_plot(fmat=ex3, tip_label,node_label)

p0 <- p0 + ggtitle(expression(F^0)) 
p1 <- p1 + ggtitle(expression(F^1)) 
p2 <- p2 + ggtitle(expression(F^2)) 
p3 <- p3 + ggtitle(expression(F^3)) 

library(patchwork)
(p0 + p1 + p2 + p3) + plot_layout(ncol = 4)&
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 12)),
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5)
  )

