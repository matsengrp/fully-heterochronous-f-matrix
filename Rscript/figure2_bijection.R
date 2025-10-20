#plot(dend, yaxt = "n") #no side-bar of height
library(ape)
library(dendextend)
library(ggtree)
library(ggplot2)

set.seed(1)
hc <- hclust(dist(1:6), method = "complete")
hc$height = c(3, 1, 2, 4, 5)
dend <- as.dendrogram(hc)

par(mfrow = c(1, 2), oma = c(0,0,0,0))  # no outer margins

## left panel: small right margin
par(mar = c(0.5, 0.5, 0.5, 0.05))

d1= dend %>% set("labels_col", "white") %>% 
  
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>% set("leaves_col","gray") %>% 
  set("by_labels_branches_col", value = c(1), TF_values = c("white",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(2), TF_values = c("white",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(3), TF_values = c("white",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(4), TF_values = c("white",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(5), TF_values = c("white",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(6), TF_values = c("white",Inf),type="all") %>%
  
  
  set("by_labels_branches_lwd", value = c(1), TF_values = c(0,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(2), TF_values = c(0,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(3), TF_values = c(0,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(4), TF_values = c(0,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(5), TF_values = c(0,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(6), TF_values = c(0,Inf),type="all") %>%
  
  set("nodes_pch", 19) %>% set("nodes_cex", 0.5) %>% set("nodes_col","black") %>%
  set("leaves_col","white") 


plot(d1, yaxt = "n")

## Get node coordinates
xy <- dendextend::get_nodes_xy(d1)      
ints <- xy[xy[,2] != 0,]  
points(ints[,1], ints[,2],pch = 19, cex = 0.4)
n= dim(ints)[1]
ints[ints[,2]!= n,1] = ints[ints[,2]!= n,1]+0.2

# labels:
lab_int <- c(0,3,1,2,4) 
text(ints[,1], ints[,2], labels = lab_int, pos = 3, cex = 0.8, xpd = NA) $plot

## right panel: small left margin
par(mar = c(0.5, 0.05, 0.5, 0.5))

d2 = dend %>% set("labels_col", "white") %>% 
  set("nodes_pch", 19) %>% set("nodes_cex", 0.5) %>% set("nodes_col","black") %>% 
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>% set("leaves_col","gray") %>% 
  set("by_labels_branches_col", value = c(1), TF_values = c("gray",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(2), TF_values = c("gray",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(3), TF_values = c("gray",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(4), TF_values = c("gray",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(5), TF_values = c("gray",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(6), TF_values = c("gray",Inf),type="all") %>%
  
  set("by_labels_branches_lwd", value = c(1), TF_values = c(2,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(2), TF_values = c(2,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(3), TF_values = c(2,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(4), TF_values = c(2,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(5), TF_values = c(2,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(6), TF_values = c(2,Inf),type="all") %>%
  
  set("by_labels_branches_lty", value = c(1), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(2), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(3), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(4), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(5), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(6), TF_values = c(3,Inf),type="all") 

plot(d2, yaxt = "n")

## Get node coordinates
xy <- dendextend::get_nodes_xy(d2)      
ints <- xy[xy[,2] != 0,]  
points(ints[,1], ints[,2],pch = 19, cex = 0.4)
n= dim(ints)[1]
ints[ints[,2]!= n,1] = ints[ints[,2]!= n,1]+0.2

# labels:
lab_int <- c(0,3,1,2,4) 
text(ints[,1], ints[,2], labels = lab_int, pos = 3, cex = 0.8, xpd = NA) $plot



