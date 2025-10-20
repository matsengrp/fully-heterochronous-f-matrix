# my plot(hc) looks like figure 1, but my plot(dend) looks like figure 2, is there a way to make (5,6) stay where it is? 
library(ape)
library(dendextend)
library(ggtree)
library(ggplot2)

library(dendextend)
library(latex2exp)
## Steps / heights wanted:
## 1: (1,2) @ 1
## 2: (3,4) @ 2
## 3: ((1,2),(3,4)) @ 3
## 4: (5,6) @ 5
## 5: root = ((5,6), ((1,2),(3,4))) @ 6

merge <- matrix(c(
  -1, -2,   # row 1 -> cluster 1
  -3, -4,   # row 2 -> cluster 2
  1,  2,   # row 3 -> cluster 3 (merge clusters from rows 1 & 2)
  -5, -6,   # row 4 -> cluster 4
  4,  3    # row 5 -> root (merge cluster 4 with cluster 3)
), ncol = 2, byrow = TRUE)

height <- c(1, 2, 3, 4, 5)     # nondecreasing merge heights
ord    <- c( 1, 2, 3, 4, 5, 6)  # desired leaf order L->R

hc <- structure(
  list(merge = merge, height = height, order = ord,
       method = "user", dist.method = "none"),
  class = "hclust"
)

plot(hc)                 # base plot
# or
dend= as.dendrogram(hc)
dend = dendextend::rotate(dend, order = as.character(c(1,2,3,4,5,6)))
plot(dend)  # dendrogram plot

d = dend %>% set("labels_col", "white") %>% 
  set("nodes_pch", 19) %>% set("nodes_cex", 0.5) %>% set("nodes_col","black") %>% 
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) %>% set("leaves_col","gray") %>% 
  set("by_labels_branches_col", value = c(1), TF_values = c("gray40",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(2), TF_values = c("gray40",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(3), TF_values = c("gray40",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(4), TF_values = c("gray40",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(5), TF_values = c("gray40",Inf),type="all") %>%
  set("by_labels_branches_col", value = c(6), TF_values = c("gray40",Inf),type="all") %>%
  
  set("by_labels_branches_lwd", value = c(1), TF_values = c(1,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(2), TF_values = c(1,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(3), TF_values = c(1,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(4), TF_values = c(1,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(5), TF_values = c(1,Inf),type="all") %>%
  set("by_labels_branches_lwd", value = c(6), TF_values = c(1,Inf),type="all") %>%
  
  set("by_labels_branches_lty", value = c(1), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(2), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(3), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(4), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(5), TF_values = c(3,Inf),type="all") %>%
  set("by_labels_branches_lty", value = c(6), TF_values = c(3,Inf),type="all") 

#p  = plot(d, yaxt = "n")
#line_y <- sort(unique(xy[xy[,2] > 0, 2]))  # internal-node heights only

#abline(h = line_y, lty = 3, col = "grey70")  # dotted guidelines
#axis(2, at = line_y)  

gd = as.ggdend(d) 

xy <- dendextend::get_nodes_xy(d)

line_y <- sort(unique(xy[,2] ))
lines_df <- data.frame(y = line_y)

labels_df = data.frame(x = xy[xy[,2] != 0,1]+0.1  ,
                       y = xy[xy[,2] != 0,2]+0.2, 
                       label = c(0,2,4,3,1)   )
points_df = data.frame(x = xy[xy[,2] != 0,1]  ,
                       y = xy[xy[,2] != 0,2]   )

states_tex <- c(
  "A[5] == paste('(', 6, ',', paste('{', '   ', '}', sep=''), ')')",
  "A[4] == paste('(', 4, ',', paste('{', 4, '}', sep=''), ')')",
  "A[3] == paste('(', 2, ',', paste('{', 3, ',', 4, '}', sep=''), ')')",
  "A[2] == paste('(', 2, ',', paste('{', 2, '}', sep=''), ')')",
  "A[1] == paste('(', 0, ',', paste('{', 1, ',', 2, '}', sep=''), ')')",
  "A[0] == paste('(', 0, ',', paste('{', '   ', '}', sep=''), ')')",  # ∅ (empty set)
  "'States'"
)

states_df = data.frame(x=rep(7.5,7), 
                       y= c(0,1,2,3,4,5,5.7)+0.3,
                       state = states_tex)
#states_df$state <- latex2exp::TeX(states_tex)

prob_df = data.frame (x=c(rep(12,6),11,11),
                      y= c(line_y,6.15,5.73)+0.3,
                      prob = c(1,1,1,"1/2",1,1,"Transition","probability") )

empty_df <- data.frame(
  x = c(9.7, 9.7),
  y = c(0.3, 5.3),
  lab = c("symbol('\\306')", "symbol('\\306')" )
)

ggplot(gd) +
  geom_segment(data = lines_df,
               aes(y = y, yend = y),
               x = 0, xend = 7,               # start & end in x
               linetype = "dotted",
               linewidth = 0.7,
               color = "grey", 
               inherit.aes = FALSE
  ) +
  geom_text(data=labels_df,aes(x = x, y = y, label = label), parse = TRUE, hjust = -0.1, size = 5,
            inherit.aes = FALSE)+
  geom_point(data = points_df, aes(x=x, y=y),size=2, inherit.aes = FALSE)+
  geom_text(data=states_df,aes(x=x, y=y,label=state),parse=TRUE,hjust = -0.1, size = 7, inherit.aes = FALSE )+ 
  geom_text(data = prob_df,aes(x=x,y=y,label=prob), hjust = -0.1, size = 7, inherit.aes = FALSE)+
  geom_text(data = empty_df, aes(x, y, label = lab),
            parse = TRUE, size = 7, inherit.aes = FALSE)+
  coord_cartesian(xlim = c(0, 13), ylim = c(0, 7.5), clip = "off") # +
scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  theme(plot.margin = margin(10, 80, 10, 10))  # big right margin

#add internal node labels
