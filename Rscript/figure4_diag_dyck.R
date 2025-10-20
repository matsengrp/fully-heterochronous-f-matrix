library(ggplot2) 
library(gridExtra)
library(cowplot)

dt_dyck <- data.frame(
  x = c(0,1,1,2,3,3,4,4,4,5,5),
  y=c(0,0,1,1,1,2,2,3,4,4,5)
)

dt_diag <- data.frame(
  x=seq(0,9,1),
  y=c(2,1,2,3,2,3,2,1,2,1)
)

diag_grid_h <- data.frame(
  x=c(0,0,0,0),
  xend=c(9,9,9,9),
  y=c(0,1,2,3),
  yend=c(0,1,2,3)
)
diag_grid_v <- data.frame(
  x=seq(0,9,1),
  xend=seq(0,9,1),
  y=rep(0,10),
  yend=rep(3.5,10)
)

dyck_grid_v <- data.frame(
  x=seq(0,5,1),
  y=rep(0,6),
  xend=seq(0,5,1),
  yend=rep(5,6)
)

dyck_grid_h <- data.frame(
  x=rep(0,6),
  y=seq(0,5,1),
  xend=rep(5,6),
  yend=seq(0,5,1)
)

diag_lab <- data.frame(
  x= c(diag_grid_v$x,diag_grid_h$x-0.3),
  y= c(diag_grid_v$y+0.2,diag_grid_h$y+0.5),
  label = c(seq(0,9),seq(0,3))
)

diag_text <- data.frame(
  x=c(4.5,-1),
  y=c(-0.2,2.5),
  label=c("j","F[list(j,j)]")
)


par(mfrow = c(1, 2), oma = c(0,0,0,0)) 

p1 = ggplot()+
  geom_path(data=dt_diag,aes(x=x,y=y+0.5))+
  coord_cartesian(xlim=c(-1.5,10),ylim=c(-1,5),expand=FALSE)+
  geom_segment(data=diag_grid_h,
               aes(x=x,y=y+0.5,xend=xend,yend=yend+0.5),
               linetype='dotted',color="grey60")+
  geom_segment(data=diag_grid_v,
               aes(x=x,y=y+0.5,xend=xend,yend=yend+0.5),
               linetype='dotted',color="grey60")+
  geom_text(data=diag_lab,aes(x=x,y=y,label=label),color='grey40')+
  geom_text(data=diag_text,aes(x=x,y=y,label=label),parse=TRUE)+
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

p1

dyck_lab <- data.frame(
  x=c(seq(0,5,1),rep(0,6)-0.2),
  y=c(rep(0,6)-0.2,seq(0,5,1)),
  label=c(seq(0,5,1),seq(0,5,1))
)


p2 = ggplot()+
  geom_path(data=dt_dyck,aes(x=x,y=y))+
  coord_cartesian(xlim=c(-1,5.5),ylim=c(-1,5.5),expand=FALSE)+ 
  geom_segment(data=dyck_grid_h,aes(x=x,y=y,xend=xend,yend=yend),linetype='dotted',color="grey60")+
  geom_segment(data=dyck_grid_v,aes(x=x,y=y,xend=xend,yend=yend),linetype='dotted',color="grey60")+
  geom_text(data=dyck_lab,aes(x=x,y=y,label=label),color='grey40')+ 
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

p2 
  
 
#grid.arrange(p1, p2, ncol = 2)
 plot_grid(p1,p2,align="h",ncol=2,rel_widths=c(3/5,2/5),
           labels =c("F-matrix diagonal", "Dyck path"),
           label_size = 16, label_fontface = "bold",
           label_x = c(0.5, 0.5),   # center
           hjust   = 0.5)
