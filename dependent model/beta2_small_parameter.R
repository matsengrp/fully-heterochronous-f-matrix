source("beta2.R")
#beta2_5_1.5_1.1=beta2_sample(sample_size = 1000,tree_size = 5,alpha=1.5,beta=1.1)
beta2_20_1.5_1.1=beta2_sample(sample_size = 1000,tree_size = 20,alpha=1.5,beta=1.1)
#beta2_50_1.5_1.1=beta2_sample(sample_size = 1000,tree_size = 50,alpha=1.5,beta=1.1)


#beta2_5_1.1_1.5=beta2_sample(sample_size = 1000,tree_size = 5,beta=1.1,alpha=1.5)
beta2_20_1.1_1.5=beta2_sample(sample_size = 1000,tree_size = 20,alpha=1.1,beta=1.5)
#beta2_50_1.1_1.5=beta2_sample(sample_size = 1000,tree_size = 50,beta=1.1,alpha=1.5)


#beta2_5_0.3_0.7=beta2_sample(sample_size = 1000,tree_size = 5,alpha =0.3, beta=0.7)
beta2_20_0.3_0.7=beta2_sample(sample_size = 1000,tree_size = 20,alpha=0.3, beta=0.7)
#beta2_50_0.3_0.7=beta2_sample(sample_size = 1000,tree_size = 50,alpha=0.3,beta=0.7)

beta2_20_0.7_0.3=beta2_sample(sample_size = 1000,tree_size = 20,alpha=0.7, beta=0.3)

beta2_20_0.5_0.5=beta2_sample(sample_size = 1000,tree_size = 20,alpha=0.5, beta=0.5)










# alpha=1.1 beta=1.5
#lightblue in colored version
p1 <- ggplot(data.frame(int_length = beta2_20_1.1_1.5$int_length),
             aes(x = int_length)) +
  #geom_histogram(bins = 6, color = "gray80", fill = "gray80") +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(y="count")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7))+
  xlim(20,140)



p2 <- ggplot(data.frame(total_length = beta2_20_1.1_1.5$total_length),
             aes(x = total_length)) +
  #geom_histogram(bins = 14, color = "gray80", fill = "gray80") +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "total tree length", title=expression(bold("Beta2 (" ~alpha == 1.1~","~ beta==1.5 ~ ")")))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.ticks = element_blank() )+
  xlim(56,400)



p3 <- ggplot(data.frame(cherry = beta2_20_1.1_1.5$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank() )+
  xlim(0,11)

# alpha =1.5,beta=1.1
p4 <- ggplot(data.frame(int_length = beta2_20_1.5_1.1$int_length),
             aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "internal tree length",y="count")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7) )+
  xlim(20,140)

p5 <- ggplot(data.frame(total_length = beta2_20_1.5_1.1$total_length),
             aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "total tree length", title=expression(bold("Beta2 (" ~alpha == 1.5~","~ beta==1.1 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.ticks = element_blank() )+
  xlim(56,400)

p6 <- ggplot(data.frame(cherry = beta2_20_1.5_1.1$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,11)


# alpha=0.3 beta=0.7

p7 <- ggplot(data.frame(int_length = beta2_20_0.3_0.7$int_length),
             aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Internal tree length",y="count")+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_blank())+
  xlim(20,140)



p8 <- ggplot(data.frame(total_length = beta2_20_0.3_0.7$total_length),
             aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Total tree length", title=expression(bold("Beta2 (" ~alpha == 0.3~","~ beta==0.7 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(56,400)



p9 <- ggplot(data.frame(cherry = beta2_20_0.3_0.7$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") + 
  theme_minimal() +
  labs(x = "Number of cherries")+ 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,11)

#alpha = 0.7, beta=0.3
p10 <- ggplot(data.frame(int_length = beta2_20_0.7_0.3$int_length),
              aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Internal tree length",y="count")+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_blank())+
  xlim(20,140)



p11 <- ggplot(data.frame(total_length = beta2_20_0.7_0.3$total_length),
              aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Total tree length", title=expression(bold("Beta2 (" ~alpha == 0.7~","~ beta==0.3 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(56,400)



p12 <- ggplot(data.frame(cherry = beta2_20_0.7_0.3$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") + 
  theme_minimal() +
  labs(x = "Number of cherries")+ 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,11)


#alpha = 0.5, beta=0.5
p13 <- ggplot(data.frame(int_length = beta2_20_0.5_0.5$int_length),
              aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Internal tree length",y="count")+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7))+
  xlim(20,140)



p14 <- ggplot(data.frame(total_length = beta2_20_0.5_0.5$total_length),
              aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Total tree length", title=expression(bold("Beta2 (" ~alpha == 0.5~","~ beta==0.5 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size=7),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(56,400)



p15 <- ggplot(data.frame(cherry = beta2_20_0.5_0.5$cherry),
              aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") + 
  theme_minimal() +
  labs(x = "Number of cherries")+ 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.ticks = element_blank() )+
  xlim(0,11)




(p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9)/(p10|p11|p12)/(p13|p14|p15) & theme(plot.margin = margin(2,2,2,2))

#min(beta2_20_1.1_1.5$cherry,beta2_20_1.5_1.1$cherry,beta2_20_0.3_0.7$cherry,beta2_20_0.7_0.3$cherry,beta2_20_0.5_0.5$cherry)
#max(beta2_20_1.1_1.5$cherry,beta2_20_1.5_1.1$cherry,beta2_20_0.3_0.7$cherry,beta2_20_0.7_0.3$cherry,beta2_20_0.5_0.5$cherry)
