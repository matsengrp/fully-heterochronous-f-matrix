source("beta2.R")


#bernoulli_1000_5.10_1=beta2_sample(sample_size = 1000,tree_size = 5,alpha=10,beta=1)
bernoulli_1000_20.10_1=beta2_sample(sample_size = 1000,tree_size = 20,alpha=10,beta=1)
#bernoulli_1000_50.10_1=beta2_sample(sample_size = 1000,tree_size = 50,alpha=10,beta=1)


#bernoulli_1000_5.10_10=beta2_sample(sample_size = 1000,tree_size = 5,beta=10,alpha=10)
bernoulli_1000_20.10_10=beta2_sample(sample_size = 1000,tree_size = 20,beta=10,alpha=10)
#bernoulli_1000_50.10_10=beta2_sample(sample_size = 1000,tree_size = 50,beta=10,alpha=10)


#bernoulli_1000_5.1_10=beta2_sample(sample_size = 1000,tree_size = 5,alpha =1, beta=10)
bernoulli_1000_20.1_10=beta2_sample(sample_size = 1000,tree_size = 20,alpha=1, beta=10)
#bernoulli_1000_50.1_10=beta2_sample(sample_size = 1000,tree_size = 50,alpha=1,beta=10)












# alpha=10 beta=1
#lightblue in colored version
p1 <- ggplot(data.frame(int_length = bernoulli_1000_20.10_1$int_length),
             aes(x = int_length)) +
  #geom_histogram(bins = 6, color = "gray80", fill = "gray80") +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(y="count")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7),)+
  xlim(0,132)



p2 <- ggplot(data.frame(total_length = bernoulli_1000_20.10_1$total_length),
             aes(x = total_length)) +
  #geom_histogram(bins = 14, color = "gray80", fill = "gray80") +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "total tree length", title=expression(bold("Beta2 (" ~alpha == 10~","~ beta==1 ~ ")")))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.ticks = element_blank() )+
  xlim(56,400)



p3 <- ggplot(data.frame(cherry = bernoulli_1000_20.10_1$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank() )+
  xlim(0,11)

# alpha =10, beta=10 
p4 <- ggplot(data.frame(int_length = bernoulli_1000_20.10_10$int_length),
             aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "internal tree length",y="count")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7), )+
  xlim(0,132)

p5 <- ggplot(data.frame(total_length = bernoulli_1000_20.10_10$total_length),
             aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(x = "total tree length", title=expression(bold("Beta2 (" ~alpha == 10~","~ beta==10 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.ticks = element_blank() )+
  xlim(56,400)

p6 <- ggplot(data.frame(cherry = bernoulli_1000_20.10_10$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "number of cherries",y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(0,11)


# beta =1 alpha =10

p7 <- ggplot(data.frame(int_length = bernoulli_1000_20.1_10$int_length),
             aes(x = int_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Internal tree length",y="count")+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),)+
  xlim(0,132)

p8 <- ggplot(data.frame(total_length = bernoulli_1000_20.1_10$total_length),
             aes(x = total_length)) +
  geom_bar(color = "gray80", fill = "gray80") +
  theme_minimal() +
  labs(x = "Total tree length", title=expression(bold("Beta2 (" ~alpha == 1~","~ beta==10 ~ ")")))+
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size=7),
        axis.title.y = element_blank(),
        axis.ticks = element_blank() )+
  xlim(56,400)



p9 <- ggplot(data.frame(cherry = bernoulli_1000_20.1_10$cherry),
             aes(x = cherry)) +
  geom_bar(color = "gray80", fill = "gray80") + 
  theme_minimal() +
  labs(x = "Number of cherries")+ 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.ticks = element_blank() )+
  xlim(0,11)




(p1|p2|p3)/(p4|p5|p6)/(p7|p8|p9) & theme(plot.margin = margin(2,2,2,2))


