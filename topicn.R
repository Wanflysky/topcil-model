library(ggplot2)
library(gridExtra)
library(factoextra)
library(cluster)
library(psych)
library(corrplot)
library(BiocManager)
library(tidyverse)
library(phyloseq)
library(topicmodels)
library(tidytext)
library(ldatuning)
library(cowplot)
library(dplyr)
library(tibble)
library(MicrobiomeStat)
library(readxl)
meta_df <- read_excel("Chinese.xlsx")
mydata <- meta_df %>%
  select(Sample_ID,Sample_part,Age,sex,Porphyrin,Skincare)%>%
  na.omit(.)

otu_skin<-read_excel("skinspecies_c.xlsx")%>%
  tibble::column_to_rownames(var = "sampleId")
otu_skin_n<- sweep(otu_skin, 2, colSums(otu_skin), `/`)
otu_skin_n <- sweep(otu_skin_n,2,1000000,'*')
tax_skin_n<-read_excel("tax_c.xlsx")%>%
  tibble::column_to_rownames(var = "species")

mydata_n = column_to_rownames(mydata,var = "Sample_ID")

ps = phyloseq(sample_data(mydata_n),
                 otu_table(otu_skin_n, taxa_are_rows = TRUE),
                 tax_table(as.matrix(tax_skin_n)))
minTotRelAbun <- 1e-5           
x <- taxa_sums(ps)
keepTaxa <- (x / sum(x)) > minTotRelAbun
(ps <- prune_taxa(keepTaxa, ps))

(ps_g1 <- tax_glom(ps, taxrank = "genus"))
view(tax_table(ps_g1))
a<-as.data.frame(tax_table(ps_g1))
a$genus=ifelse(a$genus=='unknown',rownames(a),a$genus)

taxa_names(ps_g1) <- a[, 6]
count_matrix_g1 <- data.frame(t(data.frame(otu_table(ps_g1))))
rng = rownames(count_matrix_g1)
count_matrix_g1 <- lapply(count_matrix_g1, function(x) round(x,0))
count_matrix_g1<-as.data.frame(count_matrix_g1)
rownames(count_matrix_g1) = rng

result_g1 <- FindTopicsNumber(
  count_matrix_g1,
  topics = seq(from = 2, to = 25, by = 1),
  metrics = c("CaoJuan2009","Arun2010"),
  method = "VEM",
  control = list(seed = 18),
  mc.cores = 4,
  verbose = TRUE
)
(pn = my_plot(result_g1))
lda_k9 <- LDA(count_matrix_g1, k = 9, method = "VEM", control = list(seed = 243))

b_n_df <- data.frame(tidy(lda_k9, matrix = "beta"))

g_n_df <- data.frame(tidy(lda_k9, matrix = "gamma")) %>%
  arrange(document, topic)

lib_n_size_df <- data.frame(sample_sums(ps_g1)) %>%
  dplyr::rename("read_count" = "sample_sums.ps_g1.") %>%
  rownames_to_column(var = "document")

tm_n_df <- left_join(lib_n_size_df, g_n_df) %>%
  mutate(topic_count = read_count * gamma,
         topic_count = round(topic_count, 0)) %>%
  dplyr::select(-read_count, -gamma) %>%
  pivot_wider(names_from = topic, values_from = topic_count) %>%
  dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
  column_to_rownames(var = "document") %>%
  t(.) %>%
  data.frame(.,check.names = FALSE)

###组成图
b_plot1 <- b_n_df %>%
  dplyr::filter(topic == 1) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot1 <- b_plot1[order(b_plot1$beta, decreasing = TRUE), ]
temp1 = b_plot1[,2:3]
temp1$beta = round(temp1$beta,3) 

b_plot2 <- b_n_df %>%
  dplyr::filter(topic == 2) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot2 <- b_plot2[order(b_plot2$beta, decreasing = TRUE), ]
temp2 = b_plot2[,2:3]
temp2$beta = round(temp2$beta,3) 

b_plot3 <- b_n_df %>%
  dplyr::filter(topic == 3) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot3 <- b_plot3[order(b_plot3$beta, decreasing = TRUE), ]
temp3 = b_plot3[,2:3]
temp3$beta = round(temp3$beta,3) 

b_plot4 <- b_n_df %>%
  dplyr::filter(topic == 4) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot4 <- b_plot4[order(b_plot4$beta, decreasing = TRUE), ]
temp4 = b_plot4[,2:3]
temp4$beta = round(temp4$beta,3) 

b_plot5 <- b_n_df %>%
  dplyr::filter(topic == 5) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot5 <- b_plot5[order(b_plot5$beta, decreasing = TRUE), ]
temp5 = b_plot5[,2:3]
temp5$beta = round(temp5$beta,3) 

b_plot6 <- b_n_df %>%
  dplyr::filter(topic == 6) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot6 <- b_plot6[order(b_plot6$beta, decreasing = TRUE), ]
temp6 = b_plot6[,2:3]
temp6$beta = round(temp6$beta,3) 

b_plot7 <- b_n_df %>%
  dplyr::filter(topic == 7) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot7 <- b_plot7[order(b_plot7$beta, decreasing = TRUE), ]
temp7 = b_plot7[,2:3]
temp7$beta = round(temp7$beta,3) 

b_plot8 <- b_n_df %>%
  dplyr::filter(topic == 8) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot8 <- b_plot8[order(b_plot8$beta, decreasing = TRUE), ]
temp8 = b_plot8[,2:3]
temp8$beta = round(temp8$beta,3) 

b_plot9 <- b_n_df %>%
  dplyr::filter(topic == 9) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 5) %>%
  arrange(desc(term)) %>%
  mutate(term = factor(term, levels = term))
b_plot9 <- b_plot9[order(b_plot9$beta, decreasing = TRUE), ]
temp9 = b_plot9[,2:3]
temp9$beta = round(temp9$beta,3) 

bar_plot1 <- ggplot(temp1, aes(x = term, y = beta, fill = term)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = beta), vjust = 1, size = 3, color = "black") +
  labs(x = "\nTopic1", y = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 5,0,face = "bold", hjust = 0.5),legend.position = "none")

pn1 <- ggplot(data = temp1, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic1", y = "") +
  theme(legend.position = "none")
pn2 <- ggplot(data = temp2, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic2", y = "") +
  theme(legend.position = "none")
pn3 <- ggplot(data = temp3, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic3", y = "") +
  theme(legend.position = "none")
pn4 <- ggplot(data = temp4, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic4", y = "") +
  theme(legend.position = "none")
pn5 <- ggplot(data = temp5, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic5", y = "") +
  theme(legend.position = "none")
pn6 <- ggplot(data = temp6, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic6", y = "") +
  theme(legend.position = "none")
pn7 <- ggplot(data = temp7, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic7", y = "") +
  theme(legend.position = "none")
pn8 <- ggplot(data = temp8, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic8", y = "") +
  theme(legend.position = "none")
pn9 <- ggplot(data = temp9, aes(x = beta, y = term, color = term)) +
  geom_point(aes(size = beta)) +
  labs(x = "\nTopic9", y = "") +
  theme(legend.position = "none")


gg <- grid.arrange(bar_plot1,bar_plot2,bar_plot3,bar_plot4,bar_plot5,bar_plot6,bar_plot7,bar_plot8,bar_plot9, ncol=3, nrow=3)
gn = grid.arrange(pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pn9,ncol=3, nrow=3)

print(bar_plot1)

###
g_n_plot <- data.frame(tidy(lda_k9, matrix = "gamma"))

topic_g1<-subset(g_n_plot, topic == 1)%>%
  dplyr::select(document, gamma)%>%
  rename(topic1 = gamma)  
topic_g2<-subset(g_n_plot, topic == 2)%>%
  dplyr::select(document, gamma)%>%
  rename(topic2 = gamma)
topic_g3<-subset(g_n_plot, topic == 3)%>%
  dplyr::select(document, gamma)%>%
  rename(topic3 = gamma)
topic_g4<-subset(g_n_plot, topic == 4)%>%
  dplyr::select(document, gamma)%>%
  rename(topic4 = gamma)
topic_g5<-subset(g_n_plot, topic == 5)%>%
  dplyr::select(document, gamma)%>%
  rename(topic5 = gamma)
topic_g6<-subset(g_n_plot, topic == 6)%>%
  dplyr::select(document, gamma)%>%
  rename(topic6 = gamma)
topic_g7<-subset(g_n_plot, topic == 7)%>%
  dplyr::select(document, gamma)%>%
  rename(topic7 = gamma)
topic_g8<-subset(g_n_plot, topic == 8)%>%
  dplyr::select(document, gamma)%>%
  rename(topic8 = gamma)
topic_g9<-subset(g_n_plot, topic == 9)%>%
  dplyr::select(document, gamma)%>%
  rename(topic9 = gamma)


topic_n <- left_join(topic_g1,topic_g2,by = "document")
topic_n <- left_join(topic_n,topic_g3,by = "document")
topic_n <- left_join(topic_n,topic_g4,by = "document")
topic_n <- left_join(topic_n,topic_g5,by = "document")
topic_n <- left_join(topic_n,topic_g6,by = "document")
topic_n <- left_join(topic_n,topic_g7,by = "document")
topic_n <- left_join(topic_n,topic_g8,by = "document")
topic_n <- left_join(topic_n,topic_g9,by = "document")



topic_n = column_to_rownames(topic_n,var = "document")
topic_n <- as.data.frame(t(topic_n))

result_n <- dist(topic_n, method = "euclidean")
result_nhc <- hclust(d = result_n, method = "ward.D2")
plot(result_nhc)
hcn <- agnes(topic_n, method = "ward")
pltree(hcn, cex = 0.5, hang = -1,main = "Dendrogram of topics")
fviz_dend(result_nhc,k=4,
          cex = 0.5, 
          k_colors = c("#ef1828","#f88421","#ffbc14",
                       "#00bdcd","#006b7b"),
          color_labels_by_k = TRUE, 
          rect = TRUE, 
          main = "                                                                   cluster of  soil topics "
)
bn<-scale(topic_n)
set.seed(123)
fviz_nbclust(bn,kmeans,method="wss",k.max=8) 
res<-kmeans(bn,4)
res1<-cbind(topic_n,res$cluster)
fviz_cluster(res,data=topic_n[,1:ncol(topic_n)-1])
pca_result = prcomp(as.data.frame(t(topic_n)))
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca_result, col.var = "black")
###
Topics_matrixn <- t(topic_n)
t1 = rownames_to_column(as.data.frame(Topics_matrixn),var="Id")

dndata  = column_to_rownames(b,var = "Id")

f_id <- mydata$Sample_ID
t2 = t1%>%filter(Id %in% f_id)
t2 = t2%>%
  rename(Sample_ID = Id)
t3 = left_join(t2,mydata,by ="Sample_ID")
t3 = t3%>% 
  mutate(Sample_Part = case_when(Sample_part == "forehead" ~ 1,
                                 Sample_part == "nose" ~ 2,
                                 Sample_part == "cheek" ~ 3,
  ))
tdata = t3[,-11]
tdata = column_to_rownames(tdata, var = "Sample_ID")

cnskin = corr.test(tdata,method = "pearson",adjust = "none")
cnskin1 = corr.test(tdata,method = "spearman",adjust = "none")

cnskin_c = cnskin$r
cnskin_c1 = cnskin1$r
cnskin_c = cnskin_c[c(10,12:17),1:9]
cnskin_c1 = cnskin_c1[c(11,18),1:9]
cnskin_r = cnskin$p
cnskin_r1 = cnskin1$p
cnskin_r = cnskin_r[c(10,12:17),1:9]
cnskin_r1 = cnskin_r1[c(11,18),1:9]

cnskin_c_all = rbind(cnskin_c,cnskin_c1)
cnskin_r_all = rbind(cnskin_r,cnskin_r1)

col1 <- colorRampPalette(c("blue","white","red"))
corrplot(cnskin_c_all,method ="color",col=col1(100),
         p.mat = cnskin_r_all,  diag = T, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.9, insig = 'label_sig', pch.col = 'black',tl.col="black")
###
