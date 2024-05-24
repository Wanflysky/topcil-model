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
library(RColorBrewer)
 metadata= as.data.frame(all_genus1)
otu_n<- sweep(otu_n, 2, colSums(otu_n), `/`)
otu_n <- sweep(otu_n,2,1000000,'*')
mydata_n<-otu_n[-c(10648:10648),]%>%na.omit(.)
otu_tab_n <- mydata_n %>%
  rownames_to_column()%>%
  tidyr::separate(rowname, into = c("domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "\\|") %>%
  dplyr::select(-domain,-Kingdom, -Phylum, -Class, -Order, -Family) %>%
  dplyr::mutate(Genus = gsub("g__", "", Genus)) %>%
  distinct(Genus,.keep_all = T) %>%
  tibble::column_to_rownames(var = "Genus")

tax_tab_n <- mydata_n %>%
  rownames_to_column()%>%
  dplyr::select(rowname) %>%
  tidyr::separate(rowname, into = c("domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "\\|") %>%
  dplyr::mutate(domain = gsub("d__", "", domain),
                Kingdom = gsub("k__", "", Kingdom),
                Phylum = gsub("p__", "", Phylum),
                Class = gsub("c__", "", Class),
                Order = gsub("o__", "", Order),
                Family = gsub("f__", "", Family),
                Genus = gsub("g__", "", Genus)) %>%
  dplyr::mutate(spec_row = Genus) %>%
  distinct(Genus,.keep_all = T) %>%
  tibble::column_to_rownames(var = "spec_row")


(ps_n <- phyloseq(sample_data(metadata),
                otu_table(otu_tab_n, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax_tab_n ))))
minTotRelAbun <- 1e-5           
x <- taxa_sums(ps_n)
keepTaxa <- (x / sum(x)) > minTotRelAbun
(ps_n <- prune_taxa(keepTaxa, ps_n))
(ps_n)
count_matrix_n <- data.frame(t(data.frame(otu_table(ps_n))))
rng = rownames(count_matrix_n)
count_matrix_n <- lapply(count_matrix_n, function(x) round(x,0))
count_matrix_n<-as.data.frame(count_matrix_n)
rownames(count_matrix_n) = rng

result <- FindTopicsNumber(
  count_matrix_n,
  topics = seq(from = 2, to = 25, by = 1),
  metrics = c("CaoJuan2009", "Arun2010"),
  method = "VEM",
  control = list(seed = 243),
  mc.cores = 4,
  verbose = TRUE
)


lda_k9_n <- LDA(count_matrix_n, k = 9, method = "VEM", control = list(seed = 243))
b_n_df <- data.frame(tidy(lda_k9, matrix = "beta"))

g_n_df <- data.frame(tidy(lda_k9, matrix = "gamma")) %>%
  arrange(document, topic)

lib_n_size_df <- data.frame(sample_sums(ps_n)) %>%
  dplyr::rename("read_count" = "sample_sums.ps_n.") %>%
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

data_list = list(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9)
par(mfrow=c(3, 3))


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

gn = grid.arrange(pn1,pn2,pn3,pn4,pn5,pn6,pn7,pn8,pn9,ncol=3, nrow=3)

print(bar_plot1)
###

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

rainbow_colors <- rainbow(9)

fviz_dend(result_nhc,k=9,
          cex = 0.5, 
          k_colors = rainbow_colors,
          color_labels_by_k = TRUE, 
          rect = TRUE, 
          main = "                                                                   cluster of  skin topics "
)

bn<-scale(topic_n)
set.seed(123)
fviz_nbclust(bn,kmeans,method="wss",k.max=8) 
res<-kmeans(bn,9)
res1<-cbind(topic_n,res$cluster)
fviz_cluster(res,data=topic_n[,1:ncol(topic_n)-1])
###
pca_result = prcomp(as.data.frame(t(topic_n)))
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca_result, col.var = "black")
###
violin_plot = g_n_plot
violin_plot$topic = factor(violin_plot$topic)
violin_plot$gamma = violin_plot$gamma*100
ggplot(violin_plot, aes(x = topic , y = gamma, fill = topic)) +
  geom_violin(alpha = 1,trim = T,scale = "width")+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
  theme_bw()
###


Topics_matrixn <- t(topic_n)
t1 =   rownames_to_column(as.data.frame(Topics_matrixn),var="Id")%>%
       dplyr::mutate(Id = gsub("X", "", Id))
metadata = rownames_to_column(metadata,var = "Id")
f_id <- metadata$Id
t2 = t1%>%filter(Id %in% f_id)

t3 = left_join(t2,metadatan,by ="Id")

sdata = column_to_rownames(t3,var = "Id") 

sdata_blank = sdata %>%
            dplyr::filter(Material == "Blank")
sdata_blank_ng = sdata_blank%>%
                dplyr::filter(Compost == "NG")
sdata_blank_gl = sdata_blank%>%
                dplyr::filter(Compost == "GL")

sdata_blank_ng_1 = sdata_blank_ng%>%
  dplyr::filter(Replication == 1)
sdata_blank_ng_2 = sdata_blank_ng%>%
  dplyr::filter(Replication == 2)
sdata_blank_ng_3 = sdata_blank_ng%>%
  dplyr::filter(Replication == 3)
blank_ng_1 = sdata_blank_ng_1[,c("Day","Degradation")]
blank_ng_2 = sdata_blank_ng_2[,c("Day","Degradation")]
blank_ng_3 = sdata_blank_ng_3[,c("Day","Degradation")]
blank_ng <- left_join(blank_ng_1,blank_ng_2,by = "Day")
blank_ng <- left_join(blank_ng,blank_ng_3,by = "Day")
blank_ng$mean =  (blank_ng$Degradation+blank_ng$Degradation.x+blank_ng$Degradation.y)/3
blank_ngf = blank_ng[,c("Day","mean")]
blank_ngf = blank_ngf[order(blank_ngf$Day),]


sdata_PET = sdata %>%
  dplyr::filter(Material == "PET")
sdata_PET_ng = sdata_PET%>%
  dplyr::filter(Compost == "NG")
sdata_PET_gl = sdata_PET%>%
  dplyr::filter(Compost == "GL")

sdata_PET_ng_1 = sdata_PET_ng%>%
  dplyr::filter(Replication == 1)
sdata_PET_ng_2 = sdata_PET_ng%>%
  dplyr::filter(Replication == 2)
sdata_PET_ng_3 = sdata_PET_ng%>%
  dplyr::filter(Replication == 3)
PET_ng_1 = sdata_PET_ng_1[,c("Day","Degradation")]
PET_ng_2 = sdata_PET_ng_2[,c("Day","Degradation")]
PET_ng_3 = sdata_PET_ng_3[,c("Day","Degradation")]
PET_ng <- left_join(PET_ng_1,PET_ng_2,by = "Day")
PET_ng <- left_join(PET_ng,PET_ng_3,by = "Day")
PET_ng$mean =  (PET_ng$Degradation+PET_ng$Degradation.x+PET_ng$Degradation.y)/3
PET_ng = PET_ng[order(PET_ng$Day),]
PET_ngf = PET_ng[,c("Day","mean")]
PET_ngf = PET_ngf[order(PET_ngf$Day),]



sdata_PBAT = sdata %>%
  dplyr::filter(Material == "PBAT")
sdata_PBAT_ng = sdata_PBAT%>%
  dplyr::filter(Compost == "NG")
sdata_PBAT_gl = sdata_PBAT%>%
  dplyr::filter(Compost == "GL")

sdata_PBAT_ng_1 = sdata_PBAT_ng%>%
  dplyr::filter(Replication == 1)
sdata_PBAT_ng_2 = sdata_PBAT_ng%>%
  dplyr::filter(Replication == 2)
sdata_PBAT_ng_3 = sdata_PBAT_ng%>%
  dplyr::filter(Replication == 3)
PBAT_ng_1 = sdata_PBAT_ng_1[,c("Day","Degradation")]
PBAT_ng_2 = sdata_PBAT_ng_2[,c("Day","Degradation")]
PBAT_ng_3 = sdata_PBAT_ng_3[,c("Day","Degradation")]
PBAT_ng <- left_join(PBAT_ng_1,PBAT_ng_2,by = "Day")
PBAT_ng <- left_join(PBAT_ng,PBAT_ng_3,by = "Day")
PBAT_ng$mean =  (PBAT_ng$Degradation+PBAT_ng$Degradation.x+PBAT_ng$Degradation.y)/3
PBAT_ng = PBAT_ng[order(PBAT_ng$Day),]
PBAT_ngf = PBAT_ng[,c("Day","mean")]
PBAT_ngf = PBAT_ngf[order(PBAT_ngf$Day),]
x3 = PBAT_ng$Day
y3 = PBAT_ng$mean
plot(x,y, xlab ="Day",ylab = "Degradation")
fit <- lm(y ~ x)
abline(fit, col = "red")



sdata_cel = sdata %>%
  dplyr::filter(Material == "Cellulose")
sdata_cel_ng = sdata_cel%>%
  dplyr::filter(Compost == "NG")
sdata_cel_gl = sdata_cel%>%
  dplyr::filter(Compost == "GL")

sdata_cel_ng_1 = sdata_cel_ng%>%
  dplyr::filter(Replication == 1)
sdata_cel_ng_2 = sdata_cel_ng%>%
  dplyr::filter(Replication == 2)
sdata_cel_ng_3 = sdata_cel_ng%>%
  dplyr::filter(Replication == 3)
cel_ng_1 = sdata_cel_ng_1[,c("Day","Degradation")]
cel_ng_2 = sdata_cel_ng_2[,c("Day","Degradation")]
cel_ng_3 = sdata_cel_ng_3[,c("Day","Degradation")]
cel_ng <- left_join(cel_ng_1,cel_ng_2,by = "Day")
cel_ng <- left_join(cel_ng,cel_ng_3,by = "Day")
cel_ng$mean =  (cel_ng$Degradation+cel_ng$Degradation.x+cel_ng$Degradation.y)/3
cel_ng = cel_ng[order(cel_ng$Day),]
cel_ngf = cel_ng[,c("Day","mean")]
cel_ngf = cel_ngf[order(cel_ngf$Day),]
x4 = cel_ng$Day
y4 = cel_ng$mean
plot(x,y, xlab ="Day",ylab = "Degradation")
fit <- lm(y ~ x)
abline(fit, col = "red")

points(x2, y2, col = "red", pch = 16)
points(x3, y3, col = "green", pch = 16)
points(x4, y4, col = "blue", pch = 16)
ngf = merge(blank_ngf,PET_ngf,all = TRUE,by = "Day")
ngf = merge(ngf,PBAT_ngf,all = TRUE,by = "Day")
ngf = merge(ngf,cel_ngf,all = TRUE,by = "Day")
colnames(ngf) = c("Day","Blank","PET","PBAT","Cellulose")
data1 = ngf[,c("Day","Blank")]
colnames(data1) = c("x","y")
data2 = ngf[,c("Day","PET")]
colnames(data2) = c("x","y")
data3 = ngf[,c("Day","PBAT")]
colnames(data3) = c("x","y")
data4 = ngf[,c("Day","Cellulose")]
colnames(data4) = c("x","y")
data2$Group <- "PET"
data3$Group <- "PBAT"
data4$Group <- "Cellulose"
pd1= ggplot(data1,aes(x = x,y = y))+
  geom_point()+
  labs(x = "Day", y = "CO2排放/NG")
all_data <- rbind(data2, data3, data4)
all_data$Group = factor(all_data$Group)
pd2 = ggplot(all_data, aes(x=x, y=y,color =Group)) +
    geom_point()+
    labs(x = "Day", y = "Degradation/NG")


###
sdata_blank_gl_1 = sdata_blank_gl%>%
  dplyr::filter(Replication == 1)
sdata_blank_gl_2 = sdata_blank_gl%>%
  dplyr::filter(Replication == 2)
sdata_blank_gl_3 = sdata_blank_gl%>%
  dplyr::filter(Replication == 3)
blank_gl_1 = sdata_blank_gl_1[,c("Day","Degradation")]
blank_gl_2 = sdata_blank_gl_2[,c("Day","Degradation")]
blank_gl_3 = sdata_blank_gl_3[,c("Day","Degradation")]
blank_gl <- left_join(blank_gl_1, blank_gl_2, by = "Day")
blank_gl <- left_join(blank_gl, blank_gl_3, by = "Day")
blank_gl$mean =  (blank_gl$Degradation + blank_gl$Degradation.x + blank_gl$Degradation.y) / 3
blank_gl = blank_gl[order(blank_gl$Day),]
blank_glf = blank_gl[,c("Day","mean")]
blank_glf = blank_glf[order(blank_glf$Day),]

sdata_PBAT_gl_1 = sdata_PBAT_gl%>%
  dplyr::filter(Replication == 1)
sdata_PBAT_gl_2 = sdata_PBAT_gl%>%
  dplyr::filter(Replication == 2)
sdata_PBAT_gl_3 = sdata_PBAT_gl%>%
  dplyr::filter(Replication == 3)
PBAT_gl_1 = sdata_PBAT_gl_1[,c("Day","Degradation")]
PBAT_gl_2 = sdata_PBAT_gl_2[,c("Day","Degradation")]
PBAT_gl_3 = sdata_PBAT_gl_3[,c("Day","Degradation")]
PBAT_gl <- left_join(PBAT_gl_1,PBAT_gl_2,by = "Day")
PBAT_gl <- left_join(PBAT_gl,PBAT_gl_3,by = "Day")
PBAT_gl$mean =  (PBAT_gl$Degradation+PBAT_gl$Degradation.x+PBAT_gl$Degradation.y)/3
PBAT_gl = PBAT_gl[order(PBAT_gl$Day),]
PBAT_glf = PBAT_gl[,c("Day","mean")]
PBAT_glf = PBAT_glf[order(PBAT_glf$Day),]

sdata_PET_gl_1 = sdata_PET_gl%>%
  dplyr::filter(Replication == 1)
sdata_PET_gl_2 = sdata_PET_gl%>%
  dplyr::filter(Replication == 2)
sdata_PET_gl_3 = sdata_PET_gl%>%
  dplyr::filter(Replication == 3)
PET_gl_1 = sdata_PET_gl_1[,c("Day","Degradation")]
PET_gl_2 = sdata_PET_gl_2[,c("Day","Degradation")]
PET_gl_3 = sdata_PET_gl_3[,c("Day","Degradation")]
PET_gl <- left_join(PET_gl_1,PET_gl_2,by = "Day")
PET_gl <- left_join(PET_gl,PET_gl_3,by = "Day")
PET_gl$mean =  (PET_gl$Degradation+PET_gl$Degradation.x+PET_gl$Degradation.y)/3
PET_gl = PET_gl[order(PET_gl$Day),]
PET_glf = PET_gl[,c("Day","mean")]
PET_glf = PET_glf[order(PET_glf$Day),]

sdata_cel_gl_1 = sdata_cel_gl%>%
  dplyr::filter(Replication == 1)
sdata_cel_gl_2 = sdata_cel_gl%>%
  dplyr::filter(Replication == 2)
sdata_cel_gl_3 = sdata_cel_gl%>%
  dplyr::filter(Replication == 3)
cel_gl_1 = sdata_cel_gl_1[,c("Day","Degradation")]
cel_gl_2 = sdata_cel_gl_2[,c("Day","Degradation")]
cel_gl_3 = sdata_cel_gl_3[,c("Day","Degradation")]
cel_gl <- left_join(cel_gl_1,cel_gl_2,by = "Day")
cel_gl <- left_join(cel_gl,cel_gl_3,by = "Day")
cel_gl$mean =  (cel_gl$Degradation+cel_gl$Degradation.x+cel_gl$Degradation.y)/3
cel_gl = cel_gl[order(cel_gl$Day),]
cel_glf = cel_gl[,c("Day","mean")]
cel_glf = cel_glf[order(cel_glf$Day),]

glf = merge(blank_glf,PET_glf,all = TRUE,by = "Day")
glf = merge(glf,PBAT_glf,all = TRUE,by = "Day")
glf = merge(glf,cel_glf,all = TRUE,by = "Day")
colnames(glf) = c("Day","Blank","PET","PBAT","Cellulose")
data1 = glf[,c("Day","Blank")]
colnames(data1) = c("x","y")
data2 = glf[,c("Day","PET")]
colnames(data2) = c("x","y")
data3 = glf[,c("Day","PBAT")]
colnames(data3) = c("x","y")
data4 = glf[,c("Day","Cellulose")]
colnames(data4) = c("x","y")
data2$Group <- "PET"
data3$Group <- "PBAT"
data4$Group <- "Cellulose"
pd3= ggplot(data1,aes(x = x,y = y))+
  geom_point()+
  labs(x = "Day", y = "CO2排放/GL")
all_data <- rbind(data2, data3, data4)
all_data$Group = factor(all_data$Group)
pd4 = ggplot(all_data, aes(x=x, y=y,color =Group)) +
  geom_point()+
  labs(x = "Day", y = "Degradation/GL")

pd = grid.arrange(pd1,pd2,pd3,pd4,ncol=2, nrow=2)

###



sdata_cel_ng_1 = sdata_cel_ng%>%
  dplyr::filter(Replication == 1)


sdata_a  = sdata_blank_ng[c(1:9,12,14)]

sdataf = sdata%>% 
  mutate(Compost = case_when(Compost == "NG" ~ 2,
                             Compost == "GL" ~ 1,
  ))%>% 
  mutate(Material = case_when(Material == "Blank" ~ 0,
                              Material == "PET" ~ 1,
                              Material == "Cellulose" ~ 2,
                              Material == "PBAT" ~ 3,
  ))


cor_post = corr.test(sdataf,method = "pearson",adjust = "none")
com_c = cor_post$r
com_c = com_c[c(10),1:9]
com_p = cor_post$p
com_p = com_p[c(10),1:9]

b_g = sdataf%>%
  dplyr::filter(Material == "0")%>%
  dplyr::filter(Compost == "1")
cor_bg = corr.test(b_g,method = "pearson",adjust = "none")
bg_c = cor_bg$r
bg_c = bg_c[c(12,14),1:9]
rownames(bg_c) = c("Day(Blank_GL)","Degradation(Blank_GL)")
bg_p = cor_bg$p
bg_p = bg_p[c(12,14),1:9]
rownames(bg_p) = c("Day(Blank_GL)","Degradation(Blank_GL)")

b_n = sdataf%>%
  dplyr::filter(Material == "0")%>%
  dplyr::filter(Compost == "2")
cor_bn = corr.test(b_n,method = "pearson",adjust = "none")
bn_c = cor_bn$r
bn_c = bn_c[c(12,14),1:9]
rownames(bn_c) = c("Day(Blank_NG)","Degradation(Blank_NG)")
bn_p = cor_bn$p
bn_p = bn_p[c(12,14),1:9]
rownames(bn_p) = c("Day(Blank_NG)","Degradation(Blank_NG)")

pe_g = sdataf%>%
  dplyr::filter(Material == "1")%>%
  dplyr::filter(Compost == "1")
cor_peg = corr.test(pe_g,method = "pearson",adjust = "none")
peg_c = cor_peg$r
peg_c = peg_c[c(12,14),1:9]
rownames(peg_c) = c("Day(PET_GL)","Degradation(PET_GL)")
peg_p = cor_peg$p
peg_p = peg_p[c(12,14),1:9]
rownames(peg_p) = c("Day(PET_GL)","Degradation(PET_GL)")

pe_n = sdataf%>%
  dplyr::filter(Material == "1")%>%
  dplyr::filter(Compost == "2")
cor_pen = corr.test(pe_n,method = "pearson",adjust = "none")
pen_c = cor_pen$r
pen_c = pen_c[c(12,14),1:9]
rownames(pen_c) = c("Day(PET_NG)","Degradation(PET_NG)")
pen_p = cor_pen$p
pen_p = pen_p[c(12,14),1:9]
rownames(pen_p) = c("Day(PET_NG)","Degradation(PET_NG)")

c_g = sdataf%>%
  dplyr::filter(Material == "2")%>%
  dplyr::filter(Compost == "1")
cor_cg = corr.test(c_g,method = "pearson",adjust = "none")
cg_c = cor_cg$r
cg_c = cg_c[c(12,14),1:9]
rownames(cg_c) = c("Day(Cellulose_GL)","Degradation(Cellulose_GL)")
cg_p = cor_cg$p
cg_p = cg_p[c(12,14),1:9]
rownames(cg_p) = c("Day(Cellulose_GL)","Degradation(Cellulose_GL)")

c_n = sdataf%>%
  dplyr::filter(Material == "2")%>%
  dplyr::filter(Compost == "2")
cor_cn = corr.test(c_n,method = "pearson",adjust = "none")
cn_c = cor_cn$r
cn_c = cn_c[c(12,14),1:9]
rownames(cn_c) = c("Day(Cellulose_NG)","Degradation(Cellulose_NG)")
cn_p = cor_cn$p
cn_p = cn_p[c(12,14),1:9]
rownames(cn_p) = c("Day(Cellulose_NG)","Degradation(Cellulose_NG)")

pb_g = sdataf%>%
  dplyr::filter(Material == "3")%>%
  dplyr::filter(Compost == "1")
cor_pbg = corr.test(pb_g,method = "pearson",adjust = "none")
pbg_c = cor_pbg$r
pbg_c = pbg_c[c(12,14),1:9]
rownames(pbg_c) = c("Day(PBAT_GL)","Degradation(PBAT_GL)")
pbg_p = cor_pbg$p
pbg_p = pbg_p[c(12,14),1:9]
rownames(pbg_p) = c("Day(PBAT_GL)","Degradation(PBAT_GL)")

pb_n = sdataf%>%
  dplyr::filter(Material == "3")%>%
  dplyr::filter(Compost == "2")
cor_pbn = corr.test(pb_n,method = "pearson",adjust = "none")
pbn_c = cor_pbn$r
pbn_c = pbn_c[c(12,14),1:9]
rownames(pbn_c) = c("Day(PBAT_NG)","Degradation(PBAT_NG)")
pbn_p = cor_pbn$p
pbn_p = pbn_p[c(12,14),1:9]
rownames(pbn_p) = c("Day(PBAT_NG)","Degradation(PBAT_NG)")

c_all = rbind(cncom_c,bg_c,bn_c,peg_c,pen_c,cg_c,cn_c,pbg_c,pbn_c)
rownames(c_all)[1] = "Compost"
p_all = rbind(cncom_p,bg_p,bn_p,peg_p,pen_p,cg_p,cn_p,pbg_p,pbn_p)
rownames(p_all)[1] = "Compost"

cor_post0 = corr.test(sdataf,method = "pearson",adjust = "none")

c_alld = c_all[c(1,2,4,6,8,10,12,14,16),]
p_alld = p_all[c(1,2,4,6,8,10,12,14,16),]
col1 <- colorRampPalette(c("blue","white","red"))
corrplot(c_alld,method ="color",col=col1(100),
         p.mat = p_alld,  diag = T, sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.9, insig = 'label_sig', pch.col = 'black',tl.col="black")

data1 = sdata_blank_gl[,c("topic2","Day")] 
data2 = sdata_blank_ng[,c("topic2","Day")]
data3 = sdata_PET_gl[,c("topic2","Day")]
data4 = sdata_PET_ng[,c("topic2","Day")]
data5 = sdata_cel_gl[,c("topic2","Day")]
data6 = sdata_cel_ng[,c("topic2","Day")]
data7 = sdata_PBAT_gl[,c("topic2","Day")]
data8 = sdata_PBAT_ng[,c("topic2","Day")]
data1$group = "blank_gl"
data2$group = "blank_ng"
data3$group = "PET_gl"
data4$group = "PET_ng"
data5$group = "Cellulose_gl"
data6$group = "Cellulose_ng"
data7$group = "PBAT_gl"
data8$group = "PBAT_ng"

all_data <- rbind(data1, data2 ,data3,data4,data5,data6,data7,data8)
all_data$group = factor(all_data$group)
pdt2 = ggplot(all_data, aes(x=Day, y=topic2,color =group)) +
  geom_point()+
  labs(x = "Day", y = "topic2")

data1 = sdata_blank_gl[,c("topic8","Day")] 
data2 = sdata_blank_ng[,c("topic8","Day")]
data3 = sdata_PET_gl[,c("topic8","Day")]
data4 = sdata_PET_ng[,c("topic8","Day")]
data5 = sdata_cel_gl[,c("topic8","Day")]
data6 = sdata_cel_ng[,c("topic8","Day")]
data7 = sdata_PBAT_gl[,c("topic8","Day")]
data8 = sdata_PBAT_ng[,c("topic8","Day")]
data1$group = "blank_gl"
data2$group = "blank_ng"
data3$group = "PET_gl"
data4$group = "PET_ng"
data5$group = "Cellulose_gl"
data6$group = "Cellulose_ng"
data7$group = "PBAT_gl"
data8$group = "PBAT_ng"

all_data <- rbind(data1, data2 ,data3,data4,data5,data6,data7,data8)
all_data$group = factor(all_data$group)
pdt9 = ggplot(all_data, aes(x=Day, y=topic8,color =group)) +
  geom_point()+
  scale_color_manual(values=brewer.pal(8,"Dark2"))
  labs(x = "Day", y = "topic8")