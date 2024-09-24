library(tidyverse)
library(mosaic)
library(ggpubr)
library(ggsci)
library(caret)

#import the data
df<-as_tibble(read.table(file = "pearson.input.csv",header = T,sep = ","))
glimpse(df)
tmp<-as_tibble(read.table("03112024/3_preprocess_zerovar_bagimputed_dataset.tsv",sep = "\t",header = T))
tmp

glimpse(tmp)
# DATA WRANGLING
#remove the last column that doesnt have anything in it
df<-df%>%select(-c(X))
#add an index column
df<-df%>%mutate(pt_id=paste0("pt_",rownames(df)))
df$Birth_weight
glimpse(df)
colnames(df)
#convert the appropriate columns to factor class
#ordered class
df_ordered<-df%>%
  select(pt_id,APGAR_1min,APGAR_5min,SDI_score)%>%
  mutate(APGAR_1min=as.ordered(APGAR_1min),
         APGAR_5min=as.ordered(APGAR_5min),
         SDI_score=as.ordered(SDI_score))
#factor class
df_factor<-df%>%
  select(pt_id,Race,Ethnicity,Offspring_Biological_Sex,Delivery_Route,Superfund_Status)%>%
  mutate(Superfund_Status=gsub(1,TRUE,Superfund_Status),
         Superfund_Status=gsub(2,FALSE,Superfund_Status),
         Race_Black=gsub(4,TRUE,Race),
         Race_Black=gsub(5,FALSE,Race_Black),
         Race_White=gsub(4,TRUE,Race),
         Race_White=gsub(5,FALSE,Race_White),
         Ethnicity_Hispanic=gsub(1,TRUE,Ethnicity),
         Ethnicity_Hispanic=gsub(2,FALSE,Ethnicity_Hispanic),
         Sex_Male=gsub(1,FALSE,Offspring_Biological_Sex),
         Sex_Male=gsub(2,TRUE,Sex_Male),
         Sex_Female=gsub(2,FALSE,Offspring_Biological_Sex),
         Sex_Female=gsub(1,TRUE,Sex_Female),
         Delivery_Spontaneous=gsub(1,TRUE,Delivery_Route),
         Delivery_Spontaneous=gsub(2,FALSE,Delivery_Spontaneous),
         Delivery_Cesarean=gsub(2,TRUE,Delivery_Route),
         Delivery_Cesarean=gsub(1,FALSE,Delivery_Cesarean))%>%
  select(-c(Race,Ethnicity,Offspring_Biological_Sex,Delivery_Route))%>%
  mutate(across(-c(pt_id),as.logical))
#numeric class
glimpse(df)
df_numeric<-df%>%select(pt_id,BaP:Maternal_Age,Length,Birth_weight,GA_w)
df_ordered$Superfund_Status

class(df$Birth_weight)

df_numeric$Birth_weight
#join them together
df2<-full_join(df_numeric,df_factor)
df2<-full_join(df2,df_ordered)
#lets remove them for now
# TODO ask about  SumFTIR and SumPy.GC.MS
#lets remove them for now
df3<-df2%>%select(-c(SumPy.GC.MS,SumPAHs))
#bring back in the FTIR column
df3
###################### Plot the data ###########################
# #make a composite score of the PAH to use as a target
# df3_long<-df3%>%pivot_longer(cols = c(BaP,BbF,DBA),names_to = "PAH",values_to = "PAH_count")
# df3_long
# gghistogram(data = df3_long,x = "PAH_count",color = "PAH",fill = "PAH",bins = 80)#+xscale("log10")
#preprocess the data
#convert it to a data frame with rownames
df3.5<-data.frame(df3%>%select(-pt_id))
df3.5
rownames(df3.5)<-df3$pt_id
df2$Birth_weight
preproc<-preProcess(df3.5,
                    method = c("center","scale","knnImpute", "corr", "zv", "nzv"),
                    thresh = 0.95,
                    cutoff = 0.95,
                    verbose = T)
summary(preproc)
df3.5
df4<-predict(preproc,df3.5)
df4
# vec<-data.frame("a"=colnames(df3.5))
# vec2<-data.frame("a"=colnames(df4))
# vec3<-anti_join(vec,vec2)
# vec3
# we dropped these three variables:
# a
# 1 Styrene_Butadiene
# 2     Polycarbonate
# 3      Polyurethane
# > 

df4.5<-df4%>%
  mutate(PAH_composite=BaP+BbF+DBA,
         microplastics_composite=Polyethylene+
           Polyvinyl_Chloride+
           Nylon+
           Polyethylene_Terephthalate+
           Polypropylene+
           Polystyrene+
           Acrylonitrile_Butadiene_Styrene+
           PMMA+
           N6)
glimpse(df4.5)
colnames(df4.5)
colnames(df5)
df5.1<-df4.5%>%
  dplyr::select(-c(Sex_Male,
                   Delivery_Spontaneous,
                   Delivery_Cesarean,
                   Ethnicity_Hispanic,
                   Race_Black,
                   Length,
                   Maternal_Age,
                   Birth_weight,
                   SDI_score,
                   APGAR_5min,
                   GA_w,
                   Race_White,
                   Sex_Female,
                   SumFTIR,))
colnames(df5.1)
df5$BaP
cor()
d<-cor(df5.1%>%select(where(is.numeric)))
d2<-cor.mtest(df5.1%>%select(where(is.numeric)))

range(d)
corrplot(d,
         method = "square",
         p.mat = d2$p,
         insig = "label_sig",
         tl.col = "black", 
         tl.cex = 1,
         pch.col = 'black',
         is.corr = T,
         order = "AOE",
         type = "lower",
         pch.cex = 1)
         
# write.table(d,"correlation_values.tsv",sep = "\t",row.names = F)
# df5<-df4%>%mutate(across(where(is.ordered),as.numeric))


# write.table(d,"pearson_results.tsv",sep = "\t")

# ######### FIGURE THYME ####
# ggscatter(df4.5,x = "PAH_composite",y = "Birth_weight",color = "Birth_weight",
#           cor.coef = T,
#           cor.coeff.args = list(label.x.npc="center"),
#           conf.int = T,
#           add = "reg.line",
#           # cor.coeff.args = "pearson"
# )+
#    # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# 
# ggsave(filename = "PAH_composite_birthweight_not_fixed.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# ggscatter(df,x = "DBA",y = "Birth_weight",color = "Birth_weight",
#           cor.coef = T,
#           cor.coeff.args = list(label.x.npc="center"),
#           conf.int = T,
#           add = "reg.line",
#           # cor.coeff.args = "pearson"
# )+
#   # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# 
# ggsave(filename = "DBA_birthweight.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# ggscatter(df,x = "BaP",y = "Birth_weight",color = "Birth_weight",
#           cor.coef = T,
#           cor.coeff.args = list(label.x.npc="center"),
#           conf.int = T,
#           add = "reg.line",
#           # cor.coeff.args = "pearson"
# )+
#    # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# 
# ggsave(filename = "BaP_birthweight.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# 
# ggscatter(df,x = "DBA",
#           y = "GA_w",color = "GA_w",
#           cor.coef = T,
#           cor.coeff.args = list(label.x.npc="center"),
#           conf.int = T,
#           add = "reg.line",
#           # cor.coeff.args = "pearson"
# )+
#   # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# 
# 
# 
# ggsave(filename = "DBA_Ga_W.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# 
# ggscatter(df,x = "Nylon",
#           y = "APGAR_1min",
#           color = "APGAR_1min",
#           cor.coef = T,margin.plot = "density",
#           cor.coeff.args = list(label.x.npc="center",method="spearman"),
#           conf.int = T,
#           add = "reg.line",
#           palette = "gsea",
#            # cor.coeff.args = "spearman"
# )+
#   # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# 
# 
# 
# ggsave(filename = "Nylon_APGAR1min.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# 
# ggscatter(df,x = "Polypropylene",
#           y = "APGAR_1min",
#           color = "APGAR_1min",
#           cor.coef = T,margin.plot = "density",
#           cor.coeff.args = list(label.x.npc="center",method="spearman"),
#           conf.int = T,
#           add = "reg.line",
#           palette = "gsea",
#           # cor.coeff.args = "spearman"
# )+
# # coord_fixed()+
# scale_color_viridis_c()+theme_pubr(legend = "right")
# ggsave(filename = "Polypropylene_APGAR1min.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# 
# df$SDI_score
# ggscatter(df,x = "Nylon",
#           y = "SDI_score",
#           color = "SDI_score",
#           cor.coef = T,margin.plot = "density",
#           cor.coeff.args = list(label.x.npc="center",method="spearman"),
#           conf.int = T,
#           add = "reg.line",
#           palette = "gsea",
#           # cor.coeff.args = "spearman"
# )+
#   # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# ggsave(filename = "Nylon_SDI_score.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# 
# ggscatter(df,x = "Polyethylene",
#           y = "SDI_score",
#           color = "SDI_score",
#           cor.coef = T,margin.plot = "density",
#           cor.coeff.args = list(label.x.npc="center",method="spearman"),
#           conf.int = T,
#           add = "reg.line",
#           palette = "gsea",
#           # cor.coeff.args = "spearman"
# )+
#   # coord_fixed()+
#   scale_color_viridis_c()+theme_pubr(legend = "right")
# ggsave(filename = "Polyethylene_SDI_score.png",
#        dpi = 300,
#        device = "png",
#        width = 10,
#        height = 6,units = "in")
# 
# cor.test(df$Polyethylene,df$SDI_score,method="spearman")
################ ITS ON THE WRONG END OF p=0.05 #######################
# Spearman's rank correlation rho
# 
# data:  x and y
# S = 19639, p-value = 0.05018
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2937348 

##################### create the dissimilarity matrix  ##########################

library(cluster)

#' Compute Gower distance
# df4<-df4%>%mutate(Superfund_Status=as.logical(Superfund_Status))
colnames(df5)
df5.1
gower_dist <- daisy(df5.1, metric = "gower")
gower_dist
gower_mat <- as.matrix(gower_dist)
gower_mat2<-as_tibble(gower_mat)%>%mutate(pt=rownames(gower_mat))
corrplot(gower_mat,)
M<-gower_mat
corrplot(M, order = 'hclust',type="full",method="square", hclust.method = 'ward.D2', addrect = 4)

#' Print most similar patients
most_disparate<-df[which(gower_mat2 == min(gower_mat[gower_mat != min(gower_mat)]), arr.ind = TRUE)[1, ], ]


#' Print most dissimilar patients
most_similar<-df[which(gower_mat2 == max(gower_mat[gower_mat != max(gower_mat)]), arr.ind = TRUE)[1, ], ]
glimpse(most_similar)

# lets look for a number of clusters that are both meaningful and easy to remember, 
#i.e. 2 to 8 maximum. The silhouette figure helps us identify the best option(s).
set.seed(314)
#
sil_width <- c(NA)
for(i in 2:8){  
  pam_fit <- pam(gower_dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}

dev.off()
plot(1:8, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:8, sil_width)

k <- 4
pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- df %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
pam_results[2,2]
######### Visualize your clusters
library(Rtsne)
tsne_obj2 <- Rtsne(X = gower_dist, is_distance = TRUE, perplexity = 15,max_iteration=1000)
saveRDS(tsne_obj2,"tsne_obj2.RDS")
tsne_data <- tsne_obj2$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))
# library(corrplot)
# distmat<-as.matrix(gower_dist)
# library(pheatmap)
# distmat

# df5.5<-df5%>%
#   mutate(PAH_composite=BaP+BbF+DBA,
#          microplastics_composite=Polyethylene+
#            Polyvinyl_Chloride+
#            Nylon+
#            Polyethylene_Terephthalate+
#            Polypropylene+
#            Polystyrene+
#            Acrylonitrile_Butadiene_Styrene+
#            PMMA+
#            N6)

# pheatmap(distmat,
#          cutree_cols = 5,
#          cutree_rows = 5,
#          cellwidth = 12,
#          cellheight = 12,
#          fontsize = 14,
#          border_color = NA)#,
#          # legend_breaks = c(10,20,30,40,500),
#          annotation_row = df4.5[,c(1:3,26)],
#          annotation_col = df4.5[,c(4:14,27)])
# dev.off()  
# annot_col<-df5.5%>%
#   mutate(across(where(is.logical),as.character))%>%
#   select(c(where(is.ordered),where(is.character)))
# 
# dim(annot_col)
# dim(annot_col)
# 
# 
# pheatmap(distmat,
#          # cluster_rows = F,
#          # cluster_cols = F,
#           cutree_cols = 5,
#          cutree_rows = 5,
#          cellwidth = 12,
#          cellheight = 12,
#          fontsize = 14,
#          border_color = NA,
#          annotation_row = annot_col[6:11],
#          annotation_col =annot_col[1:5])
dev.off()
df_tmp<-cbind(df5,tsne_data)

df_tmp<-df_tmp%>%arrange(desc(PAH_composite))

# df_tmp<-NULL
df_tmp
# ggballoonplot(data = df_tmp,x = "pt_id",y = "PAH_composite")
ggpar(ggscatter(data = df_tmp,
                x = "X",
                y = "Y",
                shape = "cluster",
                color = "cluster",
                palette = "igv",
                # ggtheme = theme_minimal(),
                star.plot = T,
                # margin.plot = "hist",
                # fill = "cluster"),
                legend = "right"))
colnames(df5.1)

ggsave(filename = "tsne_PAM_cluster12123.png",
       dpi = 300,
       device = "png",
       width = 10,
       height = 6,units = "in")
################################################################################

vec<-colnames(df_tmp)
vec
lapply(vec, function(VAR) {
  ggsave(filename = paste0(VAR,"_PAM_121223.png"),
         plot = ggpar(ggscatter(data = df_tmp,
                  x = "X",
                  y = "Y",
                  shape = "cluster",
                  color = VAR,
                  fill = "cluster"),legend = "right"), dpi = 300,
         device = "png",
         width = 10,
         height = 6,units = "in")
})


#### 03122024.R

setwd("/Users/enricobarrozo/Library/CloudStorage/Box-Box/BarrozoLab/1_Mentoring and Publications/Projects/Microplastics/Microplastics-Pilot-Metadata_Research/2.6.24-NewMetadataSearch/3.11.24-Michael_ML-reanalysis/03112024/")

set.seed(seed=314)


library(tidyverse) # install.packages("tidyverse")
library(mosaic) # install.packages("mosaic")
library(ggpubr) # install.packages("ggpubr")
library(ggsci) # install.packages("ggsci")
library(caret) # install.packages("caret")
library(vegan) # install.packages("vegan")
library(MASS) # install.packages("MASS")
library(cluster) # install.packages("cluster")
library(Rtsne) # install.packages("Rtsne")
library(corrplot) # install.packages("corrplot")
library(pheatmap) # install.packages("pheatmap")

setwd("/Users/enricobarrozo/Library/CloudStorage/Box-Box/BarrozoLab/1_Mentoring and Publications/Projects/Microplastics/Microplastics-Pilot-Metadata_Research/2.6.24-NewMetadataSearch/3.11.24-Michael_ML-reanalysis/03112024")
## load michael's data ## erb 3.12.24
load("03112024.RData")
load("/Users/enricobarrozo/Library/CloudStorage/Box-Box/BarrozoLab/1_Mentoring and Publications/Projects/Microplastics/Microplastics-Pilot-Metadata_Research/2.6.24-NewMetadataSearch/3.11.24-Michael_ML-reanalysis/03112024/03112024.RData")


## Try uploading to the server
setwd("/home/ebarrozo/Microplastics/03112024")
load("/home/ebarrozo/Microplastics/03112024/03112024.RData")

load("tsne_obj2.RDS")

data <- readRDS("tsne_obj2.RDS")
#import the data
df<-as_tibble(read.table(file = "pearson.input.csv",header = T,sep = ","))
glimpse(df)
tmp<-as_tibble(read.table("03112024/3_preprocess_zerovar_bagimputed_dataset.tsv",sep = "\t",header = T))
tmp
tmp2<-tmp%>%dplyr::select(-c(Superfund_Site))
tmp2

tmp2<-tmp2%>%mutate(across(where(is.character),as.factor))%>%
  mutate(ZCTA5_population=factor(ZCTA5_population))%>%
  mutate(ZIP=factor(ZIP))%>%
  mutate(APGAR_at_1_minute=ordered(APGAR_at_1_minute))

tmp3<-tmp2%>%dplyr::select(-c(SubjectID,PregnancyID,ZIP))%>%
  rename(fetal_sex=Infant_gender)%>%mutate(across(where(is.numeric),as.double))
tmp3


# 
tmp_num<-tmp3%>%dplyr::select(where(is.numeric))
tmp_notnum<-tmp3%>%dplyr::select(!where(is.numeric))

rownames(tmp_num)<-tmp3$PregnancyID
# 
# 
# 
# tmp_num<-tmp_num%>%mutate(across(where(is.numeric),as.double))
# tmp_num$Maternal_Age
tmp3


library(caret)
preproc<-preProcess(tmp_num,
                     method = c( "zv", "nzv","bagImpute"),
                     thresh = 0.95,
                     cutoff = 0.95,
                     verbose = T)

preproc2<-preProcess(tmp_num,
                    method = c("center","scale","zv", "nzv","bagImpute"),
                    thresh = 0.95,
                    cutoff = 0.95,
                    verbose = T)


# df3.5
# tmp_num2<-stats::predict(preproc,tmp_num)
tmp_num2<-stats::predict(preproc2,tmp_num)
tmp_num2


tmp4<-cbind(tmp_notnum,tmp_num2)%>%dplyr::select(-c(ZCTA5_population, SDI_score))
as_tibble(tmp4)
glimpse(as_tibble(tmp4))

library(cluster)

contaminants<-tmp4[,2:16]
contaminants

#' Compute Gower distance
# df4<-df4%>%mutate(Superfund_Status=as.logical(Superfund_Status))
gower_dist <- daisy(contaminants, metric = "gower")
euclid_dist <- daisy(contaminants, metric = "euclidean")


gower_mat <- as.matrix(gower_dist)
euclid_mat <- as.matrix(euclid_dist)


gower_mat2<-as_tibble(gower_mat)%>%mutate(pt=rownames(gower_mat))
gower_mat2

euclid_mat2<-as_tibble(euclid_mat)%>%mutate(pt=rownames(euclid_mat))
library(pheatmap)
pheatmap(euclid_mat2%>%
           dplyr::select(where(is.numeric)),
         cutree_rows = 3,cutree_cols = 3)


library(corrplot)
corrplot(gower_mat,)
corrplot(euclid_mat2,)

M<-gower_mat
corrplot(M, order = 'hclust',type="full",method="square", hclust.method = 'ward.D2', addrect = 4)

#' Print most similar patients
most_similar<-tmp4[which(gower_mat2 == min(gower_mat[gower_mat != min(gower_mat)]), arr.ind = TRUE)[1, ], ]
glimpse(most_similar)

#' Print most dissimilar patients
most_disparate<-tmp4[which(gower_mat2 == max(gower_mat[gower_mat != max(gower_mat)]), arr.ind = TRUE)[1, ], ]


# lets look for a number of clusters that are both meaningful and easy to remember, 
#i.e. 2 to 8 maximum. The silhouette figure helps us identify the best option(s).
set.seed(314)
#
sil_width <- c(NA)
for(i in 2:8){  
  pam_fit <- pam(euclid_dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}

dev.off()
plot(1:8, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:8, sil_width)

k <- 4
pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- df %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary
pam_results[2,2]
######### Visualize your clusters
library(Rtsne)
tsne_obj2 <- Rtsne(X = gower_dist, 
                   is_distance = TRUE, 
                   perplexity = 25,
                   max_iteration=1000)
saveRDS(tsne_obj2,"03112024/tsne_obj2.RDS")

tsne_obj2 <- readRDS(file = "tsne_obj2.RDS")

## load michael's data? ## erb 3.12.24
load("human_spatial_data-integrated-annotated_v1.RData")
## redo the tsne plot vis at 600 dpi
## redo the env vec fits with 4 clusters at 600 dpi at the 3 p val filters
## redo the vlnplots with 4 clusters
## make the boxplots with 4 clusters for the non-toxicants



tsne_data <- tsne_obj2$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))

tsne_data
df_tmp<-cbind(tmp4,tsne_data)
df_tmp

# df_tmp<-NULL
df_tmp
# ggballoonplot(data = df_tmp,x = "pt_id",y = "PAH_composite")
a<-ggpar(ggscatter(data = df_tmp,
                x = "X",
                y = "Y",
                shape = "cluster",
                color = "cluster",
                palette = "igv",
                # ggtheme = theme_minimal(),
                star.plot = T,
                
                # ellipse = T,ellipse.alpha = 0,ellipse.type = "euclid",
                # margin.plot = "hist",
                # fill = "cluster"),
                legend = "right")+coord_fixed())

a

ggsave(filename = "03112024/tsne_start.plot.corrd_fixed.png",
       device = "png",
       dpi = 600, 
       height=8,
       width = 8,
       units = "in")
################## ENVFIT TIME ##########

library(vegan)
library(MASS)

fit <- envfit(euclid_dist, 
              tmp4,
              perm = 9999, na.rm=TRUE)

fit

vec<-scores(fit, "vectors")
fac<-scores(fit, "factors")
vec
fac

library(ggsci)
colors<-pal_igv()(4)
map=setNames(pal_igv()(4),c(1,2,3,4))

tsne_data<-tsne_data%>%
  mutate(colors=cluster)
tsne_data$colors=map[tsne_data$cluster]
tsne_data

plot(tsne_data$X,tsne_data$Y,col=tsne_data$colors)
plot(fit, col ="black", p.max=0.05)

tsne_data

fit
dev.off()
a<-plot(tsne_data$X,tsne_data$Y,col=tsne_data$colors)+
  plot(fit, p.max = 0.07, col ="black",cex=0.75)
a

ggsave(plot = a,device = "png",
       filename = "03112024/0001.envfit.png",
       width = 8,
       height = 8,
       dpi = 600,
       units = "in")

# ################################################################################
 # #### plot the features ######


############## CUSTOM PLOTS FOR ENRICO #########################################

tmp5<-tmp%>%mutate(cluster=tsne_data$cluster,
                    X=tsne_data$X,
                    Y=tsne_data$Y)
dim(tmp)

annot_col<-data.frame(tmp)%>%
  dplyr::select(APGAR_at_1_minute, Type_of_labor_Spontaneous_Augmented)#%>%
  mutate(Type_of_labor_Spontaneous_Augmented=as.logical(Type_of_labor_Spontaneous_Augmented))

  
  dim(tmp5)
dim(annot_col)

annot_col<-cbind(annot_col,tmp5%>%dplyr::select(cluster))
annot_col
dim(annot_col)
dim(tmp4)
tmp
glimpse(tmp[,4:18])
mat<-as.matrix( tmp[,4:18])
mat
tmp4


dim(mat)
rownames(annot_col)
rownames(mat)<-tmp$PregnancyID
rownames(annot_col)<-tmp$PregnancyID


rownames()
annot_col<-annot_col%>%mutate()

pheatmap(mat = t(mat),
         cellwidth = 6.5,
         cellheight = 10,
         border_color = "black",
         cluster_rows = F,cluster_cols = F,
         cutree_rows = 2,
         cutree_cols = 4,
         annotation_col =annot_col)
glimpse(tmp5[,4:18])

tmp$Type_of_labor_Spontaneous_Augmented
spot_perc<-data.frame(tally(cluster~Type_of_labor_Spontaneous_Augmented,tmp5,"perc"))
spot_perc

spot<-mosaic::tally(cluster~Type_of_labor_Spontaneous_Augmented,tmp5,"data.frame")
spot
ggbarplot(data = spot_perc,
          # color = "cluster",
          x = "cluster",
          y = "Freq",
          color = "Type_of_labor_Spontaneous_Augmented",
          fill = "Type_of_labor_Spontaneous_Augmented",
          position = position_dodge(), facet.by ="Type_of_labor_Spontaneous_Augmented", 
          palette = "igv")



ggbarplot(data = spot_perc,
          color = "cluster",
          fill = "cluster",
          y = "Freq",
          x = "cluster",
          # fill = "Type_of_labor_Spontaneous_Augmented",
          position = position_dodge(), facet.by ="Type_of_labor_Spontaneous_Augmented", 
          palette = "igv")+stat_compare_means()

mosaicplot(spot)

df_plot2<-tmp4%>%
  dplyr::select(Delivery_route_Cesarean,Delivery_route_Spontaneous_Vaginal,Nylon,DBA,cluster)%>%
  pivot_longer(cols = -cluster)

plot=ggpar(p = ggviolin(data = df_plot2,
                        x = "cluster",
                        y = "value",
                        color = "cluster",
                        add = "jitter",
                        palette = "igv",fill="cluster",
                        trim = F,
                        alpha = 0.09,
                        ylab = "",ggtheme = theme_pubr())+
             facet_wrap(~name,nrow = 1,scales = "free",strip.position = "top")+
             # stat_compare_means(method = "anova",label.x.npc = "center",label.y.npc = "top")+
             # stat_compare_means(label.x.npc = "center",label.y.npc = "middle")+
             stat_compare_means(comparisons = my_comparisons,
                                method = "wilcox",
                                aes(label = paste0("p = ", after_stat(p.format)))),legend = "none")
plot

# 
##############


########################################################################

vec<-colnames(tmp4)

df_plot2<-tmp4%>%dplyr::select(where(is.numeric))
vec<-colnames(df_plot2)

vec[1:15]
my_comparisons<-list(c("1","2"),c("2","3"),c("3","4"),c("1","3"),c("2","4"),c("1","4"))


lapply(vec[1:15], function(VAR) {
  ggsave(filename = paste0("03112024/",VAR,"_BOX_03112024.png"),
     
         plot=ggboxplot(data = tmp5,
                       x = "cluster",
                       y = VAR,
                       color = "cluster",
                       add = "jitter",
                       palette = "igv")+
           stat_compare_means(method = "anova",label.x.npc = "center",label.y.npc = "top")+
           # stat_compare_means(label.x.npc = "center",label.y.npc = "middle")+
           stat_compare_means(comparisons = my_comparisons,
                              method = "wilcox",
                              aes(label = paste0("p = ", after_stat(p.format)))),
         dpi = 300,
         device = "png",
         width = 6,
         height = 11,units = "in")
})




############## CUSTOM PLOTS FOR ENRICO #########################################

## data_overview_visualizations_06132024.R
tmp5$PregnancyID
as_tibble(df_vec)

df_vec$PregnancyID
df_vec<-df_tmp2%>%dplyr::select(Pregnancy_ID,cluster,X,Y)%>%
  mutate(PregnancyID=as.numeric(Pregnancy_ID))%>%dplyr::select(-c(Pregnancy_ID))

p<-as_tibble(full_join(tmp5,df_vec))
p
library(tidyverse)
df3_long<-p%>%pivot_longer(cols = c(BaP,BbF,DBA),names_to = "PAH",values_to = "PAH_count")
df3_long$PregnancyID
df3_long_PAH<-df3_long%>%dplyr::select(PregnancyID,PAH,PAH_count)%>%distinct_all()
df3_long
df3_long_PAH
gghistogram(data = df3_long,x = "PAH_count",color = "PAH",fill = "PAH",bins = 80)#+xscale("log10")
#preprocess the data
library(ggpubr)
library(ggsci)
df3_long$PAH
df3_long_tmp<-df3_long%>%
  mutate(PregnancyID=factor(PregnancyID, levels=unique(PregnancyID)))%>%
  group_by(PregnancyID)%>%
  mutate(tot=sum(PAH_count))%>%
  ungroup()%>%
  filter(tot>0)
df3_long_tmp$PregnancyID
df3_long_tmp<-df3_long_tmp%>%arrange(desc(tot))%>%mutate(PregnancyID=as.numeric(as.character(PregnancyID)))#,levels(unique(PregnancyID))))
df3_long_tmp$PregnancyID
df3_long_tmp<-  df3_long_tmp%>%mutate(PregnancyID=factor(PregnancyID, levels =unique(df3_long_tmp$PregnancyID) ))
df3_long_tmp$PregnancyID
df3_long_tmp
ggbarplot(data = df3_long_tmp,
          x = "PregnancyID",
          y = "PAH_count",
          color = "PAH",
          fill = "PAH",
          xlab = "Pregnancy ID",
          ylab = "PAH levels (ng/mg)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()
ggsave(filename = "diversity_of_PAH_06132024.png",
       dpi = 300,
       device = "png",
       width = 18,
       height = 6,units = "in")
df3_long

df3_long_microplastics<-df3_long%>%
  mutate(PregnancyID=factor(PregnancyID))%>%
  pivot_longer(cols = c(Polyethylene:N6),
               names_to = "Microplastics",values_to = "microplastics_count")%>%
  dplyr::select(PregnancyID,Microplastics,microplastics_count)%>%distinct_all()
length(unique(df3_long_microplastics$PregnancyID))

df3_long_microplastics_wide<-df3_long_microplastics%>%pivot_wider(names_from = Microplastics,values_from = microplastics_count)

length(unique(df3_long_microplastics_wide$PregnancyID))

library(pheatmap)
pheatmap(df3_long_microplastics_wide[,2:13],cluster_cols = F,cluster_rows = F)

df3_long_microplastics$PregnancyID



df3_long_microplastics_tmp<-df3_long_microplastics


df3_long_microplastics_tmp[is.na(df3_long_microplastics_tmp)]<-0
df3_long_microplastics_tmp

df3_long_microplastics_tmp<-df3_long_microplastics_tmp%>%
  group_by(PregnancyID)%>%
  mutate(tot=sum(microplastics_count))%>%
  ungroup()%>%
  filter(tot>0)
df3_long_microplastics_tmp
length(unique(df3_long_microplastics_tmp$PregnancyID))


#descending order
df3_long_microplastics_tmp<-df3_long_microplastics_tmp%>%
  mutate(PregnancyID=as.numeric(as.character(PregnancyID)))%>%
  arrange(desc(tot))
df3_long_microplastics_tmp$PregnancyID
df3_long_microplastics_tmp<-df3_long_microplastics_tmp%>%
  mutate(PregnancyID=factor(PregnancyID,levels=unique(df3_long_microplastics_tmp$PregnancyID)))
df3_long_microplastics_tmp$PregnancyID


ggbarplot(data = df3_long_microplastics_tmp,
          x = "PregnancyID",
          y = "microplastics_count",
          color = "Microplastics",
          fill = "Microplastics",
          xlab = "Pregnancy ID",
          ylab = "Microplastics levels (mg/g)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()
length(unique(df3_long_microplastics_tmp$PregnancyID))
ggsave(filename = "diversity_of_microplastics_06132024.png",
       dpi = 300,
       device = "png",
       width = 12,
       height = 6,units = "in")

unique(df3_long_microplastics$microplastics_count)

df3_long_microplastics_tmp

df3_wide_microplastics<-df3_long_microplastics_tmp%>%
  # mutate(microplastics_count=gsub(NA,0,microplastics_count))%>%
  pivot_wider(names_from = Microplastics,values_from = microplastics_count,values_fill = 0)

df3_wide_microplastics
# # Replace NAs with 0s
# df3_wide_microplastics[is.na(df3_wide_microplastics)] <- 0
# glimpse(df3_wide_microplastics%>%filter(pt_id=="pt_6"))

# df3_wide_microplastics
# df3_wide_microplastics<-df3_wide_microplastics%>%
#   group_by(pt_id)%>%
#   mutate(tot=Polyethylene+
#            Polyvinyl_Chloride+
#            Nylon+
#            Polyethylene_Terephthalate+ 
#            Styrene_Butadiene+
#            Polycarbonate+ 
#            Polypropylene+ 
#            Polystyrene+ 
#            Acrylonitrile_Butadiene_Styrene+ 
#            PMMA +
#            Polyurethane+ N6)
# glimpse(df3_wide_microplastics%>%filter(pt_id=="pt_6"))
df3_wide_microplastics
df3_wide_microplastics2<-df3_wide_microplastics%>%
  mutate(across(-c(tot,PregnancyID),function(X) X/tot))
df3_wide_microplastics2$PregnancyID
#convert NA to zero
df3_long_microplastics_perc<-df3_wide_microplastics2%>%
  dplyr::select(-tot)%>%
  pivot_longer(cols = -c(PregnancyID),names_to = "Microplastics")%>%mutate(value=value*100)

ggbarplot(data = df3_long_microplastics_perc,
          x = "PregnancyID",
          y = "value",
          color = "Microplastics",
          fill = "Microplastics",
          xlab = "Placenta",
          ylab = "Proportional Microplastics levels (%)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()

ggsave(filename = "diversity_of_microplastics_perc_06132024.png",
       dpi = 300,
       device = "png",
       width = 12,
       height = 6,units = "in")

########### PAH PROP ###############
unique(df3_long_tmp$PregnancyID)
unique(df3_long_PAH$PregnancyID)
df3_wide_PAH<-df3_long_tmp%>%dplyr::select(PregnancyID,PAH,PAH_count, tot)%>%distinct_all()
df3_wide_PAH

df3_wide_PAH<-df3_wide_PAH%>%
  # mutate(microplastics_count=gsub(NA,0,microplastics_count))%>%
  pivot_wider(names_from = PAH,values_from = PAH_count,values_fill = 0)

df3_wide_PAH
# df3_wide_PAH<-df3_wide_PAH%>%
#   group_by(pt_id)%>%
#   mutate(tot=BaP+BbF+DBA)
# 
# glimpse(df3_wide_PAH%>%filter(pt_id=="pt_6"))
df3_wide_PAH
df3_wide_PAH<-df3_wide_PAH%>%
  mutate(across(-c(PregnancyID,tot),function(X) X/tot))
# glimpse(df3_wide_PAH%>%filter(pt_id=="pt_6"))
#convert NA to zero
df3_wide_PAH
df3_long_PAH_perc<-df3_wide_PAH%>%
  dplyr::select(-tot)%>%
  pivot_longer(cols = -c(PregnancyID),names_to = "PAH")%>%mutate(value=value*100)

df3_long_PAH_perc

ggbarplot(data = df3_long_PAH_perc,
          x = "PregnancyID",
          y = "value",
          color = "PAH",
          fill = "PAH",
          xlab = "Placenta",
          ylab = "Proportional PAH levels (%)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()
ggsave(filename = "diversity_of_PAH_perc06132024.png",
       dpi = 300,
       device = "png",
       width = 20,
       height = 6,units = "in")


###################### BY CLUSTER PAH ##############################################

length(unique(df_vec%>%filter(PregnancyID%in%df3_long_PAH_perc$PregnancyID)))
df3_long_PAH_perc
vec2<-df3_long_PAH_perc%>%dplyr::select(PregnancyID)%>%distinct_all()%>%arrange(PregnancyID)
data.frame(vec2$PregnancyID)
df_vec%>%filter(PregnancyID%in%vec2$PregnancyID)
unique(df_vec$PregnancyID)
df3_long_PAH_perc<-df3_long_PAH_perc%>%mutate(PregnancyID=as.character(PregnancyID))
df3_long_PAH_perc
df3_long_PAH_perc


df_vec<-df_vec%>%mutate(PregnancyID=as.character(PregnancyID))
df3_long_PAH_perc_cluster<-right_join(df_vec,df3_long_PAH_perc%>%mutate(PregnancyID=as.character(PregnancyID)))
length(unique(df3_long_PAH_perc_cluster$PregnancyID))
df3_long_PAH_perc
data.frame(df_vec)

df3_long_PAH_perc_cluster<-df3_long_PAH_perc_cluster%>%group_by(PregnancyID)%>%mutate(tot=sum(value),
                                                                                      rel=(value/tot)*100,
                                                                                      rel_tot=sum(rel))%>%arrange(cluster,PregnancyID,value)

data.frame(df3_long_PAH_perc_cluster)

bar_vec<-as.vector(unique(df3_long_PAH_perc_cluster$PregnancyID))
ggbarplot(data = df3_long_PAH_perc_cluster,
          x = "PregnancyID",
          y = "rel",
          color = "PAH",
          fill = "PAH",
          xlab = "Pregnancy ID",sort.by.groups = T,order =bar_vec,
          ylab = "Proportional PAH levels (%)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()

ggsave(filename = "diversity_of_PAH_perc_ordered_by_cluster_06132024.png",
       dpi = 300,
       device = "png",
       width = 12,
       height = 6,units = "in")

###################### BY CLUSTER MICROPLASTICS##############################################

length(unique(df_vec%>%filter(PregnancyID%in%df3_long_microplastics_perc$PregnancyID)))
df3_long_microplastics_perc
vec2<-df3_long_microplastics_perc%>%dplyr::select(PregnancyID)%>%distinct_all()%>%arrange(PregnancyID)
data.frame(vec2$PregnancyID)
df_vec%>%filter(PregnancyID%in%vec2$PregnancyID)
unique(df_vec$PregnancyID)
df3_long_microplastics_perc<-df3_long_microplastics_perc%>%mutate(PregnancyID=as.character(PregnancyID))
df3_long_microplastics_perc
df3_long_microplastics_perc


df_vec<-df_vec%>%mutate(PregnancyID=as.character(PregnancyID))
df3_long_microplastics_perc_cluster<-right_join(df_vec,df3_long_microplastics_perc%>%mutate(PregnancyID=as.character(PregnancyID)))
length(unique(df3_long_microplastics_perc_cluster$PregnancyID))
df3_long_microplastics_perc
data.frame(df_vec)

df3_long_microplastics_perc_cluster<-df3_long_microplastics_perc_cluster%>%group_by(PregnancyID)%>%mutate(tot=sum(value),
                                                                                      rel=(value/tot)*100,
                                                                                      rel_tot=sum(rel))%>%arrange(cluster,PregnancyID,value)

data.frame(df3_long_microplastics_perc_cluster)

bar_vec<-as.vector(unique(df3_long_microplastics_perc_cluster$PregnancyID))
ggbarplot(data = df3_long_microplastics_perc_cluster,
          x = "PregnancyID",
          y = "rel",
          color = "microplastics",
          fill = "microplastics",
          xlab = "Pregnancy ID",sort.by.groups = T,order =bar_vec,
          ylab = "Proportional microplastics levels (%)",
          position = position_stack(),palette = "igv",
          ggtheme = theme_pubr(legend = "right"))+
  rotate_x_text()

ggsave(filename = "diversity_of_microplastics_perc_ordered_by_cluster_06132024.png",
       dpi = 300,
       device = "png",
       width = 12,
       height = 6,units = "in")
