#load data-----
options(scipen = 999)
rm(list=ls())
graphics.off()
setwd("c:\\users\\de\\dropbox\\lehre\\nbt_manuscript\\nbt_dataanalysis")
data <- read.csv("nbt_rawdata.csv", header=T, sep=";")[2:62] #load data
questions <- read.csv("nbt_questions.csv", header = T, sep = ";") #load questions
library(ggplot2)
library(reshape)
library(svglite)
library(multcomp)
library(e1071)
library(PMCMR)
library(effsize)
library(ppcor)
library(MASS)
library(gplots)
library(multcomp)
library(coin)
library(factoextra)

colnames(data)[7:28] <- sapply(13:34, function(i){paste("active_pt_sem", i, sep="_" )}) #rename columns, number of semester in correspondance to begin of studies
colnames(data)[29:61] <- questions$sorting #rename columns, number questions
colnames(data)[1:6] <- c("sex", "age", "begin_stud", "ana", "ptp", "nep") 
data <- data[ ,c(names(data)[1:28],sort(as.numeric(names(data)[29:61])))] #sort questions

data[7:61] <- as.data.frame(sapply(data[7:61], as.numeric)) #change to numeric for likert scale items
data[1:6] <- as.data.frame(sapply(data[1:6], as.factor)) #factor for other variables
data[2] <- as.data.frame(sapply(data[2], as.numeric)) #age numeric

#semesters of active participation
data$sem_active_participation <- sapply(1:nrow(data), function(i){sum(data[i, 7:28])})

#sem first pt experience
data$sem_first_pt <- sapply(1:nrow(data), function(i){which(as.numeric(data[i, 7:28]) == 1)[1] + 13 - as.numeric(data[i, 3])})
data$sem_first_pt_abs <- as.numeric(sapply(1:nrow(data), function(i){which(as.numeric(data[i, 7:28]) == 1)[1]})) #overall, for cohort comparison

#calculate age at first pt experience
data$age_first_pt <- sapply(1:nrow(data), function(i){data[i, 2] - ((21 - which(as.numeric(data[i, 7:28]) == 1))[1]/2)})

#assign cases to groups (1=ana, 2=ptp, 3=nep, 4=mul)-----
data$group <- sapply(1:nrow(data), function(i){
  if (data[i, 4]==1 & data[i, 5]==0 & data[i, 6]==0){
    1
  } else if (data[i, 4]==0 & data[i, 5]==1 & data[i, 6]==0){
    2
  } else if (data[i, 4]==0 & data[i, 5]==0 & data[i, 6]==1){
    3
  } else
    4
})

#missing values
missing = apply(data[29:61], 1, function(data){sum(is.na(data))/33 * 100}) #missing items
table(missing)
data = subset(data, missing <= 10) #exclude if >= 10% of items are missing
missing_col = apply(data[29:61], 2, function(data){sum(is.na(data[29:61]))/33 * 100})
table(missing_col)
rm(missing, missing_col)

#data evaluation----
questions <- questions[order(questions$sorting),]
questions$median <- apply(data[29:61], 2, median, na.rm=TRUE)
questions$kurtosis <- apply(data[29:61], 2, kurtosis, na.rm=TRUE)
questions$skewness <- apply(data[29:61], 2, skewness, na.rm=TRUE)
questions$lowerquantile <- apply(data[29:61], 2, quantile, na.rm=TRUE)[2, ]
questions$upperquantile <- apply(data[29:61], 2, quantile, na.rm=TRUE)[4, ]
questions$max <- apply(data[29:61], 2, max, na.rm=TRUE)
questions$min <- apply(data[29:61], 2, min, na.rm=TRUE)

#plot
plot(melt(data[29:61])$variable, melt(data[29:61])$value)

#group analyses (1=ana, 2=comed, 3=physiol, 4=multiple combinations)-----

#compare age
summary(aov(data$age ~ data$group)) #anova
if (summary(aov(data$age ~ data$group))[[1]][[1,"Pr(>F)"]] <=.05){
  TukeyHSD(aov(data$age ~ data$group)) #post-hoc tukey
  TukeyHSD(aov(data$age ~ data$group))[[1]][[1,"p adj"]] #extract p value
} else {
  print("Difference in age between groups is not significant")
}

#compare duration of active participation
summary(aov(data$sem_active_participation ~ data$group))
if (summary(aov(data$sem_active_participation ~ data$group))[[1]][[1,"Pr(>F)"]] <= .05){
  TukeyHSD(aov(data$sem_active_participation ~ data$group))
  pairwise.t.test(data$sem_active_participation, data$group, p.adj = "bonf")
} else {
  print("Duration of active participation does not differ significanntly among groups")
}

#first pt comparison (age)
summary(aov(data$age_first_pt ~ data$group))
if (summary(aov(data$age_first_pt ~ data$group))[[1]][[1, "Pr(>F)"]] <= .05){
  TukeyHSD(aov(data$age_first_pt ~ data$group))
} else {
  print("No significant difference in the age of first pt experience")
}

#first pt comparison (semester)
summary(aov(data$sem_first_pt ~ data$group))
if (summary(aov(data$sem_first_pt ~ data$group))[[1]][[1, "Pr(>F)"]] <= .05){
  TukeyHSD(aov(data$sem_first_pt ~ data$group))
  pairwise.t.test(data$sem_first_pt, data$group, p.adj = "bonf")
} else {
  print("No significant difference in the semester of first pt experience")
}

#compare question results between groups (kruskal-wallis-test) with dunn's post-hoc if p < .05
no_mult <- subset(data, group != 4)
questions$pkw <- sapply(29:61, function(i){round(kruskal.test(no_mult[, i] ~ no_mult$group)$p.value, 4)}) #p values from kruskal wallis test
questions$med_ana <- sapply(29:61, function(i){median(no_mult[no_mult$group == 1, i], na.rm = T)})
questions$iqr_ana <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 1, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 1, i], na.rm = T)[4], sep=";")})
questions$med_ptp <- sapply(29:61, function(i){median(no_mult[no_mult$group == 2, i], na.rm = T)})
questions$iqr_ptp <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 2, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 2, i], na.rm = T)[4], sep=";")})
questions$med_nep <- sapply(29:61, function(i){median(no_mult[no_mult$group == 3, i], na.rm = T)})
questions$iqr_nep <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 3, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 3, i], na.rm = T)[4], sep=";")})
questions$ana_ptp_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][1], 5)}) #1 vs 2
questions$ana_nep_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][2], 5)}) #1 vs 3
questions$ptp_nep_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][4], 5)}) #2 vs 3
rm(no_mult)

#effect size (cliff's delta, doi:10.1037/0033-2909.114.3.494)
i <- 29
effsize1 <- c()
sup_ci1 <- c()
inf_ci1 <- c()
effsize2 <- c()
sup_ci2 <- c()
inf_ci2 <- c()
effsize3 <- c()
sup_ci3 <- c()
inf_ci3 <- c()
while(i <= 61){
  effsize1 <- c(effsize1, round(cliff.delta(data[data$group==1, i], data[data$group==2, i])[[1]], 3))
  sup_ci1 <- c(sup_ci1, round(cliff.delta(data[data$group==1, i], data[data$group==2, i])[[2]][[1]], 3))
  inf_ci1 <- c(inf_ci1, round(cliff.delta(data[data$group==1, i], data[data$group==2, i])[[2]][[2]], 3))
  effsize2 <- c(effsize2, round(cliff.delta(data[data$group==1, i], data[data$group==3, i])[[1]], 3))
  sup_ci2 <- c(sup_ci2, round(cliff.delta(data[data$group==1, i], data[data$group==3, i])[[2]][[1]], 3))
  inf_ci2 <- c(inf_ci2, round(cliff.delta(data[data$group==1, i], data[data$group==3, i])[[2]][[2]], 3))
  effsize3 <- c(effsize3, round(cliff.delta(data[data$group==2, i], data[data$group==3, i])[[1]], 3))
  sup_ci3 <- c(sup_ci3, round(cliff.delta(data[data$group==2, i], data[data$group==3, i])[[2]][[1]], 3))
  inf_ci3 <- c(inf_ci3, round(cliff.delta(data[data$group==2, i], data[data$group==3, i])[[2]][[2]], 3))
  i <- i+1
}

questions <- cbind(questions,
                   effsize1, sup_ci1, inf_ci1,
                   effsize2, sup_ci2, inf_ci2,
                   effsize3, sup_ci3, inf_ci3
                   )
names(questions)[23:31] <- c("effsize_ana_vs_ptp", "effsize_ana_vs_ptp_sup_ci", "effsize_ana_vs_ptp_inf_ci",
                             "effsize_ana_vs_nep", "effsize_ana_vs_nep_sup_ci", "effsize_ana_vs_nep_inf_ci",
                             "effsize_ptp_vs_nep", "effsize_ptp_vs_nep_sup_ci", "effsize_ptp_vs_nep_inf_ci")



#partial correlation
data_complete_naomit <- na.omit(data[, c(29:61, 66)])
data_complete_naomit[, 34] <- as.numeric(data_complete_naomit[, 34])
correlmatrix <- cbind(
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 1], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 2], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 3], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 4], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 5], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})
)
rownames(correlmatrix) <- 6:33
colnames(correlmatrix) <- 1:5
pmatrix <- cbind(
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 1], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[2]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 2], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[2]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 3], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[2]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 4], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[2]]}),
  sapply(6:33, function(i){pcor.test(data_complete_naomit[, 5], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[2]]})
)
rownames(pmatrix) <- 6:33
colnames(pmatrix) <- 1:5
heatmap.2(correlmatrix, dendrogram = "row", tracecol = NA, col=redgreen(1000),
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="complete"),
          Colv = F)

#agglomerated data evaluation-----
#calculating intra-individual medians
data_complete_naomit$mot <- apply(data_complete_naomit[1:5], 1, median)
data_complete_naomit$suppothers <- apply(data_complete_naomit[6:10], 1, median)
data_complete_naomit$selfimpr <- apply(data_complete_naomit[11:21], 1, median)
data_complete_naomit$fin <- apply(data_complete_naomit[22], 1, median)
data_complete_naomit$feed <- apply(data_complete_naomit[23:29], 1, median)
data_complete_naomit$posrew <- apply(data_complete_naomit[6:29], 1, median)
data_complete_naomit$negrew <- apply(data_complete_naomit[30:33], 1, median)

#motivation vs (positive) rewards
median(data_complete_naomit$mot)
quantile(data_complete_naomit$mot)
median(data_complete_naomit$posrew)
quantile(data_complete_naomit$posrew)
wilcox.test(data_complete_naomit$mot, data_complete_naomit$posrew, paired=T)
wilcoxsign_test(data_complete_naomit$mot ~ data_complete_naomit$posrew, paired=T)
-2.9235/sqrt(218*2) #effect size (z/sqrt(n[a]+n[b]))

#comaprison of rewards
median(data_complete_naomit$posrew)
quantile(data_complete_naomit$posrew)
median(data_complete_naomit$negrew)
quantile(data_complete_naomit$negrew)
wilcox.test(data_complete_naomit$posrew, data_complete_naomit$negrew, paired=T)
wilcoxsign_test(data_complete_naomit$posrew ~ data_complete_naomit$negrew, paired=T)
-12.739/sqrt(218*2) #effect size (z/sqrt(n[a]+n[b]))

#suppothers vs selfimprov
median(data_complete_naomit$suppothers)
quantile(data_complete_naomit$suppothers)
median(data_complete_naomit$selfimpr)
quantile(data_complete_naomit$selfimpr)
wilcox.test(data_complete_naomit$suppothers, data_complete_naomit$selfimpr, paired=T)
wilcoxsign_test(data_complete_naomit$suppothers ~ data_complete_naomit$selfimpr, paired=T)
-11.108/sqrt(218*2) #effect size (z/sqrt(n[a]+n[b]))

#plot data
long <- melt(data_complete_naomit[35:41])
long$variable <- as.factor(long$variable)
boxplot(long$value ~ long$variable)
points(long$value ~ long$variable)

library(ggplot2)
ggplot(long, aes(x=variable, y=value)) + 
  geom_jitter(alpha = .5, color = "black", size = .7, stroke = 0) +
  geom_boxplot(
    color = "black", fill = "grey", outlier.shape = NA,
    size = .3, alpha = .8, width = .7, position=position_dodge(1)) +
  stat_summary(fun.y=median, geom="line", size=2, color = "black") +
  theme_classic() +
  theme(legend.title=element_text(color="black", size = 5)) + 
  labs(y = "Likert scale value", x = "Question", size = 5 ) +
  theme(panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust = 0.5,size = 7, color="black"),
        axis.text.y = element_text(size = 7, hjust = 1,vjust = .4, color="black"),
        axis.ticks = element_line(color = "black", size = .3),
        axis.line = element_line(size = .3),
        aspect.ratio = 1) +
  scale_x_discrete(labels= labels)
rm(long)

#PCA
pcamatrix <- read.csv("nbt_pcamatrix.csv", sep=";") #load pca matrix with cliff's deltas and re-coded p values from kruskal-wallis testing
rownames(pcamatrix) <- pcamatrix[, 1]
pca <- prcomp(pcamatrix[2:7])

summary(pca)
fviz_pca_ind(pca,
             axes = c(1, 2),
             #col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
) +
  scale_x_continuous(breaks=c(seq(-1, 1, 0.5)), limits=c(-1, 1)) +
  scale_y_continuous(breaks=c(seq(-1, 1, 0.5)), limits=c(-1,1)) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

fviz_eig(pca)

res.pca <- prcomp(data[2:4])
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)