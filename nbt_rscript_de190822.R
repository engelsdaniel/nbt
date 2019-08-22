#load data-----
options(scipen = 999)
rm(list=ls())
graphics.off()
setwd("c:\\users\\de\\dropbox\\lehre\\nbt_manuscript\\nbt_dataanalysis")
data <- read.csv("nbt_rawdata.csv", header=T, sep=";")[2:62]
questions <- read.csv("nbt_manuscript_questions.csv", header = T, sep = ";")
library(ggplot2)
library(reshape)
library(svglite)
library(multcomp)
library(e1071)
library(PMCMR)
library(effsize)
library(ppcor)
library(MASS)

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

#data evaluation----
questions_analysis <- questions[order(questions$sorting),]
questions_analysis$median <- apply(data[29:61], 2, median, na.rm=TRUE)
questions_analysis$kurtosis <- apply(data[29:61], 2, kurtosis, na.rm=TRUE)
questions_analysis$skewness <- apply(data[29:61], 2, skewness, na.rm=TRUE)
questions_analysis$lowerquantile <- apply(data[29:61], 2, quantile, na.rm=TRUE)[2, ]
questions_analysis$upperquantile <- apply(data[29:61], 2, quantile, na.rm=TRUE)[4, ]
questions_analysis$max <- apply(data[29:61], 2, max, na.rm=TRUE)
questions_analysis$min <- apply(data[29:61], 2, min, na.rm=TRUE)

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
quest_analysis$pkw <- sapply(29:61, function(i){round(kruskal.test(no_mult[, i] ~ no_mult$group)$p.value, 4)}) #p values from kruskal wallis test
quest_analysis$med_ana <- sapply(29:61, function(i){median(no_mult[no_mult$group == 1, i], na.rm = T)})
quest_analysis$iqr_ana <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 1, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 1, i], na.rm = T)[4], sep=";")})
quest_analysis$med_ptp <- sapply(29:61, function(i){median(no_mult[no_mult$group == 2, i], na.rm = T)})
quest_analysis$iqr_ptp <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 2, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 2, i], na.rm = T)[4], sep=";")})
quest_analysis$med_nep <- sapply(29:61, function(i){median(no_mult[no_mult$group == 3, i], na.rm = T)})
quest_analysis$iqr_nep <- sapply(29:61, function(i){paste(quantile(no_mult[no_mult$group == 3, i], na.rm = T)[2], quantile(no_mult[no_mult$group == 3, i], na.rm = T)[4], sep=";")})
quest_analysis$ana_ptp_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][1], 5)}) #1 vs 2
quest_analysis$ana_nep_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][2], 5)}) #1 vs 3
quest_analysis$ptp_nep_posthoc <- sapply(29:61, function(i){round(posthoc.kruskal.dunn.test(x=no_mult[, i], g=no_mult$group, p.adjust.method="bonf")[[3]][4], 5)}) #2 vs 3
rm(no_mult)

#effect size (cliff's delta, doi:10.1037/0033-2909.114.3.494)
effectsize <- as.data.frame(matrix(rep(NA, 3*33), nrow=33, ncol=3))
row.names(effectsize) <- quest_names
colnames(effectsize) <- c("delta_ana_vs_ptp", "delta_ana_vs_nep", "delta_ptp_vs_nep")
i <- 29
while (i <= 61){
  if(quest_analysis[i-28, "ana_vs_ptp"] <= .05){
    effectsize[i-28, 1] <- round(as.numeric(cliff.delta(data[data$group==1, i], data[data$group==2, i], return.dm=TRUE)[1]), 3)
  } else {
    effectsize[i-28, 1] <- "NA"
  }
  if(quest_analysis[i-28, "ana_vs_nep"] <= .05){
    effectsize[i-28, 2] <- round(as.numeric(cliff.delta(data[data$group==1, i], data[data$group==3, i],return.dm=TRUE)[1]), 3)
  } else {
    effectsize[i-28, 2] <- "NA"
  }
  if(quest_analysis[i-28, "ptp_vs_nep"] <= .05){
    effectsize[i-28, 3] <- round(as.numeric(cliff.delta(data[data$group==2, i], data[data$group==3, i],return.dm=TRUE)[1]), 3)
  } else {
    effectsize[i-28, 3] <- "NA"
  }
  i <- i+1
}
quest_analysis <- cbind(quest_analysis, effectsize)
rm(effectsize, i)

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
effectsize <- matrix(rep(NA, 3*33), nrow = 33, ncol = 3)
effectsize[, 1] <- paste(effsize1, paste("[", sup_ci1, ";", inf_ci1, "]", sep = ""), sep = " ")
effectsize[, 2] <- paste(effsize2, paste("[", sup_ci2, ";", inf_ci2, "]", sep = ""), sep = " ")
effectsize[, 3] <- paste(effsize3, paste("[", sup_ci3, ";", inf_ci3, "]", sep = ""), sep = " ")

quest_analysis <- cbind(quest_analysis, effectsize)
names(quest_analysis)[24:26] <- c("ana_vs_ptp", "ana_vs_nep", "ptp_vs_nep")

#partial correlation
data_complete_naomit <- na.omit(data[, c(29:61, 66)])
data_complete_naomit[, 34] <- as.numeric(data_complete_naomit[, 34])

q01 <- sapply(6:33, function(i){pcor.test(data_complete_naomit[, 1], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})
q02 <- sapply(6:33, function(i){pcor.test(data_complete_naomit[, 2], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})
q03 <- sapply(6:33, function(i){pcor.test(data_complete_naomit[, 3], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})
q04 <- sapply(6:33, function(i){pcor.test(data_complete_naomit[, 4], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})
q05 <- sapply(6:33, function(i){pcor.test(data_complete_naomit[, 5], data_complete_naomit[, i], data_complete_naomit[, 34], method = "spearman")[[1]]})

correlmatrix <- cbind(q01, q02, q03, q04, q05)

heatmap.2(correlmatrix, dendrogram = "row", tracecol = NA, col=redgreen(1000),
          distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="complete"))

#agglomerated data evaluation-----
library(coin)
library(multcomp)
#calculating intra-individual medians
data_complete_naomit$intr <- apply(data_complete_naomit[1:5], 1, median)
data_complete_naomit$teach <- apply(data_complete_naomit[6:10], 1, median)
data_complete_naomit$selfdir <- apply(data_complete_naomit[11:21], 1, median)
data_complete_naomit$frew <- apply(data_complete_naomit[22], 1, median)
data_complete_naomit$vrew <- apply(data_complete_naomit[23:29], 1, median)
data_complete_naomit$nrew <- apply(data_complete_naomit[30:33], 1, median)
data_complete_naomit$id <- seq.int(nrow(data_complete_naomit))
data_complete_naomit$mot <- apply(data_complete_naomit[1:21], 1, median)
data_complete_naomit$rew <- apply(data_complete_naomit[22:29], 1, median)
data_complete_naomit$extr <- apply(data_complete_naomit[6:21], 1, median)
#motivation vs rewards
median(data_complete_naomit$mot)
quantile(data_complete_naomit$mot)
median(data_complete_naomit$rew)
quantile(data_complete_naomit$rew)
wilcox.test(data_complete_naomit$mot, data_complete_naomit$rew, paired=T)
wilcoxsign_test(data_complete_naomit$mot ~ data_complete_naomit$rew, paired=T)
-6.9146/sqrt(218*2) #effect size (z/sqrt(n[a]+n[b]))
cliff.delta(data_complete_naomit$mot, data_complete_naomit$rew, paired=T)
#intrinsic vs extrinsic motivation
median(data_complete_naomit$intr)
quantile(data_complete_naomit$intr)
median(data_complete_naomit$extr)
quantile(data_complete_naomit$extr)
wilcox.test(data_complete_naomit$intr, data_complete_naomit$extr, paired=T)
cliff.delta(data_complete_naomit$intr, data_complete_naomit$extr, paired=T)
#intrinsic vs extrinsic: teach and selfdir
long_intrteachselfdir <- melt(data_complete_naomit[, c("intr", "teach", "selfdir")])
pairwise.wilcox.test(long_intrteachselfdir$value, long_intrteachselfdir$variable, paired=T, p.adj="bonferroni")
cliff.delta(data_complete_naomit$intr, data_complete_naomit$teach, paired=T)
cliff.delta(data_complete_naomit$intr, data_complete_naomit$selfdir, paired=T)
#extrinsic motivation: teach vs selfdir
median(data_complete_naomit$teach)
quantile(data_complete_naomit$teach, type=1)
median(data_complete_naomit$selfdir)
quantile(data_complete_naomit$selfdir, type=1)
wilcox.test(data_complete_naomit$teach, data_complete_naomit$selfdir, paired=T)
cliff.delta(data_complete_naomit$teach, data_complete_naomit$selfdir, paired=T)
#extrinsic rewards: verbals vs financial
median(data_complete_naomit$frew)
quantile(data_complete_naomit$frew)
median(data_complete_naomit$vrew)
quantile(data_complete_naomit$vrew)
wilcox.test(data_complete_naomit$frew, data_complete_naomit$vrew, paired=T)
cliff.delta(data_complete_naomit$frew, data_complete_naomit$vrew, paired=T)
#extrinsic rewards: negative vs verbald and financial
median(data_complete_naomit$rew)
quantile(data_complete_naomit$rew)
median(data_complete_naomit$nrew)
quantile(data_complete_naomit$nrew)
long_rewcomp <- melt(data_complete_naomit[, c("frew", "vrew", "nrew")])
pairwise.wilcox.test(long_rewcomp$value, long_rewcomp$variable, paired=T, p.adj="bonferroni")
cliff.delta(data_complete_naomit$nrew, data_complete_naomit$vrew, paired=T)
cliff.delta(data_complete_naomit$nrew, data_complete_naomit$frew, paired=T)

#plot data
library(reshape)
long <- melt(data_complete_naomit[34:41], id.vars = c("id", "group"))
long$cat <- apply(long, 1, function(x) {
  ifelse(long[3] == "intr", "intr",
         ifelse(long[3] == "teach", "extr",
                ifelse(long[3] == "selfdir", "extr",
                       ifelse(long[3] == "nrew", "nrew", "rew"))))})
library(ggplot2)
boxplot <- ggplot(long, aes(x=variable, y=value)) + 
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
print(boxplot)