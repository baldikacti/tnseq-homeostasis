metaData = read.table("/Users/tingtingzhao/Documents/Research/TnqData/FinalData/RESULTS_FDR_LASTLEVELONLY.TAB",
                      header=TRUE)
head(metaData)

uniqueValues = unique(metaData$Annotation_TC)
print(uniqueValues)

# Take a subset of data that only has annotation "Conditionally_Beneficial" 
## and "Conditionally_Detrimental"
newData = metaData[ which(metaData$Annotation_TC=="Conditionally_Beneficial"
                         | metaData$Annotation_TC =="Conditionally_Detrimental" |
                           metaData$Annotation_TC=="Conditionally_Essential"), ]
dim(newData)
#19674   14
################################################################################
## get all the gene names from this subset of data "newData"
genes = newData$locus_tag
length(unique(genes)) #1736

################################################################################
## TnqSeq Final version of Data Analysis
dft = read.table("/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC_mid_batchcorrected_dclpb_included.tsv",sep = '\t', header = TRUE)
dftNonZero = dft
head(dftNonZero)
## in dfUNonZero, select only the genes that show up in "genes=newData$locus_tag"
dim(dftNonZero)
# 873976      9
source("~/Documents/Research/TnqData/newdata/code/create_Yknockoff_para.R")
library(tidyverse)
dft = dftNonZero%>%select(-batch,-N)%>%group_by(locus_tag,strain,condition,slevel,rep,time)
dftSubset = dft %>% filter(locus_tag %in% newData$locus_tag)
dftSubset
dftSubset%>%summarise(n())
# 369768 × 7
dft1<-dftSubset%>%summarise(mean_tncnt=mean(tncnt))
dft2<-dft1%>%pivot_wider(names_from = locus_tag,values_from = mean_tncnt)
dft2
# 213*1741

df = dft2
genesT = as.data.frame(lapply(dft2[, 6:dim(dft2)[2]], unlist))
dim(genesT)
# 213*1736 # 1736 unique number of genes

# zeroColsU =  which(apply(genesT, 2, function(c)(sd(c)==0))==TRUE) # integer(0)
genesUNoZero = genesT

##############################################################################
##############################################################################
## Commenting out the following piece of codes as we only consider "conditional 
## beneficial" and "conditional detrimental" genes

## find the columns that at time zero, condition is none and slevel is none
## the average of the count of a gene is smaller than 3 across 3 replications.
## if a gene satisfies this condition, it is considered as an essential gene.
## find the indices of the column that satisfies this condition
#slevelNone = which(df$slevel=="none")
#conditionNone = which(df$condition=="none")
#time0 = which(df$time==0)
#condNoneSlevelNoneTime0 = intersect(intersect(slevelNone,conditionNone), time0)
#baseU = dft2[condNoneSlevelNoneTime0, ]
#baseUMtx = as.data.frame(lapply(baseU[, 6:dim(dfu)[2]], unlist))
#meanBaseUMtx = apply(baseUMtx, 2, mean)
#essential = which(meanBaseUMtx<3)
#numEssential = length(essential)
#removeindex=essential+5
## 729

###############################################################################
## start the analysis
## first subset observations with condition is none type and slevel is none
condNoneSlevelNoneInd = intersect(which(df$slevel=="none"), which(df$condition=="none"))
timeZeroInd = which(df$time==0)
condSlevelNoneTime0 = intersect(timeZeroInd, condNoneSlevelNoneInd)
basedfu = df[condSlevelNoneTime0,]
basedfu # 16*579
# basedfu2 <- basedfu%>%group_by(strain,condition,slevel)%>%summarise_all(mean)%>%filter(strain!='DCLPB')
basedfu2 <- basedfu%>%group_by(strain,condition,slevel)%>%summarise_all(mean)#%>%filter(strain!='DCLPB')
basedfu2$rep=4
library(dplyr)
basedfu3 = basedfu2 %>% slice(rep(1:n(), each = 4))
basedfu3$rep = rep(c(1,2,3,4), 5)
# basedfu3 20*579 20 comes from 4 replications multiply 4 strains such as wild, dlon, dclpa and etc
basedfufill<-rbind(basedfu,basedfu2)%>%arrange(strain,condition,slevel,rep)
basedfufill
dfuNoTime0 = df[-condSlevelNoneTime0, ]
dfuNoTime0 # 197*1741

##############################################################################################
### since the base data only have three replications, we change the base level as the average 
### across the three replications
# dat1 = basedfufill%>%pivot_longer(col=-(1:5),names_to = 'gene',values_to = 'basecount')%>%mutate(Type='base',count=basecount)
dat1 = basedfufill%>%pivot_longer(col=-(1:5),names_to = 'gene',values_to = 'basecount')%>%mutate(Type='base',count=basecount)
# dat2<-dfuNoTime0%>%select(-3324)%>%pivot_longer(col=-(1:5),names_to = 'gene',values_to = 'count')%>%mutate(Type='Condition')
dat2 = dfuNoTime0%>%pivot_longer(col=-(1:5),names_to = 'gene',values_to = 'count')%>%mutate(Type='Condition')
dat1%>%filter(strain=='wild-type',gene=='CCNA_00001',rep==1)
dat2%>%filter(strain=='wild-type',gene=='CCNA_00001',rep==1)
dat3<-dat1%>%full_join(dat2)%>%group_by(strain,gene,rep)%>%fill(basecount)
dat3%>%filter(strain=='wild-type',gene=='CCNA_00001',rep==1)
dat3 %>% mutate(basecount=as.integer(basecount))
dat3 %>% add_column(basecountVec = unlist(dat3$basecount))
## The outcome will be 
## log(average total count at 24 hours for each strain+ 1/(average total count at time 0 for
## each strain under none condition none slevel +1))
## the plus 1 is for zero adjustment
dat3$basecount = dat3$basecount + 1
dat3$count = dat3$count + 1
dat4<-dat3%>%mutate(ratio=count/basecount)%>%filter(Type=='Condition')%>%select(-Type,-count,-basecount)
dat4$logratio = log(dat4$ratio)
dat5<-dat4%>%pivot_wider(id_cols = 1:5,names_from = gene,values_from = logratio)

# subset for each strain
wildtype= dat5 %>%filter(strain=="wild-type") 
#wildtype = wildtype %>% filter(rep<4)
## find the genetic code for lon, clpa, clpb, dnak and dnaj
## CCNA_02037: lon
## CCNA_00922: clpB
## CCNA_02553: clpA
## CCNA_00010: dnaK
## CCNA_00011: dnaJ
dlon = dat5 %>% filter(strain=="DLON")
#dlon = dlon  %>% select(-"CCNA_02037")
#dlon = dlon %>% filter(rep<4)
dclpa = dat5 %>% filter(strain=="DCLPA")
#dclpa = dclpa %>% select(-"CCNA_02553")
#dclpa = dclpa %>% filter(rep<4)
dnak = dat5 %>% filter(strain=="dnak-dnaJ")
#dnak = dnak %>% select(-c("CCNA_00010", "CCNA_00011"))
#dnak = dnak %>% filter(rep<4)
dclpb = dat5 %>% filter(strain=="DCLPB")
#dclpb = dclpb %>% select(-"CCNA_00922")

wildMtx = as.matrix(wildtype[ , 6:dim(wildtype)[2]])
dlonMtx = as.matrix(dlon[, 6:dim(dlon)[2]])
dclpaMtx = as.matrix(dclpa[, 6:dim(dclpa)[2]])
dnakMtx = as.matrix(dnak[, 6:dim(dnak)[2]])
dclpbMtx = as.matrix(dclpb[, 6:dim(dclpb)[2]])

## Save the new generated data

logRatioKnockoffParamWildAdjusted = create_Yknockoff_para(wildMtx, "sdp", shrink=TRUE)
save(logRatioKnockoffParamWildAdjusted, file="~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamWildWholeAdjusted.Rdata")
#logRatioKnockoffParamWildAdjusted = create_Yknockoff_para(wildMtx, "sdp", shrink=TRUE)
#save(logRatioKnockoffParamWildAdjusted, file="~/Documents/TnqData/newdata/code/logRatioKnockoffParamWildWholeAdjusted4Rep.Rdata")
logRatioKnockoffParamDlonAdjusted = create_Yknockoff_para(dlonMtx, "sdp", shrink=TRUE)
save(logRatioKnockoffParamDlonAdjusted, file="~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamDlonWholeAdjusted.Rdata")
#save(logRatioKnockoffParamDlonAdjusted, file="~/Documents/TnqData/newdata/code/logRatioKnockoffParamDlonWholeAdjusted4Rep.Rdata")
logRatioKnockoffParamDclpaAdjusted = create_Yknockoff_para(dclpaMtx, "sdp", shrink=TRUE)
save(logRatioKnockoffParamDclpaAdjusted, file="~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamDclpaWholeAdjusted.Rdata")
#save(logRatioKnockoffParamDclpaAdjusted, file="~/Documents/TnqData/newdata/code/logRatioKnockoffParamDclpaWholeAdjusted4Rep.Rdata")
logRatioKnockoffParamdnakAdjusted = create_Yknockoff_para(dnakMtx, "sdp", shrink=TRUE)
save(logRatioKnockoffParamdnakAdjusted, file="~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamdnakWholeAdjusted.Rdata")
#save(logRatioKnockoffParamdnakAdjusted, file="~/Documents/TnqData/newdata/code/logRatioKnockoffParamdnakWholeAdjusted4Rep.Rdata")
logRatioKnockoffParamdclpbAdjusted = create_Yknockoff_para(dclpbMtx, "sdp", shrink=TRUE)
save(logRatioKnockoffParamdclpbAdjusted, file="~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamdclpbWholeAdjusted.Rdata")
#save(logRatioKnockoffParamdclpbAdjusted, file="~/Documents/TnqData/newdata/code/logRatioKnockoffParamdclpbWholeAdjusted4Rep.Rdata")

load("~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamWildWholeAdjusted.Rdata")
load("~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamDlonWholeAdjusted.Rdata")
load("~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamDclpaWholeAdjusted.Rdata")
load("~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamdnakWholeAdjusted.Rdata")
load("~/Documents/Research/TnqData/FinalData/UC/logRatioKnockoffParamdclpbWholeAdjusted.Rdata")


logRatioKnockoffWild = generateMultiKnockoff(wildMtx, 
                                             logRatioKnockoffParamWildAdjusted$mu_k, 
                                             logRatioKnockoffParamWildAdjusted$Sigma_k,n=100)
source("/Users/tingtingzhao/Documents/Research/YKnock/Deep-YKnock_Tingting/codes/stat_modelY_classification_coef.R")
library(nnet)
wildtype$condition = relevel(factor(wildtype$condition), ref="none")
dlon$condition = relevel(factor(dlon$condition), ref="none")
dclpa$condition = relevel(factor(dclpa$condition), ref="none")
dnak$condition = relevel(factor(dnak$condition), ref="none")
dclpb$condition = relevel(factor(dclpb$condition), ref="none")

resultInd = vector(mode = "list", length = 100)
resultNames =  vector(mode = "list", length = 100)

selectVar = function(subdf, resdf, knockoffdf){
  for(i in 1:100){
    print(i)
    Yfeatures = cbind(subdf, knockoffdf[[i]])
    Z = stat_modelY_classification_coef(resdf$condition, Yfeatures)
    Z = abs(Z)
    Z= Z/max(Z)
    r = ncol(Yfeatures)/2
    orig = 1:r
    W = abs(Z[orig]) - abs(Z[orig+r])
    plot(W)
    tau_max = knockoff.threshold(W, fdr = 0.1, offset = 0)
    S_max = (which(W>tau_max))
    varNames = colnames(wildMtx)[S_max]
    resultInd[[i]] = S_max
    resultNames[[i]] = varNames
  }
  result = list(resultInd=resultInd, resultNames=resultNames)
  return(result)
}

selectVarFreq = function(result){
  countInd = sort(table(unlist(lapply(result$resultInd, unique))), decreasing=TRUE)
  countNames = sort(table(unlist(lapply(result$resultNames, unique))), decreasing=TRUE)
  selectGenes = names(countNames)[as.numeric(countNames)>20]
  selectInd = countNames[as.numeric(countNames)>20]
  output = list(selectGenes=selectGenes, selectInd=selectInd)
  return(output)
}

wildResult = selectVar(wildMtx, wildtype,logRatioKnockoffWild)
wildSelect = selectVarFreq(wildResult)
wildSelect$selectGenes

selectDF = cbind(wildtype[, 1:5], wildtype[,names(wildSelect$selectInd)])
#wildHigh = selectDF %>% filter(strain=="wild-type" && slevel=="HIGH")
wildHigh = selectDF %>% filter(strain=="wild-type"&condition!="none")
wildHigh = wildHigh %>%  unite("Pert", condition:slevel, remove = FALSE)
wildhigh_long = wildHigh%>%pivot_longer(-(1:6))

library(hrbrthemes)
library(forcats)
library(glue)
library(rlang)
library(ggplot2)
library(viridis)
library(colorspace)
library(ggdendro)
library(grid)
#wildhigh_long = wildhigh_long %>% mutate(Pert = fct_relevel(Pert, 
#                                                            "canavanine_HIGH", "canavanine_MEDIUM", 
#                                                            "canavanine_LOW", "heat_HIGH", "heat_MEDIUM", "heat_LOW",
#                                                            "oxidative-peroxide_HIGH", "oxidative-peroxide_MEDIUM", 
#                                                            "oxidative-peroxide_LOW"))
wildhigh_long = wildhigh_long %>% mutate(condition = fct_relevel(condition,
                                                                 "canavanine", "heat", "oxidative-peroxide"))
wildhigh_long = wildhigh_long %>% mutate(slevel = fct_relevel(slevel,
                                                              "LOW", "MEDIUM", "HIGH"))
wildhigh_long = wildhigh_long %>% mutate(name = fct_relevel(name, (wildSelect$selectGenes))) 
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/wildhigh_long.Rata',wildhigh_long)
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/wildSelect.Rata', wildSelect)
wild.avg = wildhigh_long%>%group_by(strain,Pert,condition,slevel,time,name)%>%summarise(mean_value=mean(value))


#-----
load(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/wildhigh_long.Rata')
condition_names <- c(
  'canavanine' = "CANAVANINE",
  'heat' = "HEAT",
  'oxidative-peroxide' = "OXIDATIVE")
g = ggplot(wildhigh_long, aes(slevel, name, fill=value))+
  geom_tile()+
  scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Blue-Red 2")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"))
#panel.margin=unit(0,"pt"))

#------
# Run clustering
wildCluster <- t(wildtype[,names(wildSelect$selectInd)])
wild.dendro <- as.dendrogram(hclust(d = dist(x = wildCluster)))
# Create dendro
dendro.plot <- ggdendrogram(data = wild.dendro, rotate = TRUE)+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_text(size=12, face="bold"))
# Preview the plot
print(dendro.plot)
## putting heatmap and dendrogram together
grid.newpage()
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
### Re-order heatmap rows to match dendrogram
wild.order <- order.dendrogram(wild.dendro)
# Order the levels according to their position in the cluster
wild.avg$name = factor(x=wild.avg$name, 
                       levels = colnames(wildtype[, names(wildSelect$selectInd)])[wild.order],
                       ordered = TRUE)
## recreate the heatmap and reposition the legend
g = ggplot(wild.avg, aes(slevel, name, fill=mean_value))+
  geom_tile()+
  #scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Red-Green")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        #axis.text.y = element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "top")
grid.newpage()
# print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
## Align dendrogram tips with heatmap rows
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 0.95))
ggsave(file="~/Documents/Research/TnqData/FinalData/UC/wild2.pdf", width=10, height=6)
#----
write.csv(selectDF, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1LogRatioWholeAdjusted.csv")
## write table of countNames to csv
write.csv(wildSelect$selectInd, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1DFUFreqWholeAdjusted.csv")
############################################################################################
dfNetWild = data.frame(source=wildSelect$selectGenes, destination=rep("wild", length(wildSelect$selectGenes)))

#This is the old reference file,
#lotus = read.csv("~/Documents/Research/TnqData/newdata/locus_attribute.tab",sep = '\t', header = TRUE)
# I am going to use the new one shared by Berent on Aug 16th, 2022.
lotus = read.csv("/Users/tingtingzhao/Documents/Research/TnqData/newdata/code/Final/CCNA_ref.csv", header = TRUE)
lotus = lotus[, 1:5]
# get genes
wildLotus = lotus[which(lotus$locus_tag  %in% wildSelect$selectGenes), ]
write.csv(wildLotus, file="~/Documents/Research/TnqData/FinalData/UC/wildLotus01WholeAdjusted.csv")

library(glmnet)
res=factor(wildtype$condition)
wildDF = data.frame(res, wildMtx[, names(wildSelect$selectInd)])
class(wildDF)
colnames(wildDF)
glmnet.fit = multinom(res~.,wildDF)
wildLotus = lotus[which(lotus$locus_tag  %in% wildSelect$selectGenes), ]
write.csv(wildLotus, file="~/Documents/Research/TnqData/FinalData/UC/wildLotus01WholeAdjusted4Rep.csv")
write.csv(t(summary(glmnet.fit)$coefficients), file="~/Documents/Research/TnqData/FinalData/UC/wildcoefficient01WholeAdjusted4Rep.csv")
################################################################################
## remove lon from original and knockoff data
logRatioKnockoffDlon = generateMultiKnockoff(dlonMtx, 
                                             logRatioKnockoffParamDlonAdjusted$mu_k, 
                                             logRatioKnockoffParamDlonAdjusted$Sigma_k,n=100)
dlonResult = selectVar(dlonMtx, dlon,logRatioKnockoffDlon)
dlonSelect = selectVarFreq(dlonResult)
res=factor(dlon$condition)
dlonDF = data.frame(res, dlonMtx[, names(dlonSelect$selectInd)])
dlonFit = multinom(res~.,dlonDF)
summary(dlonFit) 
dlonLotus = lotus[which(lotus$locus_tag  %in% dlonSelect$selectGenes), ]
summary(dlonFit) 
dlonLotus = lotus[which(lotus$locus_tag  %in% dlonSelect$selectGenes), ]
# write.csv(dlonLotus, file="~/Documents/TnqData/newdata/dlonLotus.csv") corresponds to FDR=0.2
write.csv(dlonLotus, file="~/Documents/Research/TnqData/FinalData/UC/dlonLotus01WholeAdjusted4Rep.csv") # FDR=0.1
# write.csv(t(summary(dlonFit)$coefficients), file="~/Documents/TnqData/newdata/dyloncoefficient01.csv")
write.csv(t(summary(dlonFit)$coefficients), file="~/Documents/Research/TnqData/FinalData/UC/dloncoefficient01WholeAdjusted4Rep.csv")
dfNetDlon = data.frame(source=dlonSelect$selectGenes, destination=rep("dlon", length(dlonSelect$selectGenes)))

selectDF = cbind(dlon[, 1:5], dlon[,names(dlonSelect$selectInd)])
#wildHigh = selectDF %>% filter(strain=="wild-type" && slevel=="HIGH")
dlonHigh = selectDF %>% filter(strain=="DLON"&condition!="none")
dlonHigh = dlonHigh %>%  unite("Pert", condition:slevel, remove = FALSE)
dlonhigh_long = dlonHigh%>%pivot_longer(-(1:6))
dlonhigh_long = dlonhigh_long %>% mutate(condition = fct_relevel(condition,
                                                                 "canavanine", "heat", "oxidative-peroxide"))
dlonhigh_long = dlonhigh_long %>% mutate(slevel = fct_relevel(slevel,
                                                              "LOW", "MEDIUM", "HIGH"))
dlonhigh_long = dlonhigh_long %>% mutate(name = fct_relevel(name, (dlonSelect$selectGenes))) 
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dlonhigh_long.Rata',dlonhigh_long)
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dlonSelect.Rata', dlonSelect)
#-----
load(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dlonhigh_long.Rata')
g = ggplot(dlonhigh_long, aes(slevel, name, fill=value))+
  geom_tile()+
  scale_y_discrete(limits = unique(rev(dlonSelect$selectGenes)),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Blue-Red 2")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"))
#panel.margin=unit(0,"pt"))
#------
ggsave(file="~/Documents/Research/TnqData/newdata/code/Final/dlon2.pdf", width=10, height=6)
#------
# Run clustering
dlonCluster <- t(dlon[,names(dlonSelect$selectInd)])
ind = which(rownames(dlonCluster)=="CCNA_02037")
if(length(ind)>0){
  dlonCluster = dlonCluster[-ind, ]}

dlon.dendro <- as.dendrogram(hclust(d = dist(x = dlonCluster)))
# Create dendro
dendro.plot <- ggdendrogram(data = dlon.dendro, rotate = TRUE)+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_text(size=12, face="bold"))

# Preview the plot
print(dendro.plot)
## putting heatmap and dendrogram together
#grid.newpage()
#print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
### Re-order heatmap rows to match dendrogram
dlon.avg = dlonhigh_long%>%group_by(strain,Pert,condition,slevel,time,name)%>%summarise(mean_value=mean(value))
dlon.order <- order.dendrogram(dlon.dendro)
# Order the levels according to their position in the cluster
dlon.avg$name = factor(x=dlon.avg$name, 
                       levels = colnames(dlon[, names(dlonSelect$selectInd)])[dlon.order],
                       ordered = TRUE)
## recreate the heatmap and reposition the legend
g = ggplot(dlon.avg, aes(slevel, name, fill=mean_value))+
  geom_tile()+
  #scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Red-Green")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        #axis.text.y = element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "top")
grid.newpage()
# print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
## Align dendrogram tips with heatmap rows
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 0.95))

##########################################################################################
#write.csv(dlonDF, file="~/Documents/TnqData/newdata/dataFDR0.1LogRatioWDLONholeAdjusted.csv")
write.csv(dlonDF, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1LogRatioWDLONholeAdjusted4Rep.csv")
## write table of countNames to csv
# write.csv(dlonSelect$selectInd, file="~/Documents/TnqData/newdata/dataFDR0.1DFUFreqDLONWholeAdjusted.csv")
write.csv(dlonSelect$selectInd, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1DFUFreqDLONWholeAdjusted4Rep.csv")

################################################################################
## remove dclpa from dclpaMtx and knockoffs
logRatioKnockoffDclpa = generateMultiKnockoff(dclpaMtx, 
                                              logRatioKnockoffParamDclpaAdjusted$mu_k, 
                                              logRatioKnockoffParamDclpaAdjusted$Sigma_k,n=100)
dclpaResult = selectVar(dclpaMtx, dclpa,logRatioKnockoffDclpa)
dclpaSelect = selectVarFreq(dclpaResult)
res=factor(dclpa$condition)
dclpaDF = data.frame(res, dclpaMtx[, names(dclpaSelect$selectInd)])
dclpaFit = multinom(res~.,dclpaDF)
summary(dclpaFit) 
dclpaLotus = lotus[which(lotus$locus_tag  %in% dclpaSelect$selectGenes), ]
# write.csv(dclpaLotus, file="~/Documents/TnqData/newdata/dclpaLotus01.csv") FDR=0.2
#write.csv(dclpaLotus, file="~/Documents/TnqData/newdata/dclpaLotus01WholeAdjusted.csv") # FDR=0.1
write.csv(dclpaLotus, file="~/Documents/Research/TnqData/FinalData/UC/dclpaLotus01WholeAdjusted4Rep.csv") 
# write.csv(t(summary(dclpaFit)$coefficients), file="~/Documents/TnqData/newdata/dclpacoefficient.csv") #FDR=0.2
# write.csv(t(summary(dclpaFit)$coefficients), file="~/Documents/TnqData/newdata/dclpacoefficient01WholeAdjusted.csv") #FDR=0.1
write.csv(t(summary(dclpaFit)$coefficients), file="~/Documents/Research/TnqData/FinalData/UC/dclpacoefficient01WholeAdjusted4Rep.csv")
intersect(dclpaSelect$selectGenes, dclpaSelect$selectGenes)
dfNetDclpa = data.frame(source=dclpaSelect$selectGenes, destination=rep("dclpa", length(dclpaSelect$selectGenes)))

selectDF = cbind(dclpa[, 1:5], dclpa[,names(dclpaSelect$selectInd)])
#wildHigh = selectDF %>% filter(strain=="wild-type" && slevel=="HIGH")
dclpaHigh = selectDF %>% filter(strain=="DCLPA"&condition!="none")
dclpaHigh = dclpaHigh %>%  unite("Pert", condition:slevel, remove = FALSE)
dclpahigh_long = dclpaHigh%>%pivot_longer(-(1:6))
dclpahigh_long = dclpahigh_long %>% mutate(condition = fct_relevel(condition,
                                                                   "canavanine", "heat", "oxidative-peroxide"))
dclpahigh_long = dclpahigh_long %>% mutate(slevel = fct_relevel(slevel,
                                                                "LOW", "MEDIUM", "HIGH"))
dclpahigh_long = dclpahigh_long %>% mutate(name = fct_relevel(name, (dclpaSelect$selectGenes))) 
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpahigh_long.Rata',dclpahigh_long)
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpaSelect.Rata', dclpaSelect)
#-----
load(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpahigh_long.Rata')
g = ggplot(dclpahigh_long, aes(slevel, name, fill=value))+
  geom_tile()+
  scale_y_discrete(limits = unique(rev(dclpaSelect$selectGenes)),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Blue-Red 2")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        axis.text.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"))
#panel.margin=unit(0,"pt"))
#----
# Run clustering
dclpaCluster <- t(dclpa[,names(dclpaSelect$selectInd)])
dclpa.dendro <- as.dendrogram(hclust(d = dist(x = dclpaCluster)))
# Create dendro
dendro.plot <- ggdendrogram(data = dclpa.dendro, rotate = TRUE)+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_text(size=12, face="bold"))

# Preview the plot
print(dendro.plot)
## putting heatmap and dendrogram together
#grid.newpage()
#print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
### Re-order heatmap rows to match dendrogram
dclpa.avg = dclpahigh_long%>%group_by(strain,Pert,condition,slevel,time,name)%>%summarise(mean_value=mean(value))
dclpa.order <- order.dendrogram(dclpa.dendro)
# Order the levels according to their position in the cluster
dclpa.avg$name = factor(x=dclpa.avg$name, 
                        levels = colnames(dlon[, names(dclpaSelect$selectInd)])[dclpa.order],
                        ordered = TRUE)
## recreate the heatmap and reposition the legend
g = ggplot(dclpa.avg, aes(slevel, name, fill=mean_value))+
  geom_tile()+
  #scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Red-Green")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        #axis.text.y = element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "top")
grid.newpage()
# print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
## Align dendrogram tips with heatmap rows
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 0.95))

################################################################################
write.csv(dclpaDF, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1LogRatiodclpaWholeAdjusted.csv")
## write table of countNames to csv
write.csv(dclpaSelect$selectInd, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1DFUFreqdclpaWholeAdjusted.csv")

################################################################################
logRatioKnockoffdnak = generateMultiKnockoff(dnakMtx, 
                                             logRatioKnockoffParamdnakAdjusted$mu_k, 
                                             logRatioKnockoffParamdnakAdjusted$Sigma_k,n=100)
dnakResult = selectVar(dnakMtx, dnak,logRatioKnockoffdnak)
dnakSelect = selectVarFreq(dnakResult)
res=factor(dnak$condition)
dnakDF = data.frame(res, dnakMtx[, names(dnakSelect$selectInd)])
dnakFit = multinom(res~.,dnakDF)
summary(dnakFit) 
dnakLotus = lotus[which(lotus$locus_tag  %in% dnakSelect$selectGenes), ]
# write.csv(dclpaLotus, file="~/Documents/TnqData/newdata/dclpaLotus01.csv") FDR=0.2
# write.csv(dnakLotus, file="~/Documents/TnqData/newdata/dnakLotus01WholeAdjusted.csv") # FDR=0.1
write.csv(dnakLotus, file="~/Documents/Research/TnqData/FinalData/UC/dnakLotus01WholeAdjusted4Rep.csv")
# write.csv(t(summary(dclpaFit)$coefficients), file="~/Documents/TnqData/newdata/dclpacoefficient.csv") #FDR=0.2
# write.csv(t(summary(dnakFit)$coefficients), file="~/Documents/TnqData/newdata/dnakcoefficient01WholeAdjusted.csv")
write.csv(t(summary(dnakFit)$coefficients), file="~/Documents/Research/TnqData/FinalData/UC/dnakcoefficient01WholeAdjusted4Rep.csv")#FDR=0.1
intersect(dnakSelect$selectGenes, dnakSelect$selectGenes)
dfNetDnak = data.frame(source=dnakSelect$selectGenes, destination=rep("dnak", length(dnakSelect$selectGenes)))

selectDF = cbind(dnak[, 1:5], dnak[,names(dnakSelect$selectInd)])
#wildHigh = selectDF %>% filter(strain=="wild-type" && slevel=="HIGH")
dnakHigh = selectDF %>% filter(strain=="dnak-dnaJ"&condition!="none")
dnakHigh = dnakHigh %>%  unite("Pert", condition:slevel, remove = FALSE)
dnakhigh_long = dnakHigh%>%pivot_longer(-(1:6))
dnakhigh_long = dnakhigh_long %>% mutate(condition = fct_relevel(condition,
                                                                 "canavanine", "heat", "oxidative-peroxide"))
dnakhigh_long = dnakhigh_long %>% mutate(slevel = fct_relevel(slevel,
                                                              "LOW", "MEDIUM", "HIGH"))
dnakhigh_long = dnakhigh_long %>% mutate(name = fct_relevel(name, (dnakSelect$selectGenes))) 
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dnakhigh_long.Rata',dnakhigh_long)
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dnakSelect.Rata', dnakSelect)
#-----
load(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dnakhigh_long.Rata')
#----
# Run clustering
dnakCluster <- t(dnak[,names(dnakSelect$selectInd)])
dnak.dendro <- as.dendrogram(hclust(d = dist(x = dnakCluster)))
# Create dendro
dendro.plot <- ggdendrogram(data = dnak.dendro, rotate = TRUE)+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_text(size=12, face="bold"))

# Preview the plot
print(dendro.plot)
## putting heatmap and dendrogram together
#grid.newpage()
#print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
### Re-order heatmap rows to match dendrogram
dnak.avg = dnakhigh_long%>%group_by(strain,Pert,condition,slevel,time,name)%>%summarise(mean_value=mean(value))
dnak.order <- order.dendrogram(dnak.dendro)
# Order the levels according to their position in the cluster
dnak.avg$name = factor(x=dnak.avg$name, 
                       levels = colnames(dnak[, names(dnakSelect$selectInd)])[dnak.order],
                       ordered = TRUE)
## recreate the heatmap and reposition the legend
g = ggplot(dnak.avg, aes(slevel, name, fill=mean_value))+
  geom_tile()+
  #scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Red-Green")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        #axis.text.y = element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "top")
grid.newpage()
# print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
## Align dendrogram tips with heatmap rows
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 0.95))
#-------
write.csv(dnakDF, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1LogRatiodnakWholeAdjusted.csv")
## write table of countNames to csv
write.csv(dnakSelect$selectInd, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1DFUFreqdnakWholeAdjusted.csv")


################################################################################################################
logRatioKnockoffdclpb = generateMultiKnockoff(dclpbMtx, 
                                              logRatioKnockoffParamdclpbAdjusted$mu_k, 
                                              logRatioKnockoffParamdclpbAdjusted$Sigma_k,n=100)
dclpbResult = selectVar(dclpbMtx, dclpb,logRatioKnockoffdclpb)
dclpbSelect = selectVarFreq(dclpbResult)
res=factor(dclpb$condition)
dclpbDF = data.frame(res, dclpbMtx[, names(dclpbSelect$selectInd)])
dclpbFit = multinom(res~.,dclpbDF)
summary(dclpbFit) 
dclpbLotus = lotus[which(lotus$locus_tag  %in% dclpbSelect$selectGenes), ]
# write.csv(dclpaLotus, file="~/Documents/TnqData/newdata/dclpaLotus01.csv") FDR=0.2
# write.csv(dclpbLotus, file="~/Documents/TnqData/newdata/dclpbLotus01WholeAdjusted.csv") # FDR=0.1
write.csv(dclpbLotus, file="~/Documents/Research/TnqData/FinalData/UC/dclpbLotus01WholeAdjusted4Rep.csv") 
# write.csv(t(summary(dclpaFit)$coefficients), file="~/Documents/TnqData/newdata/dclpacoefficient.csv") #FDR=0.2
# write.csv(t(summary(dclpbFit)$coefficients), file="~/Documents/TnqData/newdata/dclpbcoefficient01WholeAdjusted.csv") #FDR=0.1
write.csv(t(summary(dclpbFit)$coefficients), file="~/Documents/Research/TnqData/FinalData/UC/dclpbcoefficient01WholeAdjusted4Rep.csv") #FDR=0.1
intersect(dclpbSelect$selectGenes, dclpbSelect$selectGenes)
dfNetDclpb = data.frame(source=dclpbSelect$selectGenes, destination=rep("dclpb", length(dclpbSelect$selectGenes)))

selectDF = cbind(dclpb[, 1:5], dclpb[,names(dclpbSelect$selectInd)])
#wildHigh = selectDF %>% filter(strain=="wild-type" && slevel=="HIGH")
dclpbHigh = selectDF %>% filter(strain=="DCLPB"&condition!="none")
dclpbHigh = dclpbHigh %>%  unite("Pert", condition:slevel, remove = FALSE)
dclpbhigh_long = dclpbHigh%>%pivot_longer(-(1:6))
dclpbhigh_long = dclpbhigh_long %>% mutate(condition = fct_relevel(condition,
                                                                   "canavanine", "heat", "oxidative-peroxide"))
dclpbhigh_long = dclpbhigh_long %>% mutate(slevel = fct_relevel(slevel,
                                                                "LOW", "MEDIUM", "HIGH"))
dclpbhigh_long = dclpbhigh_long %>% mutate(name = fct_relevel(name, (dclpbSelect$selectGenes))) 
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpbhigh_long.Rata',dclpbhigh_long)
save(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpbSelect.Rata', dclpbSelect)
#-----
load(file='/Users/tingtingzhao/Documents/Research/TnqData/FinalData/UC/dclpbhigh_long.Rata')
#------
# Run clustering
dclpbCluster <- t(dclpb[,names(dnakSelect$selectInd)])
dclpb.dendro <- as.dendrogram(hclust(d = dist(x = dclpbCluster)))
# Create dendro
dendro.plot <- ggdendrogram(data = dnak.dendro, rotate = TRUE)+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_text(size=12, face="bold"))

# Preview the plot
print(dendro.plot)
## putting heatmap and dendrogram together
#grid.newpage()
#print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
#print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
### Re-order heatmap rows to match dendrogram
dclpb.avg = dclpbhigh_long%>%group_by(strain,Pert,condition,slevel,time,name)%>%summarise(mean_value=mean(value))
dclpb.order <- order.dendrogram(dclpb.dendro)
# Order the levels according to their position in the cluster
dclpb.avg$name = factor(x=dclpb.avg$name, 
                        levels = colnames(dnak[, names(dclpbSelect$selectInd)])[dclpb.order],
                        ordered = TRUE)
## recreate the heatmap and reposition the legend
g = ggplot(dclpb.avg, aes(slevel, name, fill=mean_value))+
  geom_tile()+
  #scale_y_discrete(limits = unique(rev(wildSelect$selectGenes)),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_continuous_diverging("Blue-Red 3")+
  scale_fill_continuous_diverging(name="Log Ratio", "Red-Green")+
  facet_grid(~condition, 
             labeller = as_labeller(condition_names),
             scales = "free_x", # Let the x axis vary across facets.
             space = "free",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  theme_bw()+
  theme(axis.ticks.y = element_blank(),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face="bold"),
        axis.title = element_blank(),
        #axis.text.y = element_text(size=14, face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        panel.border = element_blank(),
        panel.spacing.x =unit(0,'lines'),
        panel.grid.minor.x = element_blank(),
        panel.grid= element_line(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "top")
grid.newpage()
# print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
# print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))
## Align dendrogram tips with heatmap rows
print(g, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.47, width = 0.2, height = 0.95))
################################################################################################
dclpbGeneMat = dclpbHigh[, 6:dim(dclpbHigh)[2]]
cluster = dclpbHigh$condition
write.csv(dclpbDF, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1LogRatiodclpbWholeAdjusted.csv")
## write table of countNames to csv
write.csv(dclpbSelect$selectInd, file="~/Documents/Research/TnqData/FinalData/UC/dataFDR0.1DFUFreqdclpbWholeAdjusted.csv")

###############################################################################################################
library(igraph)
dfTotal = rbind(rbind(rbind(rbind(dfNetWild, dfNetDlon), dfNetDclpa), dfNetDnak), dfNetDclpb)
## remove the following row from dfTotal
ind1 = intersect(which(dfTotal$source=='CCNA_02037'), which(dfTotal$destination=="dlon"))
ind2 = intersect(which(dfTotal$source=='CCNA_00922'), which(dfTotal$destination=="dclpb"))
ind3 = intersect(which(dfTotal$source=='CCNA_02553'), which(dfTotal$destination=="dclpa"))
ind4 = intersect(which(dfTotal$source=='CCNA_00010'), which(dfTotal$destination=="dnak"))
ind5 = intersect(which(dfTotal$source=='CCNA_00011'), which(dfTotal$destination=="dnaJ"))
# dfTotal = dfTotal[-ind1, ]
## draw a network code
sources <- dfTotal %>%
  distinct(source) %>%
  rename(label = source)

destinations <- dfTotal %>%
  distinct(destination) %>%
  rename(label = destination)

## create nodes
nodes <- full_join(sources, destinations, by = "label")
nodes <- nodes %>% rowid_to_column("id")
nodes

## create edges
per_route <- dfTotal %>%  
  group_by(source, destination) %>%
  summarise(weight = n()) %>% 
  ungroup()
per_route

edges <- per_route %>% 
  left_join(nodes, by = c("source" = "label")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)

edges <- select(edges, from, to, weight)
edges

library(igraph)
routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
#lo <- layout.fruchterman.reingold(routes_igraph, repulserad = vcount(routes_igraph)^2.8, 
#                                  area = vcount(routes_igraph)^2.3, niter = 1000)
#plot(routes_igraph, layout = lo, vertex.size = 3, vertex.frame.color = NULL, 
#     vertex.label.dist = 0.5, vertex.label.cex = 0.7, edge.width = 0.5)
plot(routes_igraph, edge.arrow.size = 0.2)
outdegree = degree(routes_igraph, mode="out")
## mannualy change
outdegree[150:154]=6
library(tidygraph)
library(ggraph)
#load("~/Documents/TnqData/newdata/routes_tidy_old.Rdata")
#load("~/Documents/TnqData/newdata/routes_tidy_df_old.Rdata")
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
#rearrange the rows in the edges tibble to list those with the highest “weight” first, 
#I could use activate() and then arrange()
routes_tidy %>% 
  activate(edges) %>% 
  arrange(desc(weight))

# dh, fr
library(RColorBrewer)
brewer.pal(n = 5, name = "Set1")
# [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00"
colorPlate = palette(brewer.pal(n = 5, name = "Set1"))
## the custom function using Color Brewer
cols_f <- colorRampPalette(RColorBrewer::brewer.pal(6, 'Set1'))

ggraph(routes_tidy, layout = 'fr') + 
  geom_node_point(aes(size = outdegree), alpha=0.7,color=outdegree,
                  #color = cols_f(vcount(routes_tidy)), # custom function for node color
                  show.legend = F) +
  geom_edge_link(aes(width = weight), alpha = 0.5, show.legend = F) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), color=outdegree,size=outdegree*2,repel = TRUE) +
  #labs(edge_width = "Letters") +
  theme_graph()
#save(routes_tidy,file="~/Documents/TnqData/newdata/routes_tidy_adjusted.Rdata")
#save(dfTotal, file="~/Documents/TnqData/newdata/routes_tidy_df_adjusted.Rdata")

save(routes_tidy,file="~/Documents/Research/TnqData/FinalData/UC/routes_tidy_adjusted4Rep.Rdata")
save(dfTotal, file="~/Documents/Research/TnqData/FinalData/UC/routes_tidy_df_adjusted4Rep.Rdata")













