library(caret)
library(pROC)
library(e1071)
library(preprocessCore)
library(kknn)
library(fastAdaboost)
library(limma)

options(stringsAsFactors = F)
set.seed(333)

#path <- "insert path here"
path <- "C:/Users/viche/Desktop/vc/code/foru"

pt6613 <- read.table(paste(path, "/GSE6613_phtype.tsv", sep=""))
gt6613 <- read.table(paste(path, "/GSE6613_gtype.tsv", sep=""))
gt7621 <- read.table(paste(path,"/GSE7621_gtype.tsv", sep=""))
pt7621 <- read.table(paste(path,"/GSE7621_phtype.tsv", sep=""))
pt8397_96 <- read.table(paste(path,"/GSE8397-96_phtype.tsv", sep=""))
pt8397_97 <- read.table(paste(path,"/GSE8397-97_phtype.tsv", sep=""))
gt8397_97 <- read.table(paste(path,"/GSE8397-97_gtype.tsv", sep=""))
gt8397_96 <- read.table(paste(path,"/GSE8397-96_gtype.tsv", sep=""))
gt20141 <- read.table(paste(path,"/GSE20141_gtype.tsv", sep=""))
pt20141 <- read.table(paste(path,"/GSE20141_phtype.tsv", sep=""))
pt20163 <- read.table(paste(path,"/GSE20163_phtype.tsv", sep=""))
gt20163 <- read.table(paste(path,"/GSE20163_gtype.tsv", sep=""))
gt20164 <- read.table(paste(path,"/GSE20164_gtype.tsv", sep=""))
pt20164 <- read.table(paste(path,"/GSE20164_phtype.tsv", sep=""))
pt20291 <- read.table(paste(path,"/GSE20291_phtype.tsv", sep=""))
gt20291 <- read.table(paste(path,"/GSE20291_gtype.tsv", sep=""))
gt20292 <- read.table(paste(path,"/GSE20292_gtype.tsv", sep=""))
pt20292 <- read.table(paste(path,"/GSE20292_phtype.tsv", sep=""))
#pt20333 <- read.table(paste(path,"/GSE20333_phtype.tsv", sep=""))
#gt20333 <- read.table(paste(path,"/GSE20333_gtype.tsv", sep=""))
#gt24378 <- read.table(paste(path,"/GSE24378_gtype.tsv", sep=""))
#pt24378 <- read.table(paste(path,"/GSE24378_phtype.tsv", sep=""))

# other nuerodegenerative stuff
pt26927 <- read.table(paste(path,"/GSE26927_phtype.tsv", sep=""))
gt26927 <- read.table(paste(path,"/GSE26927_gtype.tsv", sep=""))

# merging data with hgu133A.db
gt <- merge(gt20164, gt8397_96, by='V1')
gt <- merge(gt, gt20292, by='V1')
gt <- merge(gt, gt20163, by='V1')
gt <- merge(gt, gt6613, by='V1')
gt <- merge(gt, gt20291, by='V1')
# originally 56,000 something probes, but shared 22278/22284 probes with other hgu133a.db datasets
gt <- merge(gt, gt20141, by = "V1")
# shares 21942/22278 probes, will temporarily add for now
gt <- merge(gt, gt7621, by = "V1")



#add sample accession
pt20164[nrow(pt20164)+1,] <- "GSE20164"
pt20164[nrow(pt20164), 1] <- "code"
pt20141[nrow(pt20141)+1,] <- "GSE20141"
pt20141[nrow(pt20141), 1] <- "code"
pt20163[nrow(pt20163)+1,] <- "GSE20163"
pt20163[nrow(pt20163), 1] <- "code"
pt7621[nrow(pt7621)+1,] <- "GSE7621"
pt7621[nrow(pt7621), 1] <- "code"
pt20291[nrow(pt20291)+1,] <- "GSE20291"
pt20291[nrow(pt20291), 1] <- "code"
pt20292[nrow(pt20292)+1,] <- "GSE20292"
pt20292[nrow(pt20292), 1] <- "code"
pt8397_96[nrow(pt8397_96)+1,] <- "GSE8397_96"
pt8397_96[nrow(pt8397_96), 1] <- "code"
pt6613[nrow(pt6613)+1,] <- "GSE6613"
pt6613[nrow(pt6613), 1] <- "code"

#change rows to merge by
pt8397_96[1,1] <- "!Sample_characteristics_ch1"

#6613 does not have any gender or age :(
pt <- merge(pt20164[c(1,4, nrow(pt20164)),], pt8397_96[c(1, 2, nrow(pt8397_96)),], by='V1')
pt <- merge(pt, pt20292[c(1, 2, nrow(pt20292)),], by='V1')
pt <- merge(pt, pt20163[c(1, 4, nrow(pt20163)),], by='V1')
pt <- merge(pt, pt6613[c(1, 2, nrow(pt6613)),], by='V1')
pt <- merge(pt, pt20291[c(1, 2, nrow(pt20291)),], by='V1')
pt <- merge(pt, pt20141[c(1, 2, nrow(pt20141)),], by = "V1")
pt <- merge(pt, pt7621[c(1, 2, nrow(pt7621)),], by='V1')

# transposing the genotype data so that GSE is rownames
tgt <- t(gt)
gtds <- data.frame(tgt[2:288, 1:21941])
# setting rownames and colnames for gt and pt
rownames(gtds) <- tgt[2:288, 21942]
colnames(gtds) <- as.character(tgt[1, 1:21941])

ptds <- data.frame(pt[c(1, 3), 2:288])
rownames(ptds) <- pt$V1[c(1, 3)]
colnames(ptds) <- as.character(pt[2, 2:288])


simpleptd <- gsub(".*normal.*", replace = "ctr", ptds[1,])
simpleptd <- gsub(".*Parkinson.*", replace = "pkd", simpleptd)
simpleptd <- gsub(".*control.*", replace = "ctr", simpleptd)
simpleptd <- gsub(".*Control.*", replace = "ctr", simpleptd)

ptds2 <- data.frame(sample=colnames(ptds),group=simpleptd, code = as.character(pt[3,-1]))
ptds3 <- data.frame(ptds2[,2:3])
rownames(ptds3) <- ptds2[,1]

ds <- cbind(ptds3, gtds)


ds <- ds[,1:21907]
colnames(ds)[1] <- "group"

## drop neuro ctr
nctr <- readLines(paste(path, "/neuroctr.tsv", sep=""))
ds2 <- ds[!rownames(ds) %in% nctr,]

numds <- as.data.frame(apply(ds2[,-(1:2)], 2, as.numeric))
rownames(numds) <- rownames(ds2)
numdst <- data.frame(t(numds))



logdata = log2(numdst + 1)
logdata = logdata[rowMeans(logdata) > 3, ]
colramp = colorRampPalette(c(3,"white",2))(dim(numdst)[2])
plot(density(logdata[,1]),col=colramp[1],lwd=3,ylim=c(0,.80), 
     main = "Distribution before normalization")
for(i in 2:dim(numdst)[2]){lines(density(logdata[,i]),lwd=3,col=colramp[i])}


norm_ds = normalize.quantiles(as.matrix(logdata))
plot(density(norm_ds[,1]),col=colramp[1],lwd=3,ylim=c(0,.30), 
     main = "Distribution after normalization")
for(i in 2:dim(numdst)[2]){lines(density(norm_ds[,i]),lwd=3,col=colramp[i])}

rmbds <- removeBatchEffect(norm_ds, ds2$code,
                           design=model.matrix(~0+group, data=ds2))

# row & col names
fds <- data.frame(t(data.frame(rmbds)))
rownames(fds) <- rownames(ds2)
grouplabels <- ds2[,1:2]
findata <- cbind(grouplabels, fds)


# mds plot
d <- dist(fds) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim



# ggplot mds
rmdsp <- data.frame(dim1=fit$points[,1],dim2=fit$points[,2],
                    group=as.factor(ds2$group), study=ds2$code)

ggplot(data=rmdsp, aes(dim1, dim2, color = study, shape=group)) + 
    geom_point(size = 3) + 
    ggtitle("MDS plot after batch effect correction") + 
    theme(plot.title = element_text(hjust = 0.5))



set.seed(333)
pca <- prcomp(log2(findata[,-1:-2]+2))
#pca <- prcomp(log2(findata[,2:20725]+1))
pcadata <- data.frame(group=findata$group, study = findata$code, pca$x)


set.seed(333)


# training using train(method = svm)
set.seed(333)
train_control <- trainControl(method="repeatedcv", number = 5, repeats = 10)


svm <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
             method = "svmLinearWeights")
svm <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
             method = "svmLinear2", tuneLength = 10)
svmrad <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
                method = "svmRadial", tuneLength = 10)

# too many dependencies, not important
#set.seed(333)
#tr5 <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
#naivemod <- train(as.factor(group) ~ ., data = pcadata, 
#                  trControl = tr5, method = "nb")


## nb
#train_control <- trainControl(method="repeatedcv", number=5, repeats=5)
# train the model
#mnb <- train(as.factor(group)~., data=pcadata, trControl=train_control, method="nb")
# summarize results
#print(mnb)
#confusionMatrix(mnb)


## custom rf (VERY slow)
Sys.time()
library("randomForest")
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

set.seed(seed)
control2 <- trainControl(method="repeatedcv", number=5, repeats=10)
tunegrid2 <- expand.grid(.mtry=c(1:20), .ntree=c(1000, 1500, 2000, 2500))
mrfc <- train(as.factor(group)~., data=pcadata, method=customRF,
              metric='Accuracy', tuneGrid=tunegrid2, trControl=control2)
Sys.time()
summary(mrfc)
plot(mrfc)


treemd1 <- train(group ~ ., method = "rpart", data = pcadata,
                 trControl = control2)


# neural network
set.seed(333)
mnet <- train(group ~ ., data = pcadata, trControl = control2, 
              method = "nnet", tuningLength = 1, MaxNWts = 1317)



# boostedtrees
set.seed(333)
mtbst <- train(as.factor(group) ~., data = pcadata, 
               trControl = control2, method = "adaboost")


# model         accuracy

# svm           0.569
# svmRadial     0.573
# naivebayes    dependency issues
# randomForest  takes too long
# nnet          0.793
# boostedTrees  0.928??? overfitting?


# adaboost output

# > mtbst
# AdaBoost Classification Trees 
# 
# 254 samples
# 255 predictors
# 2 classes: 'ctr', 'pkd' 
# 
# No pre-processing
# Resampling: Cross-Validated (5 fold, repeated 10 times) 
# Summary of sample sizes: 203, 203, 203, 204, 203, 203, ... 
# Resampling results across tuning parameters:
#     
#     nIter  method         Accuracy   Kappa    
# 50    Adaboost.M1    0.9043686  0.8032668
# 50    Real adaboost  0.8772157  0.7457349
# 100    Adaboost.M1    0.9236784  0.8437119
# 100    Real adaboost  0.8850745  0.7614762
# 150    Adaboost.M1    0.9280314  0.8525355
# 150    Real adaboost  0.8933490  0.7787139
# 
# Accuracy was used to select the optimal model using the largest value.
# The final values used for the model were nIter = 150 and method = Adaboost.M1.