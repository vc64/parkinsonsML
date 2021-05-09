library("annotate")
library("hgu133a.db")
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
pt7621 <- read.table(paste(path,"path/GSE7621_phtype.tsv", sep=""))
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

## change for new gt - naming rows n cols
# gt20164d <- data.frame(gt20164[2:22284,2:12])
# rownames(gt20164d) <- gt20164$V1[2:22284]
# colnames(gt20164d) <- as.character(gt20164[1,2:12])

# mapping probe to gene
# pr20164 <- as.character(rownames(gt20164d))
# come back later to see which gene is good predictor
# sb20164 <- select(hgu133a.db, pr20164, c("SYMBOL", "GENENAME"))

#head(sort(-table(sb20164$PROBEID)))
#gt20164d[rownames(gt20164d)=="AFFX-ThrX-3_at",]
#sum(is.na(sb20164$SYMBOL))

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

# clean up pt first
# pt20164t <- t(pt20164)

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
# head(pt[1:5])
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

#ptdc <- gsub("diagnosis: normal", replace = "ctr", ptds[1,])
#simpleptd <- gsub("diagnosis: Parkinson's disease", replace = "pkd", ptdc)
#simpleptd <- gsub("Superior frontal gyrus control case 2 - A chip", replace = "pkd", ptdc)

simpleptd <- gsub(".*normal.*", replace = "ctr", ptds[1,])
simpleptd <- gsub(".*Parkinson.*", replace = "pkd", simpleptd)
# simpleptd <- gsub(".nuerological.*", replace = "nro", simpleptd)
simpleptd <- gsub(".*control.*", replace = "ctr", simpleptd)
simpleptd <- gsub(".*Control.*", replace = "ctr", simpleptd)

ptds2 <- data.frame(sample=colnames(ptds),group=simpleptd, code = as.character(pt[3,-1]))
ptds3 <- data.frame(ptds2[,2:3])
rownames(ptds3) <- ptds2[,1]

ds <- cbind(ptds3, gtds)

# write.table(ds, file = paste(path,"/normpkdata", quote = FALSE, sep = "/t")

# nulls <- apply(ds, 2, function(r) any(r %in% c("null")))
# which(nulls)
ds <- ds[,1:21907]
colnames(ds)[1] <- "group"

## drop neuro ctr
nctr <- readLines("neuroctr.tsv")
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

# plot solution 
#x <- fit$points[,1]
#y <- fit$points[,2]
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     #main="Metric	MDS",	type="n")
#text(x, y, labels = row.names(fds), cex=.7)


# ggplot mds
rmdsp <- data.frame(dim1=fit$points[,1],dim2=fit$points[,2],
                   group=as.factor(ds2$group), study=ds2$code)

ggplot(data=rmdsp, aes(dim1, dim2, color = study, shape=group)) + 
    geom_point(size = 3) + 
    ggtitle("MDS plot after batch effect correction") + 
    theme(plot.title = element_text(hjust = 0.5))



set.seed(333)
pca <- prcomp(log2(findata[,-1:-2]+1))
pcadata <- data.frame(group=findata$group, study = findata$code, pca$x)
#indata <- createDataPartition(pcadata$group, 1, 0.7, list = FALSE)
#test <- pcadata[-indata,]
#indextv <- pcadata[indata,]
#intrn <- createFolds(indextv$group, 6)
#trainindex <- c(intrn$Fold1, intrn$Fold2, intrn$Fold3, intrn$Fold4, intrn$Fold5)
#intrain <- indata[trainindex]
#ftrn <- pcadata[intrain,]

#invald <- indata[intrn$Fold6]
#fvald <- pcadata[invald,]

set.seed(333)
#svmmodel <- svm(data.matrix(ftrn[,-1]), ftrn$group, type = "C-classification")

#ctrl <- trainControl(method = "repeatedcv", repeats = 100)
#mod <- train(group ~., data=ftrn, method = "svmLinear", trControl = ctrl)

#svm_model <- svm(as.numeric(as.factor(group)) ~ ., data = ftrn)

#sv <- svm(as.factor(group) ~ ., data = ftrn, type = "C-classification")

#tuning <- tune(svm, as.numeric(as.factor(group)) ~., ftrn,
               #kernel="radial", type = "C-classification",
               #ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

#tuning <- tune(svm, as.numeric(as.factor(group)) ~., ftrn,
               #ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
               #tunecontrol = tune.control(sampling = "fix"))
               
# tuning svm
#tobj <- tune.svm(as.numeric(as.factor(group))~., data = ftrn, gamma = 2^(-1:1), cost = 2^(2:4))
#tobj <- tune.svm(as.factor(group) ~., data = ftrn, kernel = "linear"
                #, gamma = 2^(-1:1), cost = 2^(2:4))
#tobj <- tune.svm(as.factor(group) ~., data = ftrn, kernel = "sigmoid"
                 #, gamma = 2^(-1:1), cost = 2^(2:4))
# training using svm func
#svm <- svm(as.factor(group) ~ ., data = ftrn, kernel = "linear", type = "C-classification", 
           #gamma = 0.5, cost = 4)
#svm <- svm(as.factor(group) ~ ., data = ftrn, kernel = "sigmoid", type = "C-classification", 
           #gamma = 2, cost = 4)

# training using train(method = svm)
set.seed(333)
train_control <- trainControl(method="repeatedcv", number = 5, repeats = 10)

#svm <- train(as.factor(group) ~ ., data = pcadata, kernel = "sigmoid",
             #type = "C-classification", 
             #gamma = 2, cost = 4, trControl = train_control,
             #method = "svmLinear2")

svm <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
                          method = "svmLinearWeights")
svm <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
             method = "svmLinear2", tuneLength = 10)
svmrad <- train(as.factor(group) ~ ., data = pcadata, trControl = train_control,
             method = "svmRadial", tuneLength = 10)

#pred <- predict(svm, fvald)
#table(pred, fvald$group)

#fit <- vglm(group~., family=binomial, data=ftrn)

#glmmod <- glmnet(data.matrix(ftrn[,-1]), 
          #y= as.factor(ftrn$group), alpha=1, family="binomial")
#lasso.pred <- predict(glmmod, s = bestlam, newx = x[test,])


lambda <- 10^seq(10, -2, length = 100)
cv.out <- cv.glmnet(data.matrix(ftrn[,-1]), as.numeric(as.factor(ftrn$group)),
                    alpha = 0, family = "binomial", type.measure = "mse")
minlam <- cv.out$lambda.min
bestlam <- cv.out$lambda.1se

x_vald <- model.matrix(as.numeric(as.factor(group)) ~., fvald)
lassoprob <- predict(cv.out, newx = x_vald, s=bestlam, type = "response")

lasso.mod <- glmnet(data.matrix(ftrn[,-1]), as.numeric(as.factor(ftrn$group)),
                     alpha = 1, lambda = lambda)
lassomod <- train(as.factor(group) ~ ., data = pcadata, 
                  trControl = train_control, method = "lasso")

lasso.pred <- predict(lasso.mod, s = bestlam, newx = data.matrix(fvald[,-1]))
mean((lasso.pred-ytest)^2)

set.seed(333)
tr5 <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
naivemod <- train(as.factor(group) ~ ., data = pcadata, 
                  trControl = tr5, method = "nb")


## nb
train_control <- trainControl(method="repeatedcv", number=5, repeats=5)
# train the model
mnb <- train(as.factor(group)~., data=pcadata, trControl=train_control, method="nb")
# summarize results
print(mnb)
confusionMatrix(mnb)


## rf
# no tuning
mrf1 <- train(as.factor(group) ~ ., data = pcadata, 
                trControl = train_control, method = "rf")
print(mrf1)
confusionMatrix(mrf1)

# tuning: grid search
seed=333
control <- trainControl(method="repeatedcv", number=5, repeats=10, search="grid")
set.seed(seed)
tunegrid1 <- expand.grid(.mtry=c(1:20))
mrf2 <- train(as.factor(group)~., data=pcadata, method="rf",
                       metric='Accuracy', tuneGrid=tunegrid1, trControl=control)
print(mrf2)
plot(mrf2)
confusionMatrix(mrf2)


## custom rf
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

control3 <- trainControl(method="repeatedcv", number=5, repeats=3)
mnet <- train(group ~ ., data = pcadata, trControl = control3, method = "dnn")

set.seed(333)
#control3 <- trainControl(method="repeatedcv", number=5, repeats=6)
mnet <- train(group ~ ., data = pcadata, trControl = control2, 
              method = "nnet", 
              tuningGrid = data.frame(size=10, decay = -0.1))

#mnetx <- train(group ~ ., data = pcadata, trControl = control2, 
               #method = "mxnet")
#generalized linear model

mglm <- train(group ~ ., data = pcadata, trControl = control2, method = "glm")
mglmnet <- train(group ~., data = pcadata, trControl = control2, method = "glmnet")
#mglmnet <- train(group ~., data = pcadata, trControl = control2, method = "glmnet", tuningLength = 5)

# knn
mknn <- train(group ~ ., data = pcadata, trControl = control2, method = "knn")
#mknn <- train(group ~ ., data = pcadata, trControl = control2, method = "knn", tuningLength = 10)
#mknn <- train(group ~ ., data = pcadata, trControl = control2, method = "knn", 
              #tuningGrid = data.frame(k = 100))

#mkknn <- train(as.factor(group) ~., data = pcadata, trControl = control2, method = "kknn")

#ldamod <- train(as.factor(group) ~ ., data = pcadata, )

# boostedtrees
mtbst <- train(as.factor(group) ~., data = pcadata, 
               trControl = control2, method = "adaboost")


# plot mnet
require(RCurl)

root.url<-'https://gist.github.com/fawda123'
raw.fun<-paste(
    root.url,
    '5086859/raw/17fd6d2adec4dbcf5ce750cbd1f3e0f4be9d8b19/nnet_plot_fun.r',
    sep='/'
)
script<-getURL(raw.fun, ssl.verifypeer = FALSE)
eval(parse(text = script))
rm('script','raw.fun')



# roc
roc_mnet <- roc(pcadata$group, mnet$pred)



