
# List of packages for session
.packages = c("tidyverse", "dplyr", "caret", "reshape2", "gmodels", "expss", "gridExtra",
              "pROC", "ggthemes", "ROCR", "earth", "caretEnsemble", "caTools", "e1071",
              "randomForest")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

# CM table intro
cm.table <- matrix(c("True Positive (TP)", "False Positive (FP)", "False Negative (FN)", 
                     "True Negative (TN)"), ncol=2, byrow = TRUE)
rownames(cm.table) <- c("Predicted: Yes", "Predicted: No")
colnames(cm.table) <- c("Actual: Yes", "Actual: No")
cm.table <- as.table(cm.table)
cm.table

#download data
url <- "https://raw.githubusercontent.com/PiperS52/Liver-Disease-ML/master/indian_liver_patient.csv"
rawdata <- read.csv(url)
rawdata[0:20,]

str(rawdata)
head(rawdata)
rawdata <- rawdata %>% rename(TB = Total_Bilirubin, DB = Direct_Bilirubin, ALP = Alkaline_Phosphotase,
                   ALT = Alamine_Aminotransferase, AST = Aspartate_Aminotransferase,
                   TP = Total_Protiens, A = Albumin, AG_Ratio = Albumin_and_Globulin_Ratio,
                   Disease = Dataset)
# original data 1 = Disease, 2 = No disease
rawdata$Disease <- replace(rawdata$Disease, rawdata$Disease==2, "No")
rawdata$Disease <- replace(rawdata$Disease, rawdata$Disease==1, "Yes")
# create and modify gender col
rawdata$Male <- as.character(rawdata$Gender)
rawdata$Male <- replace(rawdata$Male, rawdata$Male == "Male", 1)
rawdata$Male <- replace(rawdata$Male, rawdata$Male =="Female", 0)
rawdata$Male <- as.numeric(rawdata$Male)

# addressing missing values
colSums(is.na(rawdata))
which(is.na(rawdata$AG_Ratio))
rawdata$AG_Ratio <- replace(rawdata$AG_Ratio, rawdata$AG_Ratio=="", NA)
rawdata$AG_Ratio <- as.numeric(rawdata$AG_Ratio)
rawdata$AG_Ratio
# 4 missing values in dataset
rawdata[210,]
rawdata[242,]
rawdata[254,]
rawdata[313,]
mean(rawdata$AG_Ratio, na.rm = TRUE)
# replacing the missing values with the mean of that column
rawdata$AG_Ratio[is.na(rawdata$AG_Ratio)] <- 
  mean(rawdata$AG_Ratio, na.rm = TRUE)

# checking multicollinearity - drop some variables
cor_dat <- rawdata %>% select(Age, TB, DB, ALP, ALT, AST, TP, A, AG_Ratio)
cormat <- round(cor(cor_dat),2)
cormat

melted_cormat <- melt(cormat)
melted_cormat

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
# Melt the correlation matrix

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# concerned DB and TB are highly corr, therefore create new col for DB/TB ratio 
# strong corr between A and AG_Ratio, therefore create Globulin = A/AG_Ratio
rawdata <- rawdata %>% mutate(DB_TB = DB / TB, G = A / AG_Ratio)
# check corr for A, TP and G
cor_dat <- rawdata %>% select(A, TP, G)
cormat <- round(cor(cor_dat),2)
cormat

# as A is strongly corr with TP but not G, and TP is largely composed of Albumin and Globulin TP can be dropped
rawdata <- rawdata %>% select(-TP,-TB,-DB,-AG_Ratio)
# drop DB, TB, TP and AG_Ratio)
# AST and ALT are strongly correlated - although both left in as separate predictors as concerned
# that if omitted or altered it would weaken the predictive power of the model + inflate 
# p values from inc residual variance

###### Describing the data ######
# density plot for each num variable - overlapping disease/no disease
plot1 <- rawdata %>% ggplot(aes(Age, fill = Disease, group = Disease)) + 
  geom_density(aes(x = Age), alpha = 0.2, stat = "density")
plot1
plot2 <- rawdata %>% ggplot(aes(DB_TB, fill = Disease, group = Disease)) + 
  geom_density(aes(x = DB_TB), alpha = 0.2, stat = "density")
plot2
plot3 <- rawdata %>% ggplot(aes(ALP, fill = Disease, group = Disease)) + 
  geom_density(aes(x = ALP), alpha = 0.2, stat = "density")
plot3
plot4 <- rawdata %>% ggplot(aes(ALT, fill = Disease, group = Disease)) + 
  geom_density(aes(x = ALT), alpha = 0.2, stat = "density")
plot4
plot5 <- rawdata %>% ggplot(aes(AST, fill = Disease, group = Disease)) + 
  geom_density(aes(x = AST), alpha = 0.2, stat = "density")
plot5
plot6 <- rawdata %>% ggplot(aes(A, fill = Disease, group = Disease)) + 
  geom_density(aes(x = A), alpha = 0.2, stat = "density")
plot6 # low levels of A typically indicative of disease
plot7 <- rawdata %>% ggplot(aes(G, fill = Disease, group = Disease)) + 
  geom_density(aes(x = G), alpha = 0.2, stat = "density")
plot7 # low levels of G typically indicative of disease
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, nrow = 3)

# boxplot for each num variable - disease/no disease
plot1 <- rawdata %>% ggplot(aes(x = Disease, y = Age, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot1
plot2 <- rawdata %>% ggplot(aes(x = Disease, y = DB_TB, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot2
plot3 <- rawdata %>% ggplot(aes(x = Disease, y = ALP, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot3
plot4 <- rawdata %>% ggplot(aes(x = Disease, y = ALT, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot4
plot5 <- rawdata %>% ggplot(aes(x = Disease, y = AST, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot5
plot6 <- rawdata %>% ggplot(aes(x = Disease, y = A, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot6
plot7 <- rawdata %>% ggplot(aes(x = Disease, y = G, fill = Disease, group = Disease)) + 
  geom_boxplot(aes(x = Disease), stat = "boxplot", position = "dodge2", show.legend = FALSE) 
plot7
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, nrow = 3)

# barplot of disease/no disease - col = gender as well as table
plot_gen<- rawdata %>% ggplot(aes(x = Disease, fill = Gender, group = Gender)) + 
  geom_bar(aes(x = Disease), stat = "count", position = "stack", show.legend = TRUE) +
  geom_text(stat="count", aes(label=..count..), position=position_stack(vjust=0.2), 
            colour = "white",vjust=-1)
plot_gen
# table of disease by gender and 
table <- table(rawdata$Gender, rawdata$Disease, dnn=list("Gender", "Disease"))

CrossTable(table, expected = FALSE, prop.r = TRUE, prop.c = TRUE, prop.t = FALSE, 
           chisq = FALSE, prop.chisq = FALSE, fisher = FALSE, mcnemar = FALSE,
           resid = FALSE, sresid = FALSE, aresid = FALSE)
##### Method & Analysis ######

str(rawdata)
#creating the train/test sets
set.seed(755, sample.kind = "Rounding")
test_index <- createDataPartition(y = rawdata$Disease, times = 1, p = 0.2, list = FALSE)
train_set <- rawdata[-test_index,]
test_set <- rawdata[test_index,]

# upsampling (so disease vs no disease are equal) - data_tr currently has 350 'Disease', 150 'No Disease'
'%ni%' <- Negate('%in%')
options(scipen = 999)
set.seed(100, sample.kind ="Rounding")
up_train <- upSample(x = train_set[, colnames(train_set) %ni% "Class"], y = factor(train_set$Disease))
dim(up_train)

#plot to show balanced data following upsampling
plot_gen <- up_train %>% ggplot(aes(x = Disease, fill = Gender, group = Gender)) + 
  geom_bar(aes(x = Disease), stat = "count", position = "stack", show.legend = TRUE) +
  geom_text(stat="count", aes(label=..count..), position=position_stack(vjust=0.2), 
            colour = "white",vjust=-2)
plot_gen

# data prep
up_train <- up_train %>% select(-Gender)
up_train$Disease <- as.factor(up_train$Disease)
test_set$Disease <- as.factor(test_set$Disease)

str(up_train)

#### Model 1 - Logistic regresssion 
# 1.1 Simple Logistic
set.seed(1, sample.kind = "Rounding")
fit_glm <- glm(Disease ~ Age + A + G + DB_TB + log(ALP) + log(ALT) + log(AST) + Male, 
               family = "binomial", data = up_train)
summary(fit_glm) # coef interp

# 1.2 k-fold (trControl = cv, number = 5)
ctrl_glm <- trainControl(method = "cv", number = 5) # k-fold with 5 
up_train$Disease <- as.factor(up_train$Disease)
set.seed(1, sample.kind = "Rounding")
fit_glmk <- train(Disease ~ Age + A + G + DB_TB + Male + log(ALP) + log(ALT) + log(AST), 
                  method = "glm", family = "binomial", data = up_train,
                  trControl = ctrl_glm)


# graphs
pred <- predict(fit_glmk, newdata = test_set, type = "prob") 
pred <- prediction(pred[,2], factor(test_set$Disease, levels = c("No", "Yes"), ordered = TRUE))
pred

# ROC/AUC plot
perf1 <- performance(pred, "tpr", "fpr") # ROC plot
ruc.plot <- plot(perf1, colorize = TRUE, main="ROC and AUC", ylab="Sensitivity", xlab="1 - Specificity",
     lwd=4)
abline(a=0, b=1, lty = 4)
# calc AUC
auc1 <- performance(pred, "auc")
auc1 <- unlist(slot(auc1, "y.values"))
auc1 <- round(auc1, 4)
legend(.6, .3, auc1, title = "AUC", cex = 1.25, bty = "n")

# MCC
perf2 <- performance(pred, "mat") # MCC value
rec.plot <- plot(perf2, col = "black", main="Matthew's Correlation Coefficient", ylab="MCC", xlab="Threshold",
                 lwd=2)
max <- which.max(slot(perf2, "y.values")[[1]])
mcc1 <- slot(perf2, "y.values")[[1]][max]
mcc1 <- round(mcc1, 4)
legend(.5, .2, mcc1, title = "MCC", cex = 1.25, bty = "n")

# CM
y_hat_glmk <- predict(fit_glmk, newdata = test_set) # works for CM (type = response by def)
cm_glmk <-confusionMatrix(data = factor(y_hat_glmk), reference = factor(test_set$Disease), 
                positive = "Yes", mode="everything")
cm_glmk
cm_glmk_tab <- as.matrix(cm_glmk, what = "classes")

####   Model 2 - kNN
ctrl_k <- trainControl(method = "cv", number = 5) # search = "random" to get idea of range for k
set.seed(1, sample.kind = "Rounding")
fit_knn <- train(Disease ~ Age + A + G + DB_TB + Male + log(ALP) + log(ALT) + log(AST), 
                  method = "knn", data = up_train, tuneGrid = data.frame(k = seq(3,15,2)),
                 trControl = ctrl_k) 
plot(fit_knn)
fit_knn$bestTune
varImp(fit_knn)

# graphs
pred <- predict(fit_knn, newdata = test_set, type = "prob") 
pred <- prediction(pred[,2], factor(test_set$Disease, levels = c("No", "Yes"), ordered = TRUE))
pred

# ROC/AUC plot
perf1 <- performance(pred, "tpr", "fpr") # ROC plot
ruc.plot <- plot(perf1, colorize = TRUE, main="ROC and AUC", ylab="Sensitivity", xlab="1 - Specificity",
                 lwd=4)
abline(a=0, b=1, lty = 4)
# calc AUC
auc2 <- performance(pred, "auc")
auc2 <- unlist(slot(auc2, "y.values"))
auc2 <- round(auc2, 4)
legend(.6, .3, auc2, title = "AUC", cex = 1.25, bty = "n")

# MCC
perf2 <- performance(pred, "mat") # MCC value
rec.plot <- plot(perf2, col = "black", main="Matthew's Correlation Coefficient", ylab="MCC", xlab="Threshold",
                 lwd=2)
max <- which.max(slot(perf2, "y.values")[[1]])
mcc2 <- slot(perf2, "y.values")[[1]][max]
mcc2 <- round(mcc2, 4)
legend(.6, .34, mcc2, title = "MCC", cex = 1.25, bty = "n")

# CM
y_hat_knn <- predict(fit_knn, newdata = test_set)
cm_knn <-confusionMatrix(data = factor(y_hat_knn), reference = factor(test_set$Disease), 
                          positive = "Yes", mode="everything")
cm_knn
cm_knn_tab <- as.matrix(cm_knn, what = "classes")
cm_knn_tab

#### Model 3 - Random Forest
ctrl_rf <- trainControl(method = "cv", number = 5) # k-fold with 5
set.seed(1, sample.kind = "Rounding")
fit_rf <- train(Disease ~ Age + A + G + DB_TB + Male + log(ALP) + log(ALT) + log(AST), 
                  method = "rf", data = up_train, trControl = ctrl_rf,
                tuneGrid = data.frame(mtry = seq(1,7,1))) 
fit_rf$bestTune
plot(fit_rf)
varImp(fit_rf)

# graphs
pred <- predict(fit_rf, newdata = test_set, type = "prob") 
pred <- prediction(pred[,2], factor(test_set$Disease, levels = c("No", "Yes"), ordered = TRUE))
pred
# ROC/AUC plot
perf1 <- performance(pred, "tpr", "fpr") # ROC plot
ruc.plot <- plot(perf1, colorize = TRUE, main="ROC and AUC", ylab="Sensitivity", 
                 xlab="1 - Specificity", lwd=4)
abline(a=0, b=1, lty = 4)
# calc AUC
auc3 <- performance(pred, "auc")
auc3 <- unlist(slot(auc3, "y.values"))
auc3 <- round(auc3, 4)
legend(.6, .3, auc3, title = "AUC", cex = 1.25, bty = "n")

# MCC
perf2 <- performance(pred, "mat") # MCC value
rec.plot <- plot(perf2, col = "black", main="Matthew's Correlation Coefficient", ylab="MCC", xlab="Threshold",
                 lwd=2)
max <- which.max(slot(perf2, "y.values")[[1]])
mcc3 <- slot(perf2, "y.values")[[1]][max]
mcc3 <- round(mcc3, 4)
legend(.5, .2, mcc3, title = "MCC", cex = 1.25, bty = "n")

# CM
y_hat_rf <- predict(fit_rf, newdata = test_set)
cm_rf <- confusionMatrix(data = factor(y_hat_rf), reference = factor(test_set$Disease), 
                           positive = "Yes", mode="everything")
cm_rf
cm_rf_tab <- as.matrix(cm_rf, what = "classes")
cm_rf_tab

#### Model 4 - Earth (Multivariate Adaptive Regression Spline)
ctrl_e <- trainControl(method = "cv", number = 5) # k-fold with 5
set.seed(1, sample.kind = "Rounding")
fit_e <- train(Disease ~ Age + A + G + DB_TB + Male + log(ALP) + log(ALT) + log(AST), 
                method = "earth", data = up_train, trControl = ctrl_e, 
               tuneGrid = expand.grid(.degree=1, .nprune=(1:10)*2))

fit_e$bestTune
plot(fit_e)
varImp(fit_e)

# graphs
pred <- predict(fit_e, newdata = test_set, type = "prob") 
pred <- prediction(pred[,2], factor(test_set$Disease, levels = c("No", "Yes"), ordered = TRUE))

# ROC/AUC plot
perf1 <- performance(pred, "tpr", "fpr") # ROC plot
ruc.plot <- plot(perf1, colorize = TRUE, main="ROC and AUC", ylab="Sensitivity", xlab="1 - Specificity", lwd=4)
abline(a=0, b=1, lty = 4)
# calc AUC
auc4 <- performance(pred, "auc")
auc4 <- unlist(slot(auc4, "y.values"))
auc4 <- round(auc4, 4)
legend(.6, .3, auc4, title = "AUC", cex = 1.25, bty = "n") # AUC = 0.8022

# MCC
perf2 <- performance(pred, "mat") # MCC value
rec.plot <- plot(perf2, col = "black", main="Matthew's Correlation Coefficient", ylab="MCC", xlab="Threshold",
                 lwd=2)
max <- which.max(slot(perf2, "y.values")[[1]])
mcc4 <- slot(perf2, "y.values")[[1]][max]
mcc4 <- round(mcc4, 4)
legend(.5, .2, mcc4, title = "MCC", cex = 1.25, bty = "n") # mcc = 0.5155

# CM
y_hat_e <- predict(fit_e, newdata = test_set)
cm_e <- confusionMatrix(data = factor(y_hat_e), reference = factor(test_set$Disease), 
                         positive = "Yes", mode="everything")

cm_e
cm_e_tab <- as.matrix(cm_e, what = "classes")
cm_e_tab


#### Model 5 - Ensemble
ctrl_en <- trainControl(method = "cv", number = 5, savePredictions = "final", classProbs = TRUE)
models <- c("knn", "earth")
set.seed(1, sample.kind = "Rounding")
fit_en <- caretList(Disease ~ Age + A + G + DB_TB + Male + log(ALP) + log(ALT) + log(AST), 
                    data = up_train, methodList = models,
                    trControl = ctrl_en)

stackControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
greedy_ensemble <- caretEnsemble(fit_en, trControl = stackControl, metric="ROC")
greedy_ensemble$models

varImp(greedy_ensemble)

modelCor(resamples(fit_en)) # shows only small +ve corr between knn and earth
# graphs
pred <- lapply(fit_en, predict, newdata = test_set, type = "prob")
pred <- lapply(pred, function(x) x[,"Yes"])
pred <- data.frame(pred)
ens_preds <- predict(greedy_ensemble, newdata = test_set, type = "prob")
pred$ensemble <- ens_preds
pred$ensemble # this object is a list of probs 118 entries

caTools::colAUC(pred, test_set$Disease) # reporting AUC
# graphs
p <- as.data.frame(pred$ensemble)
colnames(p)[colnames(p) == 'pred$ensemble'] <- 'Yes'
p <- p %>% mutate(No = 1-Yes)
p
pred <- prediction(p[,2], factor(test_set$Disease, levels = c("No", "Yes"), ordered = TRUE))
pred
# ROC/AUC plot
perf1 <- performance(pred, "tpr", "fpr") # ROC plot
ruc.plot <- plot(perf1, colorize = TRUE, main="ROC and AUC", ylab="Sensitivity", xlab="1 - Specificity",
                 lwd=4)
abline(a=0, b=1, lty = 4)
# calc AUC
auc5 <- performance(pred, "auc")
auc5 <- unlist(slot(auc5, "y.values"))
auc5 <- round(auc5, 4)
legend(.6, .3, auc5, title = "AUC", cex = 1.25, bty = "n")
# MCC
perf2 <- performance(pred, "mat") # MCC value
mcc.plot <- plot(perf2, col = "black", xlab="Threshold",
                 lwd=2)
max <- which.max(slot(perf2, "y.values")[[1]])
mcc5 <- slot(perf2, "y.values")[[1]][max]
mcc5 <- round(mcc5, 4)
legend(.5, .2, mcc5, title = "MCC", cex = 1.25, bty = "n")

# Precision
perf3 <- performance(pred, "prec") # accuracy for diff threshold values
prec.plot <- plot(perf3, col = "blue", xlab="Threshold",
                 lwd=2)
# Recall
perf4 <- performance(pred, "rec") # accuracy for diff threshold values
rec.plot <- plot(perf4, col = "red", xlab="Threshold",
                 lwd=2)
# F1
perf5 <- performance(pred, "f") # accuracy for diff threshold values
f1.plot <- plot(perf5, col = "green", xlab="Threshold",
                 lwd=2)
# final plot
plot(perf2, ylab="", col=1, lty=2, xlim=c(0,1), ylim=c(0,1), lwd=3)
plot(perf3, col=2, lty=3, lwd=3, add=TRUE)
plot(perf4, col=3, lty=4, lwd=3, add=TRUE)
plot(perf5, col=4, lty=5, lwd=3, add=TRUE)
legend(.8,.93,legend=c("MCC","Precision","Recall","F1"), col=c("black","blue","red","green"), lty = 2:3, cex=1.2)

# CM
y_hat_ensemble
y_hat_ensemble <- predict(greedy_ensemble, newdata = test_set)
cm_ensemble <- confusionMatrix(data = factor(y_hat_ensemble), reference = factor(test_set$Disease), 
                          positive = "Yes", mode="everything")
cm_ensemble
cm_ensemble_tab <- as.matrix(cm_ensemble, what = "classes")
cm_ensemble_tab

# Results
results <- cbind("Logistic" = cm_glmk_tab, "kNN" = cm_knn_tab, "RF" = cm_rf_tab, "Earth" =
                   cm_e_tab, "Ensemble" = cm_ensemble_tab)
colnames(results) <- c("Logistic", "kNN", "Random Forest", "Earth", "Ensemble")
results

# ROC, MCC table
auc.mcc <- matrix(c(auc1,auc2,auc3,auc4,auc5,mcc1,mcc2,mcc3,mcc4,mcc5), ncol = 5, byrow = TRUE)
rownames(auc.mcc) <- c("AUC", "MCC")
colnames(auc.mcc) <- c("Logistic", "kNN", "Random Forest", "Earth", "Ensemble")
auc.mcc <- as.table(auc.mcc)
auc.mcc
merged.results <- rbind(auc.mcc, results)
merged.results

# table intro
cm.table <- structure(c("True Positive (TP)", "False Negative (FN)", "False Positive (FP)", 
                     "True Negative (TN)"), dim = c(2,2), dimnames=structure(list(Predicted=c("Positive", "Negative"),
                                            Actual=c("Positive","Negative"))),
                                            rownames="Predicted",colnames="Actual", class="table")
cm.table