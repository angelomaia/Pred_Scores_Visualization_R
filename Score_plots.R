### ---------------------- ###
# Printing lines with symbols for statistical parameters
### ---------------------- ###

require(dplyr)
require(stringr)

# Setting wd and importing db ----
dir = ("K:/git_path/Pred_Scores_Visualization_R/")

test.Orig = readRDS(paste0(dir, "db_val.Orig.rda"))
test.MSC = readRDS(paste0(dir, "db_val.MSC.rda"))
test.SNV = readRDS(paste0(dir, "db_val.SNV.rda"))
test.SG = readRDS(paste0(dir, "db_val.SG.rda"))
test.DET = readRDS(paste0(dir, "db_val.DET.rda"))
test.CR = readRDS(paste0(dir, "db_val.CR.rda"))

models <- list.files(path=paste0(dir, 'models'), full.names = FALSE)

i = 2 # first column with elemental concentration
j = 7 # last column with elemental concentration

results = data.frame(1:9)
metrics = c("")

# Function to assess prediction results ----
goof <- function(observed, predicted){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[4] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                   MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)
  
  return(gf)
}

# forward loop for prediction and results ----
for (k in 1:length(models)){
  if (k <= 36){m = i+4}
  if (k <= 30){m = i+3}
  if (k <= 24){m = i+5}
  if (k <= 18){m = i+2}
  if (k <= 12){m = i+1}
  if (k <= 6){m = i}

  if (paste(str_sub(models[k], 4, -14)) == "Orig"){test = test.Orig}
  if (paste(str_sub(models[k], 4, -14)) == "CR"){test = test.CR}
  if (paste(str_sub(models[k], 4, -14)) == "DET"){test = test.DET}
  if (paste(str_sub(models[k], 4, -14)) == "SG"){test = test.SG}
  if (paste(str_sub(models[k], 4, -14)) == "SNV"){test = test.SNV}
  if (paste(str_sub(models[k], 4, -14)) == "MSC"){test = test.MSC}

    model <- readRDS(paste0(dir, "models/", paste(models[k])))
    results[k] = stats::predict(model, test[, j:length(test)])
  
    colnames(results)[k] = c(paste(str_sub(models[k], 1, -5)))
  
    metric = goof(observed = test[,m], 
                     predicted = results[,k])
    
    if (nchar(names(test[m])) == 2){
    metric$Variable = paste(str_sub(models[k], 1, 2))
    metric$Preproc = paste(str_sub(models[k], 4, -14))
    }
    if (nchar(names(test[m])) == 1){
      metric$Variable = paste(str_sub(models[k], 1, 1))
      metric$Preproc = paste(str_sub(models[k], 3, -14))
    }
    assign(paste("metric", str_sub(models[k], 1, -14), sep="."), metric)
    
    metrics = c(metrics, paste("metric", str_sub(models[k], 1, -14), sep="."))
}
rm(metric)

metricsList = list()

for (x in 2:length(metrics)) 
  { dat <- data.frame(get(metrics[x]))
  metricsList[[(x-1)]] <- dat
}
rm(x)

res <- do.call(rbind, metricsList)


head(res, 5)
res= res %>%
  select(c('Variable', 'Preproc'), everything())

db1 <- res
# Building db for plot ----

dbRAW <- db1 %>% filter(Preproc == 'Orig') %>% group_by(Variable) %>% slice(which.min(RMSE))
dbMSC <- db1 %>% filter(Preproc == 'MSC') %>% group_by(Variable) %>% slice(which.min(RMSE))
dbSG <- db1 %>% filter(Preproc == 'SG') %>% group_by(Variable) %>% slice(which.min(RMSE))
dbSNV <- db1 %>% filter(Preproc == 'SNV') %>% group_by(Variable) %>% slice(which.min(RMSE))
dbDET <- db1 %>% filter(Preproc == 'DET') %>% group_by(Variable) %>% slice(which.min(RMSE))
dbCR <- db1 %>% filter(Preproc == 'CR') %>% group_by(Variable) %>% slice(which.min(RMSE))

Elements <- c(dbRAW$Variable)

db <- as.data.frame(Elements)
db$RAW <- dbRAW$RPIQ
db$MSC <- dbMSC$RPIQ
db$SG <- dbSG$RPIQ
db$SNV <- dbSNV$RPIQ
db$DET <- dbDET$RPIQ
db$CR <- dbCR$RPIQ
db$Labels <- c(1:6)

# Setting fonts for plot ----
windowsFonts(
  B=windowsFont('Arial'),
  A=windowsFont('Times New Roman')
)

# RPIQ
# Plot ----
jpeg(paste0(dir, 'RPIQ_scores.jpeg'), 
     width = 16, height = 10, units = 'in', res = 300)

summary(db)

par(mar = c(6, 9, 6, 2), lwd = 3, lab = c(22,9,10)) # bottom, left, top, right
plot(db$Labels, db$RAW,
     type = 'l',
     cex = 1.2,
     xlab = NA,
     ylab = NA,
     ylim = c(-0.09, 4),
     pch = 16,
     cex.axis=2,
     las=2,
     col = 1,
     family = 'A',
     xaxt = "n",
     axes = FALSE)
title(ylab = expression(paste('RPIQ')),
      cex.lab=2.5, line = 5, family = 'A')
axis(2, ylim = c(-0.1, 1), las=1, family = 'A', 
     cex.axis=2, lwd=3)
axis(1, labels = as.character(db$Elements), 
     at = as.numeric(db$Labels), family = 'A', 
     cex.axis=2, las=1, lwd=3)
### Lines
lines(db$Labels, db$RAW, type = 'b',
      col = 1, pch = 15, cex=2)
lines(db$Labels, db$MSC, type = 'l',
      col = 'red')
lines(db$Labels, db$MSC, type = 'b',
      col = 'red', pch = 16, cex=2)
lines(db$Labels, db$CR, type = 'l',
      col = 6)
lines(db$Labels, db$CR, type = 'b',
      col = 6, pch = 0, cex=2)
lines(db$Labels, db$DET, type = 'l',
      col = 51)
lines(db$Labels, db$DET, type = 'b',
      col = 51, pch = 18, cex=3)
lines(db$Labels, db$SG, type = 'l',
      col = 77)
lines(db$Labels, db$SG, type = 'b',
      col = 77, pch = 25, cex=2, bg = 77)
lines(db$Labels, db$SNV, type = 'l',
      col = 31)
lines(db$Labels, db$SNV, type = 'b',
      col = 31, pch = 17, cex=2)

legend("top",                                       # Add legend to plot
       legend = c(colnames(db[,2:7])),
       col = c(1, "red", 6, 51, 77, 31),
       pch = c(15, 16, 0, 18, 25, 17),
       pt.bg = c(NA, NA, NA, NA, 77, NA),
       lty = 1,
       cex = 2,
       pt.cex = c(2, 2, 2, 3, 2, 2),
       ncol = 3)

dev.off()

#----
#
#
#