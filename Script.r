library(ggplot2)
library(psych)
library(stringr)
library(neuralnet)
library(stats)
library(tidyverse)

#********************************Data Loading********************************

colNames = c ("Gender", "Symptoms", "Alcohol", "Hepatitis B Surface Antigen",
              "Hepatitis B e Antigen", "Hepatitis B Core Antibody",
              "Hepatitis C Virus Antibody", "Cirrhosis", "Endemic Countries",
              "Smoking", "Diabetes", "Obesity", "Hemochromatosis",
              "Arterial Hypertension", "Chronic Renal Insufficiency",
              "Human Immunodeficiency Virus", "Nonalcoholic Steatohepatitis",
              "Esophageal Varices", "Splenomegaly", "Portal Hypertension",
              "Portal Vein Thrombosis", "Liver Metastasis", "Radiological Hallmark",
              "Age at diagnosis", "Grams of Alcohol per day", "Packs of cigarets per year",
              "Performance Status", "Encefalopathy degree", "Ascites degree",
              "International Normalised Ratio", "Alpha-Fetoprotein (ng/mL)",
              "Haemoglobin (g/dL)", "Mean Corpuscular Volume (fl)","Leukocytes(G/L)",
              "Platelets (G/L)", "Albumin (mg/dL)", "Total Bilirubin(mg/dL)",
              "Alanine transaminase (U/L)", "Aspartate transaminase (U/L)",
              "Gamma glutamyl transferase (U/L)", "Alkaline phosphatase (U/L)",
              "Total Proteins (g/dL)", "Creatinine (mg/dL)", "Number of Nodules",
              "Major dimension of nodule (cm)", "Direct Bilirubin (mg/dL)", "Iron	(mcg/dL)",
              "Oxygen Saturation (%)", "Ferritin (ng/mL)", "Class")

data = read.delim("hcc-data.txt", header = FALSE, sep = ",",
                  col.names = colNames, na.strings = "?", dec = ".")

cat(sum(is.na(data)), " values are missing.")

#***********************Preprocessing & Transformation***********************

isBinary <- function(v) {
  x <- na.omit(unique(v))
  return (length(x)==2)
}

#Filling missing values
clean_data = data
for (i in c(1:nrow(clean_data)))
  for (j in c(1:ncol(clean_data)))
    if (is.na(clean_data[i,j])){
      subData = clean_data[clean_data$Class == clean_data[i,ncol(data)],]
      if(isBinary(data[,j]))
       clean_data[i,j] = round(mean(subData[,j], na.rm = TRUE),0)
      else
       clean_data[i,j] = mean(subData[,j], na.rm = TRUE)
    }
clean_data = round(clean_data,2)

#Check if all missing values were filled
sum(is.na(clean_data))

#Some visualizations
hist(clean_data$Class, col = "red", xlab = "Class", main = "Histogram of Classes")

#Scatter plot of matrices (SPLOM)
df = data.frame(clean_data$Gender, clean_data$Direct.Bilirubin..mg.dL., clean_data$Total.Bilirubin.mg.dL.)
df = df %>%
  rename(
    Gender = clean_data.Gender,
    Direct_Bilirubin = clean_data.Direct.Bilirubin..mg.dL.,
    Total_Bilirubin = clean_data.Total.Bilirubin.mg.dL.
  )
pairs.panels(df, lm = TRUE)

boxplot(clean_data, col = "green")

getmode <- function(v) {
  uniqv <- unique(v)
  occurences <- tabulate(match(v, uniqv))
  uniqv[which(occurences == max(occurences))]
}

cat("Mode of the feature 'Packs of cigarettes per year':",
    getmode(clean_data[clean_data$Class == 0,]$Packs.of.cigarets.per.year))
summary(clean_data[clean_data$Class == 0,]$Packs.of.cigarets.per.year)

#Min-Max normalization
max = apply(clean_data, 2, max)
min = apply(clean_data, 2, min)
clean_data = as.data.frame(scale(clean_data, center = min, scale = max - min))
boxplot(clean_data, col = "green")

#Check if all data are normalized (between 0 and 1)
cat("Min =", min(clean_data))
cat("Max =", max(clean_data))

#Detect and remove outliers
correlations = round(cor(clean_data),2)
heatmap(correlations) 
diag(correlations) = 0
output_corrs = correlations[ncol(clean_data),]
impactful_var_index = which(output_corrs == max(output_corrs))
boxplot(clean_data[,impactful_var_index], col = "salmon")
sum_stats = summary(clean_data[,impactful_var_index])
no_big_outliers = length(which(clean_data[,impactful_var_index]>sum_stats[5]+1.5*(sum_stats[5]-sum_stats[2])))
no_small_outliers = length(which(clean_data[,impactful_var_index]<sum_stats[2]-1.5*(sum_stats[5]-sum_stats[2])))

#Smoothing by bin mean
sorted_impactful_var = sort(clean_data[,impactful_var_index])
indexes = sort(clean_data[,impactful_var_index],index.return =TRUE)$ix
bins_means = round(ave(sorted_impactful_var,
                   rep(1:length(sorted_impactful_var), each = 3,
                   length.out = length(sorted_impactful_var))),digits=1)
#Replacing outliers above Q3 + 1.5*IQR
for (o1 in (length(indexes)-no_big_outliers+1):(length(indexes))) {
  clean_data[indexes[o1],impactful_var_index]=bins_means[o1]
}
#Replacing outliers below Q1 - 1.5*IQR
for (o2 in (1:no_small_outliers)) {
  clean_data[indexes[o2],impactful_var_index]=bins_means[o2]
}

boxplot(clean_data[,impactful_var_index], col = "salmon")

#Look for the most correlated variables to reduce their numerosity
correlations[!lower.tri(correlations)] = 0
highest_cor = which(correlations >= 0.9, arr.ind = TRUE) #Threshold set to 90%
toDrop = c()
for (k in c(1:nrow(highest_cor))){
  if(abs(correlations[nrow(correlations),highest_cor[k,1]]) < abs(correlations[ncol(correlations),highest_cor[k,2]]))
    toDrop = c(toDrop, colNames[highest_cor[k,1]])
  else
    toDrop = c(toDrop, colNames[highest_cor[k,2]])
}

#Use regex to solve the features' misnomers (Special characters are stored as dots in R)
toDrop = str_replace_all(toDrop, c("\\("=".", "\\)"=".", " "=".", "/"="."))
clean_data = clean_data[,!(names(clean_data) %in% toDrop)]

#Remove duplicated records
clean_data = unique(clean_data)

#Principal Component Analysis
data_pca = prcomp(clean_data, center = TRUE, scale. = TRUE)
standard_deviation = data_pca$sdev^2
standard_deviation_percentage = round(standard_deviation/sum(standard_deviation)*100,0)
#Principal Components
barplot(standard_deviation_percentage, xlab = "Principal Component",
        ylab = "variation", main = "Scree Plot")
scores = abs(data_pca$rotation[,1])
ranking = names(sort(scores, decreasing = TRUE))

#********************************Classification********************************

sampleSize = 0.7 * nrow(clean_data)
set.seed(1) #Equivalent to "random_state" in Python
index = sample(seq_len(nrow(clean_data)), size = sampleSize )
data_train = clean_data[index,]
hist(data_train$Class, xlab = "Class", main = "Histogram of the training data's Classes", col = "red")
data_test = clean_data[-index,]
correlations = round(cor(clean_data),2) #We recalculate the correlations in case a feature has been dropped
most_correlated = sort(correlations[nrow(correlations),], decreasing = TRUE)
NN = neuralnet(Class ~ Iron..mcg.dL.+Smoking+Haemoglobin..g.dL.+Albumin..mg.dL., data = data_train, hidden = 3, threshold = 0.01, linear.output = TRUE)
plot(NN)
test_inputs = data_test[,1:ncol(data_test)-1]
actual_outputs = data_test[,ncol(data_test)]
predicted_outputs = round(neuralnet::compute(NN, test_inputs)$net.result,0)

#View the actual and the predicted values to spot the missclassified records
data.frame(actual_outputs, predicted_outputs)

#Confusion matrix
confusion = table(actual_outputs, predicted_outputs)

accuracy = (confusion[1,1]+confusion[2,2])/length(predicted_outputs)
error_rate = 1 - accuracy
sensitivity = confusion[2,2]/(confusion[2,1]+confusion[2,2])
specificity = confusion[1,1]/(confusion[1,1]+confusion[1,2])
