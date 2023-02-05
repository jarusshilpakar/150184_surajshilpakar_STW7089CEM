## All the library used for this task. 
library(matlib)
library(ggplot2)
library(rsample)
# After importing your data

# Importing X data
X <- read_csv("Desktop/softwarica MSC/X.csv")
x = data.matrix(X)
X = as.matrix(read.csv(file="Desktop/softwarica MSC/X.csv", header = F))
colnames(X)<-c("X1","X2","X3","X4")

# Importing Y data 
Y = as.matrix(read.csv(file="Desktop/softwarica MSC/y.csv", header = F)) 
colnames(Y)<-c("Y")
# Importing Time data
Time = read.csv("Desktop/softwarica MSC/time.csv", header = F, skip = 1) 
Time = as.matrix(rbind(0, Time))

## Task 1: Creating a time series plot
#Defining the value of X and Y against time for time series plot. 
ts.X<-ts(X,start = c(min(Time),max(Time)), frequency = 1) 
ts.Y<-ts(Y,start = c(min(Time),max(Time)), frequency = 1)
#Plotting those time series graphs
plot(ts.X, main = "Time series plot of (Input) X EEG Signal", xlab = "Time", ylab = "Input EEG signal", col = "darkcyan") 
plot(ts.Y, main = "Time series plot of Y (Output) EEG Signal", xlab = "Time", ylab = "Output EEG signal", col = "darkcyan")

##Task 1: Creating a density plot with the help of histogram
#Creating a density if X signal
den = density(X)
plot(den, main = "Density plot of input EEG signal", col = "darkcyan") 

# Creating a Histogram of X signal
hist(X, freq = FALSE, main = "Histogram and Density plot of input EEG signal")

#Adding density in the histogram 
lines(den, lwd=2, col="darkcyan")
rug(jitter(X))

#Creating a density if X1 signal
den_X1 = density(X[,"X1"])
hist(X[,"X1"], freq = FALSE, main = "Histogram and density plot of X1 Signal", xlab = "X1 Signal") 
lines(den_X1, lwd=2, col="darkcyan")

# Add the data-poins with noise in the X-axis
rug(jitter(X[,"X1"]))

#Creating a density if X2 signal
den_X2 = density(X[,"X2"])
hist(X[,"X2"], freq = FALSE, main = "Histogram and density plot of X2 Signal", xlab = "X2 Signal")
lines(den_X2, lwd=2, col="darkcyan")
rug(jitter(X[,"X2"]))

#Creating a density if X3 signal
den_X3 = density(X[,"X3"])
hist(X[,"X3"], freq = FALSE, main = "Histogram and density plot of X3 Signal", xlab = "X3 Signal")
lines(den_X3,lwd=2,col="darkcyan")
rug(jitter(X[,"X3"]))

#Creating a density if X4 signal
den_X4 = density(X[,"X4"])
hist(X[,"X4"], freq = FALSE, main = "Histogram and density plot of X4 Signal", xlab = "X4 Signal") 
lines(den_X1, lwd=2, col="darkcyan")
rug(jitter(X[,"X4"]))

#Creating a density if Y signal
den_y = density(Y)
plot(den_y, main = "Density plot of Y (output) EEG signal", xlab = "Output Signal", col = "darkcyan")
hist(Y, freq = FALSE, main = "Histogram and density plot of Y Signal", xlab = "Output Signal") 
lines(den_y, lwd=2, col="darkcyan")
rug(jitter(Y))

# Plotting X1 against Y
plot(X[,"X1"], Y, main = "Correlation betweeen X1 and Y signal", xlab = "X1 signal", ylab = "Output signal", col = "darkcyan" )
# Plotting X2 against Y
plot(X[,"X2"],Y,main = "Correlation betweeen X2 and Y signal", xlab = "X2 signal", ylab = "Output signal", col = "darkcyan" )
# Plotting X3 against Y
plot(X[,"X3"],Y,main = "Correlation betweeen X3 and Y signal", xlab = "X3 signal", ylab = "Output signal", col = "darkcyan" )
# Plotting X4 against Y
plot(X[,"X4"],Y,main = "Correlation betweeen X4 and Y signal", xlab = "X4 signal", ylab = "Output signal", col = "darkcyan" )

#Task 2: Regression – modeling the relationship between EEG signals 
#Task 2.1
# Calculating ones for binding the data 
ones = matrix(1 , length(X)/4, 1)
ones

#Calculating thetahat of Model 1
#Binding data from equation of model 1.
X_M1 <- cbind(ones,X[,"X4"],X[,"X1"]^2,X[,"X1"]^3,X[,"X2"]^4,X[,"X1"]^4) 
X_M1
#Calculating thetahat of model 1
thetahat_M1 = solve(t(X_M1) %*% X_M1) %*% t(X_M1) %*% Y
thetahat_M1

#Model 2
#Binding data from equation of model 2. 
X_M2 <- cbind(ones,X[,"X4"],X[,"X1"]^3,X[,"X3"]^4)
X_M2
#Calculating thetahat of Model 2
thetahat_M2 = solve(t(X_M2) %*% X_M2) %*% t(X_M2) %*% Y 
thetahat_M2

#For model 3
#Binding data from equation of model 3. 
X_M3 <- cbind(ones, (X[,"X3"])^3, (X[,"X3"])^4)
X_M3
#Calculating thetahat of Model 3
thetahat_M3 = solve(t(X_M3) %*% X_M3) %*% t(X_M3) %*% Y
thetahat_M3

#For model 4
#Binding data from equation of model 4. 
X_M4 <- cbind(ones,X[,"X2"],(X[,"X1"])^3,(X[,"X3"])^4)
X_M4
#Calculating thetahat of Model 4
thetahat_M4 = solve(t(X_M4) %*% X_M4) %*% t(X_M4) %*% Y 
thetahat_M4

#Binding data from equation of model 5. 
X_M5<-cbind(ones,(X[,"X4"]),(X[,"X1"])^2,(X[,"X1"])^3,(X[,"X3"])^4)
X_M5
#Calculating thetahat of Model 5
thetahat_M5=solve(t(X_M5) %*% X_M5) %*% t(X_M5) %*% Y 
thetahat_M5

#Task 2.2
#Calculating Yhat and RSS Model1
Yhat_M1 = X_M1 %*% thetahat_M1 
Yhat_M1
#Calculating RSS for model1
RSS_M1 = sum((Y-Yhat_M1)^2) 
RSS_M1

#Calculating Yhat and RSS Model2
Yhat_M2 = X_M2 %*% thetahat_M2
Yhat_M2
#Calculating RSS for model2
RSS_M2 = sum((Y-Yhat_M2)^2) 
RSS_M2

#Calculating Yhat and RSS Model3
Yhat_M3 = X_M3 %*% thetahat_M3 
Yhat_M3
#Calculating RSS for model3
RSS_M3 = sum((Y-Yhat_M3)^2) 
RSS_M3

#Calculating Yhat and RSS Model4
Yhat_M4 = X_M4 %*% thetahat_M4 
Yhat_M4
#Calculating RSS for model4
RSS_M4 = sum((Y-Yhat_M4)^2) 
RSS_M4

#Calculating Yhat and RSS Model5 
Yhat_M5 = X_M5 %*% thetahat_M5 
Yhat_M5
#Calculating RSS for model5
RSS_M5 = sum((Y-Yhat_M5)^2) 
RSS_M5

### Task 2.3 Calculating likelihood and Variance of each model
N = length(Y)
#Calculating the Variance of Model1
Variance_M1 = RSS_M1/(N-1)
Variance_M1
#Calculating the log-likelihood of Model1
likehood_M1 = -(N/2)*(log(2*pi))-(N/2)*(log(Variance_M1))-(1/(2*Variance_M1))*RSS_M1 
likehood_M1

#Calculating the Variance of Model2
Variance_M2 = RSS_M2/(N-1)
Variance_M2
#Calculating the log-likelihood of Model2
likehood_M2 = -(N/2)*(log(2*pi))-(N/2)*(log(Variance_M2))-(1/(2*Variance_M2))*RSS_M2 
likehood_M2

#Calculating the Variance of Model3
Variance_M3 = RSS_M3/(N-1)
Variance_M3
#Calculating the log-likelihood of Model2
likehood_M3 = -(N/2)*(log(2*pi))-(N/2)*(log(Variance_M3))-(1/(2*Variance_M3))*RSS_M3 
likehood_M3

#Calculating the Variance of Model4
Variance_M4 = RSS_M4/(N-1)
Variance_M4
#Calculating the log-likelihood of Model2
likehood_M4 = -(N/2)*(log(2*pi))-(N/2)*(log(Variance_M4))-(1/(2*Variance_M4))*RSS_M4
likehood_M4

#Calculating the Variance of Model5
Variance_M5 = RSS_M5/(N-1)
Variance_M5
#Calculating the log-likelihood of Model3
likehood_M5 = -(N/2)*(log(2*pi))-(N/2)*(log(Variance_M5))-(1/(2*Variance_M5))*RSS_M5
likehood_M5


### Task 2.4 Calculating AIC And BIC of Each model
##Calculating AIC and BIC of model 1 
K1 <- length(thetahat_M1)
K1 
AIC_M1 = 2*K1-2*likehood_M1
AIC_M1 
BIC_M1 = K1*log(N)-2*likehood_M1 
BIC_M1

##Calculating AIC and BIC of model 2
K2 <- length(thetahat_M2)
K2 
AIC_M2 = 2*K2-2*likehood_M2
AIC_M2 
BIC_M2 = K2*log(N)-2*likehood_M2 
BIC_M2

##Calculating AIC and BIC of model 3
K3 <- length(thetahat_M3)
K3 
AIC_M3 = 2*K3-2*likehood_M3
AIC_M3 
BIC_M3 = K3*log(N)-2*likehood_M3 
BIC_M3

##Calculating AIC and BIC of model 4
K4 <- length(thetahat_M4)
K4 
AIC_M4 = 2*K4-2*likehood_M4
AIC_M4 
BIC_M4 = K4*log(N)-2*likehood_M4 
BIC_M4

##Calculating AIC and BIC of model 5
K5 <- length(thetahat_M5)
K5 
AIC_M5 = 2*K5-2*likehood_M5
AIC_M5 
BIC_M5 = K5*log(N)-2*likehood_M5 
BIC_M5

## Task 2.5
## Error of model1
error_M1 <- Y-Yhat_M1
## Plotting the graph QQplot and QQ line of model 1 
qqnorm(error_M1, col = "blue", main = "QQ plot of model 1") 
qqline(error_M1, col = "red", lwd=1)

## Error of model2
error_M2 <- Y-Yhat_M2
## Plotting the graph QQplot and QQ line of model 2
qqnorm(error_M2, col = "blue", main = "QQ plot of model 2") 
qqline(error_M2, col = "red", lwd=1)

## Error of model3
error_M3 <- Y-Yhat_M3
## Plotting the graph QQplot and QQ line of model 3
qqnorm(error_M3, col = "blue", main = "QQ plot of model 3") 
qqline(error_M3, col = "red", lwd=1)

## Error of model4
error_M4 <- Y-Yhat_M4
## Plotting the graph QQplot and QQ line of model 4
qqnorm(error_M4, col = "blue", main = "QQ plot of model 4") 
qqline(error_M4, col = "red", lwd=1)

## Error of model5
error_M5 <- Y-Yhat_M5
## Plotting the graph QQplot and QQ line of model 5
qqnorm(error_M5, col = "blue", main = "QQ plot of model 5") 
qqline(error_M5, col = "red", lwd=1)


### Task 2.7
## Spliting the data of y into 2 form i.e. Traning and testing data set. 
Y_split <- initial_split(data = as.data.frame(Y), prop=.7)
## Traning Y data split 
training_dataset_Y<-training(Y_split) 
testing_dataset_Y<-as.matrix(testing(Y_split)) 
## Testing Y data split 
training_data_Y<-as.matrix(Y_training_set)



## Spliting the data of x into 2 form i.e. Traning and testing data set. 
X_split <- initial_split(data = as.data.frame(X), prop=.7)
## Traning X data split
training_dataset_X <- training(X_split)
## Testing X data split 
testing_dataset_X<-as.matrix(testing(X_split)) 
testing_data_X<-as.matrix(testing_dataset_X) 
traning_data_X<-as.matrix(training_dataset_X)

### Estimating model parameters using Training set
ones_training = matrix(1 , length(training_dataset_X$X1),1) 
traning_model_X <- cbind(ones_training, training_dataset_X[,"X4"], (training_dataset_X[,"X1"])^3, (training_dataset_X[ ,"X3"])^4)
traning_thetahat = solve(t(traning_model_X) %*% traning_model_X) %*% t(traning_model_X) %*% training_data_Y

### Model out/Prediction
testing_hat_Y = testing_data_X %*% traning_thetahat 
testing_hat_Y 
RSS_testing = sum((testing_dataset_Y-testing_hat_Y)^2) 
RSS_testing
t.test(training_data_Y, mu=500, alternative="two.sided", conf.level=0.95)
C_I1 = -0.1747735 
C_I2 = 0.4326150

p2 <- plot(density(training_data_Y), col="darkcyan", lwd=2, main="Distribution of Training Data")
abline(v=C_I1, col="red", lty=2) 
abline(v=C_I2, col="red", lty=2)
thetaHat_training = solve(t(traning_data_X) %*% traning_data_X) %*% t(traning_data_X) %*% training_data_Y
thetaHat_training
length(thetaHat_training)
dis_test = density(training_data_Y) 
plot((dis_test))
plot(dis_test,main = "Density plot of Y Signal")

### Calculating Confidential interval 
z=1.96 ##(95%) Confidential interval
error = ((testing_dataset_Y - testing_hat_Y)) 
error
n_len = length(testing_hat_Y)
C_I_1 = z * sqrt( (error * (1-error) ) / n_len) 
C_I_1
error
C_I_2 = z*sqrt((error*(1+error)/n_len))
C_I_2

#ploting error bars for X1 data
SD = sqrt(Variance_M1)
plot_data = data.frame(XA_Value = X_M1, YA_Value = Y)

ggplot(plot_data) +
 geom_bar( aes(x=XA_Value.1, y=Y), stat="identity", fill="red", alpha=0.7) +
  geom_errorbar( aes(x=XA_Value.1, ymin=Y-SD, ymax=Y+SD), width=0.4, colour="darkcyan", alpha=0.9, linewidth=1)

#ploting error bars for X2 data
SD = sqrt(Variance_M2)
plot_data = data.frame(XA_Value = X_M2, YA_Value = Y)

ggplot(plot_data) +
  geom_bar( aes(x=XA_Value.2, y=Y), stat="identity", fill="red", alpha=0.7) +
  geom_errorbar( aes(x=XA_Value.2, ymin=Y-SD, ymax=Y+SD), width=0.4, colour="darkcyan", alpha=0.9, linewidth=1)

#ploting error bars for X3 data
SD = sqrt(Variance_M3)
plot_data = data.frame(XA_Value = X_M3, YA_Value = Y)

ggplot(plot_data) +
  geom_bar( aes(x=XA_Value.3, y=Y), stat="identity", fill="red", alpha=0.7) +
  geom_errorbar( aes(x=XA_Value.3, ymin=Y-SD, ymax=Y+SD), width=0.4, colour="darkcyan", alpha=0.9, linewidth=1)

#ploting error bars for X4 data
SD = sqrt(Variance_M4)
plot_data = data.frame(XA_Value = X_M4, YA_Value = Y)

ggplot(plot_data) +
  geom_bar( aes(x=XA_Value.4, y=Y), stat="identity", fill="red", alpha=0.7) +
  geom_errorbar( aes(x=XA_Value.4, ymin=Y-SD, ymax=Y+SD), width=0.4, colour="darkcyan", alpha=0.9, linewidth=1)

##Task 3
## Model 2 will be used, parameter are selected and kept constant. 
array_1 = 0
array_2 = 0
f_value = 0
s_value = 0
thetahat_M2
#values from thetahat
thetebias <- 0.483065688 #choosen parameter
thetaone <- 0.143578928 # chosen prarameter
thetatwo <- 0.010038614 # constant value
thetafour <- -0.001912836 # constant value
elison <- RSS_M2 * 2 ## fixing value of elision 
num <- 100 #number of iteration

##Calculating Y-hat for performing rejection ABC 
counter <- 0
for (i in 1:num) {
  range1 <- runif(1,-0.483065688, 0.483065688) # calculating the range 
  range1
  range2 <- runif(1,-0.143578928, 0.143578928)
  new_thetahat <- matrix(c(range1,range2,thetatwo,thetafour)) 
  new_Y_Hat <- X_M2 %*% new_thetahat ## New Y hat
  new_RSS <- sum((Y-new_Y_Hat)^2) 
  new_RSS
  if (new_RSS > elison){
    array_1[i] <- range1 
    array_2[i] <- range2 
    counter = counter+1
    f_value <- matrix(array_1) 
    s_value <- matrix(array_2)
  } 
}
hist(f_value)
hist(s_value)
###ploting the graph
plot(f_value, s_value, col = c("red", "darkcyan"), main = "Joint and Marginal Posterior Distribution")
