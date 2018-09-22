#Loading necessary Libraries
library(Matching)
library(rbounds)
library(rgenoud)
library(quantreg)

#Load the data
nsw_dw <- read_dta("C:/Users/fabia/Desktop/nsw_dw.dta")
cps_controls <- read_dta("C:/Users/fabia/Desktop/cps_controls.dta")


#1 Getting difference of means and confidence interval of initial nsw_dw data
treatdata <- subset(nsw_dw, treat == 1)
controldata <- subset(nsw_dw, treat ==0)

#Calculating difference of means between treatment and control
diffofmean <- mean(treatdata$re78) - mean(controldata$re78)
diffofmean

#Running univariate regression to obtain confidence interval
glm1 <- glm(nsw_dw$re78 ~ nsw_dw$treat)
summary(glm1)
confint(glm1)


#Problem 2: Obtaining the difference of means and confidence interval of NSw_dw treat data with cps controls
cpscontroldata <- cps_controls

#Calculating difference of means between treatment and control
meandiffwcps <- mean(treatdata$re78) - mean(cpscontroldata$re78)
meandiffwcps

#Running univariate regression to obtain confidence interval
new_nsw_wcps <- rbind(data.frame(treatdata), data.frame(cpscontroldata))
glm2 <- glm(new_nsw_wcps$re78 ~ new_nsw_wcps$treat)
summary(glm2)
confint(glm2)

#3 Using propensity score matching to obtain a treatment effect and confidence interval
#Creating the propensity score model
glm3 <- glm(treat ~ age+ education +black + hispanic + married + nodegree + re74 + re75 , family = binomial, data = new_nsw_wcps)
X  <- glm3$fitted
Y  <- new_nsw_wcps$re78
Tr  <- new_nsw_wcps$treat

#Matching 
rr  <- Match(Y=Y, Tr=Tr, X=X, M=1, ties = TRUE, estimand = "ATT");
summary(rr)

#Checking the match balance
MatchBalance(treat~age + education + black +
               hispanic + married + nodegree + re74  + re75, data = new_nsw_wcps, match.out=rr, nboots=10)

#Obtaining the confidence interval and treatment effect a univariate linear regression
glm3.1 <- glm(rr$mdata$Y ~ rr$mdata$Tr )
summary(glm3.1)
confint(glm3.1)



#4 Running a multivariate matching procedure using all the covariates and estimated propensity scores
Y2 <- new_nsw_wcps$re78
Tr2 <- new_nsw_wcps$treat
X2 <- cbind(new_nsw_wcps$age,
            new_nsw_wcps$education,
            new_nsw_wcps$black == 1,
            new_nsw_wcps$black == 0,
            new_nsw_wcps$hispanic == 1,
            new_nsw_wcps$hispanic == 0,
            new_nsw_wcps$married == 1,
            new_nsw_wcps$married == 0,
            new_nsw_wcps$nodegree == 1,
            new_nsw_wcps$nodegree == 0,
            new_nsw_wcps$re74,
            new_nsw_wcps$re75,
            glm3$fitted)

#Matching
mout.mv <- Match(Y=Y2, Tr=Tr2, X=X2, M=1)
summary(mout.mv)

#Checking matching balannce
mbofmout.mv <-MatchBalance(treat~age + education + black +
               hispanic + married + nodegree + re74  + re75, data = new_nsw_wcps, match.out=mout.mv, nboots=10)

#Getting confidence intervals by running a linear regression and viewing its summary
glm4 <- glm(mout.mv$mdata$Y ~ mout.mv$mdata$Tr )
summary(glm4)
confint(glm4)


#5 Repetition of 3 and 4 with Genetic matching
#Repetition of 3 with genetic matching

#Running the genmatch
p.genout53 <- GenMatch(Tr=Tr, X=X, M=1, estimand="ATT", pop.size = 100, max.generations = 20, wait.generations = 7)

#Running the match function to obtain nearest neighbour match obtained by genmatch function
mout.genout53 <- Match(Y=Y, Tr=Tr, X=X, M=1, estimand="ATT", Weight.matrix = p.genout53)
summary(mout.genout53)

#Checking match balance
matchbalance53 <- MatchBalance(treat~age + education + black +
                                 hispanic + married + nodegree + re74  + re75, 
                               data=new_nsw_wcps, match.out=mout.genout53, nboots=500)

#Obtaining confidence interval by running a regression viewing summary
glm5_3 <- glm(mout.genout53$mdata$Y ~ mout.genout53$mdata$Tr )
summary(glm5_3)
confint(glm5_3)


#Repetition of 4 with Genetic Matching

#Running the genmatch
p.genout54 <- GenMatch(Tr=Tr2, X=X2, M=1, estimand="ATT", pop.size = 100, max.generations = 10, wait.generations = 7)

#Running the match function to obtain nearest neighbour match obtained by genmatch function
mout.genout54 <- Match( Tr=Tr2, X=X2, M=1, estimand="ATT",  Weight.matrix = p.genout54)

summary(mout.genout54)

#Checking the match balance
matchbalance54 <- MatchBalance(treat~age + education + black +
                                 hispanic + married + nodegree + re74  + re75, 
                               data=new_nsw_wcps, match.out=mout.genout54, nboots=500)

#Obtaining confidence interval by running a regression viewing summary
glm5_4 <- glm(mout.genout54$mdata$Y ~ mout.genout54$mdata$Tr )
summary(glm5_4)
confint(glm5_4)

#Calculating the treatment effect at different quartiles
n<- 0.5
Ynot <- mout.genout53$mdata
Y_t <- subset(Ynot, Tr == 1)
Y_c <- subset(Ynot, Tr == 0)
quartileatn <- quantile(Y_t$Y, n) - quantile(Y_c$Y, n)
quartileatn
