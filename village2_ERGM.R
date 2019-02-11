
rm(list=ls())

#Install all relevant packages
#install.packages('statnet')
#install.packages('latticeExtra')

library(statnet)
#library(latticeExtra)


#Load in all Village 2 data

v2.coop <- as.matrix(read.csv('./inputs/village2_cooperation.csv', header = FALSE))
v2.kin <- as.matrix(read.csv('./inputs/village2_kinship.csv', header = FALSE))
v2.attributes <- read.csv('./inputs/village2_attributes.csv', header = TRUE)

rownames(v2.coop) <- v2.attributes$pid
colnames(v2.coop) <- v2.attributes$pid
rownames(v2.kin) <- v2.attributes$pid
colnames(v2.kin) <- v2.attributes$pid

v2.cooperation <- network(v2.coop, directed = TRUE)
v2.kinship <- network(v2.kin, directed = FALSE)

#Take a look at the descriptive statistics for village 2

plot(v2.coop)
gden(v2.coop)
isolates(v2.coop)
mean(degree(v2.coop))

plot(v2.kinship)
gden(v2.kinship)
isolates(v2.kinship)
mean(degree(v2.kinship))

#Specify the attributes for village 2

v2.cooperation %v% 'id' <- v2.attributes$pid
v2.cooperation %v% 'age' <- v2.attributes$log_age
v2.cooperation %v% 'income' <- v2.attributes$income
v2.cooperation %v% 'status' <- v2.attributes$status
v2.cooperation %v% 'strength' <- v2.attributes$strength


#Check that the attributes are there
list.vertex.attributes(v2.cooperation)


#Now lets run the models

#A base model

v2.coop.model1 <- ergm(v2.cooperation ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship)
										)
summary(v2.coop.model1)

#Partial model with status & log age

v2.coop.model2 <- ergm(v2.cooperation ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship) 
										+ nodeicov('age') 
										+ nodeocov('age')
										)
summary(v2.coop.model2)

#Full model presented in Supplementary materials

v2.coop.model6 <- ergm(v2.cooperation ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship) 
										+ nodeicov('status') 
										+ nodeocov('status') 
										+ nodeicov('income') 
										+ nodeocov('income') 
										+ nodeicov('strength')  
										+ nodeocov('strength')  
										+ nodeicov('age') 
										+ nodeocov('age'),
										control = control.ergm(MCMC.burnin = 20000, 
														MCMC.samplesize = 40000,
														MCMC.interval = 800) 
										)
summary(v2.coop.model6)


#Lets check the diagnostics and fit 

mcmc.diagnostics(v2.coop.model1)
v2.coop.model1.GoF <- gof(v2.coop.model1)
plot(v2.coop.model1.GoF)

mcmc.diagnostics(v2.coop.model2)
v2.coop.model1.GoF <- gof(v2.coop.model2)
plot(v2.coop.model2.GoF)

mcmc.diagnostics(v2.coop.model6)
v2.coop.model6.GoF <- gof(v2.coop.model6)
plot(v2.coop.model6.GoF)

results <- v2.coop.model6
parameter <- results$effects$effectName
estimate <- results$theta
st.error <- sqrt(diag(results$covtheta)
norm.var <- estimate/st.error
pval <- 2*pnorm(abs(norm.var), lower.tail = F)

v2.status.coop.out <- data.frame(parameter,
                       estimate = round(estimate, 3),
                       st.error = round(st.error, 3),
                       norm.var = round(norm.var, 2),
                       pval = round(pval, 3),
                       OR = round(exp(estimate), 2),
                        CI.low = round(exp(estimate -
                                             qnorm(0.975)*st.error), 2),
                        CI.high = round(exp(estimate + qnorm(0.975)*st.error), 2)
                        )

v2.status.coop.out
write.csv(v2.status.coop.out, file = "village2-ergm-outputs.csv")



