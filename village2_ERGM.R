
rm(list=ls())

#Install all relevant packages
#install.packages('statnet')
#install.packages('latticeExtra')

library(parallel)
library(statnet)

#detect number of cores available
ncores <- detectCores()

#Load in all Village 2 data

v2.coop <- as.matrix(read.csv('./inputs/village2_cooperation.csv', header = FALSE))
v2.coop.delp <- as.matrix(read.csv('./inputs/village2_cooperation_deleted.csv', header = FALSE))


v2.kin <- as.matrix(read.csv('./inputs/village2_kinship.csv', header = FALSE))
v2.kin.delp <- as.matrix(read.csv('./inputs/village2_kinship_deleted.csv', header = FALSE))
v2.attributes <- read.csv('./inputs/village2_attributes.csv', header = TRUE)
v2.attributes.delp <- read.csv('./inputs/village2_attributes_deleted.csv', header = TRUE)

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
										+ nodeicov('status') 
										+ nodeocov('status')
										+ nodeicov('age') 
										+ nodeocov('age')
										)
summary(v2.coop.model2)

#Full model presented in Supplementary materials

v2.coop.model3 <- ergm(v2.cooperation ~ edges 
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
summary(v2.coop.model3)


#Lets check the diagnostics and fit 

mcmc.diagnostics(v2.coop.model1)
v2.coop.model1.GoF <- gof(v2.coop.model1)
plot(v2.coop.model1.GoF)

mcmc.diagnostics(v2.coop.model2)
v2.coop.model1.GoF <- gof(v2.coop.model2)
plot(v2.coop.model2.GoF)

mcmc.diagnostics(v2.coop.model3)
v2.coop.model3.GoF <- gof(v2.coop.model3)
plot(v2.coop.model3.GoF)


#Format results

results <- v2.coop.model3
estimate <- results$coef
st.error <- sqrt(diag(results$covar))
norm.var <- estimate/st.error
pval <- 2*pnorm(abs(norm.var), lower.tail = F)

v2.status.coop.out <- data.frame(
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
write.csv(v2.status.coop.out, file = "./outputs/village2-ergm-outputs.csv")


#######################
# Model for Revisions
#######################

#Assess the difference between the model with imputed data vs. removing participants with NAs 

rownames(v2.coop.delp) <- v2.attributes.delp$internal_id
colnames(v2.coop.delp) <- v2.attributes.delp$internal_id
rownames(v2.kin.delp) <- v2.attributes.delp$internal_id
colnames(v2.kin.delp) <- v2.attributes.delp$internal_id

v2.cooperation.delp <- network(v2.coop.delp, directed = TRUE)
v2.kinship.delp <- network(v2.kin.delp, directed = FALSE)

v2.cooperation.delp %v% 'id' <- v2.attributes.delp$internal_id
v2.cooperation.delp %v% 'age' <- v2.attributes.delp$log_age
v2.cooperation.delp %v% 'income' <- v2.attributes.delp$income
v2.cooperation.delp %v% 'status' <- v2.attributes.delp$status
v2.cooperation.delp %v% 'strength' <- v2.attributes.delp$strength


#A comparable base model

v2.partial.model1 <- ergm(v2.cooperation.delp ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship.delp)
										)
summary(v2.partial.model1)

#Partial model with status & log age

v2.partial.model2 <- ergm(v2.cooperation.delp ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship.delp) 
										+ nodeicov('status') 
										+ nodeocov('status')
										+ nodeicov('age') 
										+ nodeocov('age')
										)
summary(v2.partial.model2)

#Full model presented in Supplementary materials

v2.partial.model3 <- ergm(v2.cooperation.delp ~ edges 
										+ mutual 
										+ gwesp(decay=0.1, fixed=FALSE) 
										+ edgecov(v2.kinship.delp) 
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
summary(v2.partial.model3)


#Lets check the diagnostics and fit 


mcmc.diagnostics(v2.partial.model3)
v2.partial.model3.GoF <- gof(v2.partial.model3)
plot(v2.partial.model3.GoF)

#format results

results <- v2.partial.model3
estimate <- results$coef
st.error <- sqrt(diag(results$covar))
norm.var <- estimate/st.error
pval <- 2*pnorm(abs(norm.var), lower.tail = F)

v2.revision.out <- data.frame(
                       estimate = round(estimate, 3),
                       st.error = round(st.error, 3),
                       norm.var = round(norm.var, 2),
                       pval = round(pval, 3),
                       OR = round(exp(estimate), 2),
                        CI.low = round(exp(estimate -
                                             qnorm(0.975)*st.error), 2),
                        CI.high = round(exp(estimate + qnorm(0.975)*st.error), 2)
                        )

v2.revision.out
write.csv(v2.revision.out, file = "./outputs/revision-village2-ergm-outputs.csv")


