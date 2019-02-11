
rm(list=ls())

#Install all relevant packages
#install.packages('Rsiena')

library(RSiena)

#Read in cooperation networks

coop1 <- as.matrix(read.csv("./inputs/village1_cooperation_t1.csv",
	 header = FALSE, stringsAsFactors = FALSE))
coop2 <- as.matrix(read.csv("./inputs/village1_cooperation_t2.csv",
	 header = FALSE, stringsAsFactors = FALSE))
coop3 <- as.matrix(read.csv("./inputs/village1_cooperation_t3.csv", 
	header = FALSE, stringsAsFactors = FALSE))

#Read in kinship networks

kin1 <- as.matrix(read.csv("./inputs/village1_kinship_t1.csv", 
	header = FALSE, stringsAsFactors = FALSE))
kin2 <- as.matrix(read.csv("./inputs/village1_kinship_t2.csv",	
	header = FALSE, stringsAsFactors = FALSE))

#Read in attributes

att <- as.matrix(read.csv('./inputs/village1_attributes.csv',
	header = TRUE, stringsAsFactors = FALSE))

#Combine & specify networks

coop <- sienaNet(array( c(coop1,coop2, coop3), dim = c(74, 74, 3)))
cooperation <- sienaDependent(coop)

kin <- sienaNet(array( c(kin1, kin2), dim = c(74, 74, 2)))
kinship <- varDyadCovar(kin)

#Specify the attributes & behaviors

status <- sienaDependent(att[ , 2:4], type = 'behavior')
income <- varCovar(att[ , 5:6])
log.age <- varCovar(att[ , 7:8])
strength <- varCovar(att[ , 9:10])

#Specify changing network composition

comp <- rep(list(c(1,3)), 74)

comp[[20]] <- c(1,2.5)
comp[[21]] <- c(1,2.5)
comp[[52]] <- c(1,2.5)
comp[[72]] <- c(1,2.5)
comp[[5]] <- c(2,3)
comp[[8]] <- c(2,3)
comp[[12]] <- c(2,3)
comp[[32]] <- c(2,3)
comp[[41]] <- c(2,3)
comp[[42]] <- c(2,3)
comp[[43]] <- c(2,3)
comp[[45]] <- c(2,3)
comp[[51]] <- c(2,3)
comp[[55]] <- c(2,3)
comp[[59]] <- c(2,3)
comp[[70]] <- c(2,3)
comp[[71]] <- c(2,3)
comp[[74]] <- c(2,3)

comp.change <- sienaCompositionChange(comp)

#Create the Siena dataset

cooperation.data <- sienaDataCreate(cooperation, 
									kinship, 
									status, 
									strength,
									income, 
									log.age, 
									comp.change
									)

#Check the descriptive statistics

print01Report(data = cooperation.data,
	 modelname = 'cooperation.status.out')

print("Data specified and network descriptives available ")

#All looks ok

#Create the algorithms 

alg1 <- sienaAlgorithmCreate(projname = 'cooperation_status', 
								n3 = 1000, 
								seed = 54321, 
								doubleAveraging = 0, 
								diagonalize = 0.2
								)

#Increase n3 to 3000 iterations for final model

alg2 <- sienaAlgorithmCreate(projname = 'cooperation_status', 
								n3 = 3000, 
								seed = 12345, 
								doubleAveraging = 0, 
								diagonalize = 0.2
								)

#Specify the effects for kinship (base) model 

eff.1 <- getEffects(cooperation.data)
eff.1 <- includeEffects(eff.1, gwespFF, inPopSqrt, outActSqrt)
eff.1 <- includeEffects(eff.1, X, interaction1 = 'kinship')

#Run kinship (base) model

base.model <- siena07(alg1, 
						data = cooperation.data, 
						effects = eff.1,
						initC = T,
            useCluster = T,
            nbrNodes = 4,
            batch = F,
            verbose = F
						)

summary(base.model)

#Check for time heterogeneity

tt.1 <- sienaTimeTest(base.model)
summary(tt.1)

#Seems like it may be due to changing density and transitivity,
#so lets include dummies

eff.1.2 <- includeTimeDummy(eff.1, density)
eff.1.2 <- includeTimeDummy(eff.1.2, gwespFF)



base.model.2 <- siena07(alg1, 
						data = cooperation.data, 
						effects = eff.1.2,
						initC = T,
                        useCluster = T,
                        nbrNodes = 4,
                        batch = F,
                        verbose = F
                        )
summary(base.model.2)

tt.2 <- sienaTimeTest(base.model.2)
summary(tt.2)

#Now specify the full model

eff.2 <- includeEffects(eff.1.2, egoX, altX, simX, 
									interaction1 = 'status')
eff.2 <- includeEffects(eff.2, egoX, altX, simX, 
									interaction1 = 'strength')
eff.2 <- includeEffects(eff.2, egoX, altX, simX, 
									interaction1 = 'income')
eff.2 <- includeEffects(eff.2, egoX, altX, simX, 
									interaction1 = 'log.age')
eff.2 <- includeEffects(eff.2,totAlt, name = 'status', 
									interaction1 = 'cooperation')
eff.2 <- includeEffects(eff.2, effFrom, name = 'status', 
									interaction1 = 'strength')
eff.2 <- includeEffects(eff.2, effFrom, name = 'status', 
									interaction1 = 'income')
eff.2 <- includeEffects(eff.2, effFrom, name = 'status', 
									interaction1 = 'log.age')

full.model <- siena07(alg1, 
						data = cooperation.data,
						effects = eff.2,
						initC = T,
                        useCluster = T,
                        nbrNodes = 4,
                        batch = F,
                        verbose = F
                        )
summary(full.model)

tt.3 <- sienaTimeTest(full.model)
summary(tt.3)

eff.3 <- includeTimeDummy(eff.2, simX, interaction1 = "log.age")

full.model.2 <- siena07(alg1, 
                      data = cooperation.data,
                      effects = eff.3,
                      initC = T,
                      useCluster = T,
                      nbrNodes = 4,
                      batch = F,
                      verbose = F
)

summary(full.model)

tt.4 <- sienaTimeTest(full.model.2)
summary(tt.4)

eff.3.2 <- includeTimeDummy(eff.2, outActSqrt)

full.model.3.2 <- siena07(alg1, 
data = cooperation.data, 
effects = eff.3.2, 
initC = T,
useCluster = T,
nbrNodes = 4,
batch = F,
verbose = F,
prevAns = full.model
)

summary(full.model.3.2)

tt.4 <- sienaTimeTest(full.model.3)
summary(tt.4)

eff.3.3 <- includeTimeDummy(eff.3.2, simX, interaction1 = "log.age")

full.model.3.3 <- siena07(alg1, 
                          data = cooperation.data, 
                          effects = eff.3.3, 
                          initC = T,
                          useCluster = T,
                          nbrNodes = 4,
                          batch = F,
                          verbose = F,
                          prevAns = full.model.3
)
summary(full.model.3.3)

tt.5 <- sienaTimeTest(full.model.3.3)
summary(tt.5)

#Re-run the full model with larger n3

full.model.3.4 <- siena07(alg2, 
						data = cooperation.data, 
						effects = eff.3.3, 
						initC = T,
                        useCluster = T,
                        nbrNodes = 4,
                        batch = F,
                        verbose = F,
						prevAns = full.model.3.3
						)
summary(full.model.3.4)




#Get p-values & CI for all parameters

results <- full.model.3.4
parameter <- results$effects$effectName
estimate <- results$theta
st.error <- sqrt(diag(results$covtheta))
norm.var <- estimate/st.error
pval <- 2*pnorm(abs(norm.var), lower.tail = F)


coop.status.out <- data.frame(parameter,
                       estimate = round(estimate, 3),
                       st.error = round(st.error, 3),
                       norm.var = round(norm.var, 2),
                       pval = round(pval, 3),
                       OR = round(exp(estimate), 2),
                        CI.low = round(exp(estimate -
                       	 	qnorm(0.975)*st.error), 2),
                        CI.high = round(exp(estimate 
                        	+ qnorm(0.975)*st.error), 2)
                        )
                       
coop.status.out

write.csv(coop.status.out, file = "cooperation.status.outputs.csv")

#Check partial model with only log age included 

eff.3 <- includeEffects(eff.1.2, egoX, altX, simX,
										 interaction1 = 'log.age')
eff.3 <- includeEffects(eff.3, egoX, altX, simX,
										 interaction1 = 'status')
eff.3 <- includeEffects(eff.3, totAlt, name = 'status', 
										interaction1 = 'cooperation')
eff.3 <- includeEffects(eff.3, effFrom, name = 'status', 
										interaction1 = 'log.age')
eff.3 <- includeEffects(eff.3, effFrom, name = 'status', 
										interaction1 = 'strength')
eff.3 <- includeEffects(eff.3, effFrom, name = 'status', 
										interaction1 = 'income')

age.model <- siena07(alg1, 
					data = cooperation.data, 
					effects = eff.3,
					initC = T,
                    useCluster = T,
                    nbrNodes = 4,
                    batch = F,
                    verbose = F
                    )

summary(age.model)

age.model.2 <- siena07(alg2, 
						data = cooperation.data, 
						effects = eff.3, 
						initC = T,
                         useCluster = T,
                         nbrNodes = 4,
                         batch = F,
                         verbose = F,
                         prevAns = age.model
                         )

results <- age.model.2
parameter <- results$effects$effectName
estimate <- results$theta
st.error <- sqrt(diag(results$covtheta))
norm.var <- estimate/st.error
pval <- 2*pnorm(abs(norm.var), lower.tail = F)


partial.out <- data.frame(parameter,
                       estimate = round(estimate, 3),
                       st.error = round(st.error, 3),
                       norm.var = round(norm.var, 2),
                       pval = round(pval, 3),
                       OR = round(exp(estimate), 2),
                        CI.low = round(exp(estimate -
                                             qnorm(0.975)*st.error), 2),
                        CI.high = round(exp(estimate + qnorm(0.975)*st.error), 2)
                        )
                       
partial.out

write.csv(partial.out, file ="./outputs/partial-model-output.csv")


