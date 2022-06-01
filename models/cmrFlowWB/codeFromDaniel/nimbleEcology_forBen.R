

## fastest - using dCJS distribution
## very efficient for binary-valued CJS data
## if you're not doing multi-state CR, then you should use this.

library(nimbleEcology)


load('./models/cmrFlowWB/dataOut/eh_2002200320042005200620072008200920102011201220132014_wb obear.RData')

y <- eh$eh
(nCohorts <- nrow(unique(eh$cohorts)))
(nSeasons <- nrow(unique(eh$seasons)))
seasonArray <- c(3,4,1,2,3,4,1,2,3,4,1,2)

first <- eh$first #apply(y, 1, function(x) min(which(x !=0)))
last <- eh$last
cohort = ((eh$cohorts) - min(eh$cohorts) + 1)$cohort #can't be a data frame or tibble

zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
zInitsNA <- ifelse(is.na(eh$flow), NA, 1)



myConstants <- list(N = nrow(y), 
                    T = ncol(y), 
                    first = first,
                    last = last,
                    cohort = cohort, 
                    nCohorts = nCohorts,
                    season = seasonArray, #eh$seasons$season,
                    flow = eh$flow,
                    ## DT changes:
                    ## this is used by both the dCJS and dDHMM distributions
                    length = last - first + 1
                    )


## DT changes:
myData <- list(yCJS = y,    ## data for CJS distribution
               y = y + 1)   ## data for DHMM distribution


initialValues <- function() list(
                                betaInt = rnorm(1, 0, 1),
                                ## DT change:
                                ## don't give phi and p initial values;
                                ## they're deterministic nodes, so they'll be calculated
                                ## in terms of other variables.  Also, when I made changes to these
                                ## in the code, these (unnecessary) initial values were the wrong sizes
                                ##phi = array(runif((myConstants$T - 1) * myConstants$N, 0, 1),c((myConstants$T - 1), myConstants$N)),
                                ##p =   array(runif((myConstants$T - 1) * myConstants$N, 0, 1),c((myConstants$T - 1), myConstants$N)),
                                z = zInitsNA,
                                betaPhi = array(runif((myConstants$T - 1) * nCohorts, 0, 1),c((myConstants$T - 1), nCohorts)),
                                betaP =   array(runif((myConstants$T - 1) * nCohorts, 0, 1),c((myConstants$T - 1), nCohorts)),
                                betaPhiCohort = array(runif(nCohorts, 0, 1),c(nCohorts)),
                                betaPCohort =   array(runif(nCohorts, 0, 1),c(nCohorts)),
                                betaFlow = array(rnorm(2 * 4 * nCohorts, 0, 1), c(2, 4, nCohorts)),
                                betaFlowCohort = array(rnorm(2 * nCohorts, 0, 1), c(2, nCohorts)),
                                betaFlowTop = rnorm(2, 0, 1)
                            )


## if you change this FALSE to TRUE
## this makes the dataset smaller - only 200 observations,
## for quicker testing
if(FALSE) {
    newN <- 200
    oldN <- dim(y)[1]
    set.seed(0)
    indToKeep <- sample(1:oldN, size = newN, replace = FALSE)
}

## this removes the very last observation,
## since first[2376] = T = 12, which is not allowed
## to have the first observation occur on the final sampling period
## for either CJS or DHMM distributions
if(TRUE) {
    indToKeep <- which(first < 12)
}

myConstants <- list(
    N = newN,
    T = myConstants$T,
    first = myConstants$first[indToKeep],
    last = myConstants$last[indToKeep],
    cohort = myConstants$cohort[indToKeep],
    nCohorts = myConstants$nCohorts,
    season = myConstants$season,
    flow = myConstants$flow[indToKeep,],
    length = myConstants$length[indToKeep]
)

myData <- list(
    yCJS = myData$yCJS[indToKeep,],
    y = myData$y[indToKeep,]
)

zInitsNA <- zInitsNA[indToKeep,]









## code using CJS distribution
hmm.phiT_pT_cohort_flowCohortHierCJS <- nimbleCode({
    ## DT changes:
    ##delta[1] <- 1                    # Pr(alive t = 1) = 1
    ##delta[2] <- 0                    # Pr(dead t = 1) = 0
    ##
    for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
            logit(phi[t,i]) <- 
                betaInt +
                betaPhi[t,cohort[i]] + 
                betaFlow[1,season[t],cohort[i]] * flow[i,t] +
                betaFlow[2,season[t],cohort[i]] * flow[i,t] * flow[i,t]
        }
        ## DT changes:
        ## time t = first[i]:
        ## note this first value of p[] is not acually used by the dCJS distribution,
        ## but we include it for correctness
        p[first[i],i] <- 1
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
            ## DT changes:
            ## note the indexing on betaP:
            logit(p[t,i]) <- betaP[t-1,cohort[i]]             # prior detection
        }
    }
    ##    
    betaInt ~ dnorm(0,1)
    betaFlowTop[1] ~ dnorm(0,1)
    betaFlowTop[2] ~ dnorm(0,1)
    ##    
    for (c in 1:nCohorts){
        # mean values
        betaPhiCohort[c] ~ dnorm(0,1)
        betaPCohort[c] ~ dnorm(0,1)
        betaFlowCohort[1,c] ~ dnorm(betaFlowTop[1],1)
        betaFlowCohort[2,c] ~ dnorm(betaFlowTop[2],1)
        for (t in 1:(T-1)){ 
            betaPhi[t,c] ~ dnorm(betaPhiCohort[c],1)
            betaP[t,c] ~ dnorm(betaPCohort[c],1)
        }
    }
    ##    
    # back-transform for examining output
    for (c in 1:nCohorts){
        betaPhiCohortOut[c] <- 1/(1 + exp(-betaPhiCohort[c]))
        betaPCohortOut[c] <- 1/(1 + exp(-betaPCohort[c]))
        for (t in 1:(T-1)){ 
            betaPhiOut[t,c] <- 1/(1 + exp(-betaPhi[t,c]))
            betaPOut[t,c] <- 1/(1 + exp(-betaP[t,c])) 
        }
    }
    ##    
    for (s in 1:nSeasons){
        for (c in 1:nCohorts){
            betaFlow[1,s,c] ~ dnorm(betaFlowCohort[1,c],1)
            betaFlow[2,s,c] ~ dnorm(betaFlowCohort[2,c],1)
        }   
    }
    ##    
    # likelihood
    for (i in 1:N){
        ## DT changes:
        ##z[i,first[i]] ~ dcat(delta[1:2])
        ##for (j in (first[i]+1):(last[i])){
        ##    z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1, i])
        ##    y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1, i])
        ##}
        yCJS[i,first[i]:last[i]] ~ dCJS_vv(probSurvive = phi[first[i]:last[i], i],
                                           probCapture = p[first[i]:last[i], i],
                                           len = length[i])
    }
})






set.seed(0)

system.time(
    Rmodel <- nimbleModel(
        code = hmm.phiT_pT_cohort_flowCohortHierCJS,
        constants = myConstants,
        data = myData,              
        inits = initialValues(),
        calculate = FALSE
    )
)




Rmodel$calculate()
## latent state (200 ind): -2847.746
## dCJS_vv (200 individuals): -1199.098
## dCJS_vv (all but last observation): -1217.127


parametersToSave <- c("betaInt", 
                      "betaPhi", "betaP", "betaPhiCohort", "betaPCohort",
                      "betaPhiOut", "betaPOut", "betaPhiCohortOut", "betaPCohortOut", 
                      "betaFlow",
                      "betaFlowCohort", "betaFlowTop")  




system.time(
    conf <- configureMCMC(
        Rmodel,
        monitors = parametersToSave
    )
)













## something like 3x slower,
## using DHMM distribution,
## but this approach would allow for
## multi-state CR models.



## model code using DHMMo distribution
hmm.phiT_pT_cohort_flowCohortHierDHMM <- nimbleCode({
    delta[1] <- 1                    # Pr(alive t = 1) = 1
    delta[2] <- 0                    # Pr(dead t = 1) = 0
    ##
    for (i in 1:N){
        for (t in 1:(T-1)){ # loop over time
            logit(phi[t,i]) <- 
                betaInt +
                betaPhi[t,cohort[i]] + 
                betaFlow[1,season[t],cohort[i]] * flow[i,t] +
                betaFlow[2,season[t],cohort[i]] * flow[i,t] * flow[i,t]
            # prior survival
            ##
            gamma[1,1,t,i] <- phi[t,i]         # Pr(alive t -> alive t+1)
            gamma[1,2,t,i] <- 1 - phi[t,i]     # Pr(alive t -> dead t+1)
            gamma[2,1,t,i] <- 0              # Pr(dead t -> alive t+1)
            gamma[2,2,t,i] <- 1              # Pr(dead t -> dead t+1)
            ##            
            ## DT changes:
            ## definition of omega is moved below, to make it
            ## correctly condition on the first (positive) observation
            ##logit(p[t,i]) <- betaP[t,cohort[i]]             # prior detection
            ##omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
            ##omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
            ##omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
            ##omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
        ## DT changes:
        ## need to pad the gamma matrix with an extra t=T row, to ensure it's
        ## always a matrix.  This values are never actually used (except maybe for internal checking of row sums = 1),
        ## but defining them is necessary.
        gamma[1,1,T,i] <- 0
        gamma[1,2,T,i] <- 1
        gamma[2,1,T,i] <- 0
        gamma[2,2,T,i] <- 1
        ## DT changes:
        ## time period t = first[i]: guaranteed detection:
        omega[1,1,first[i],i] <- 0       # Pr(alive t -> non-detected t)
        omega[1,2,first[i],i] <- 1           # Pr(alive t -> detected t)
        omega[2,1,first[i],i] <- 1              # Pr(dead t -> non-detected t)
        omega[2,2,first[i],i] <- 0              # Pr(dead t -> detected t)
        ## DT changes:
        ## time t > first[i]:
        for(t in (first[i]+1):last[i]) {
            logit(p[t,i]) <- betaP[t-1,cohort[i]]             # prior detection
            omega[1,1,t,i] <- 1 - p[t,i]       # Pr(alive t -> non-detected t)
            omega[1,2,t,i] <- p[t,i]           # Pr(alive t -> detected t)
            omega[2,1,t,i] <- 1              # Pr(dead t -> non-detected t)
            omega[2,2,t,i] <- 0              # Pr(dead t -> detected t)
        }
    }
    ##    
    betaInt ~ dnorm(0,1)
    betaFlowTop[1] ~ dnorm(0,1)
    betaFlowTop[2] ~ dnorm(0,1)
    ##    
    for (c in 1:nCohorts){
        # mean values
        betaPhiCohort[c] ~ dnorm(0,1)
        betaPCohort[c] ~ dnorm(0,1)
        betaFlowCohort[1,c] ~ dnorm(betaFlowTop[1],1)
        betaFlowCohort[2,c] ~ dnorm(betaFlowTop[2],1)
        for (t in 1:(T-1)){ 
            betaPhi[t,c] ~ dnorm(betaPhiCohort[c],1)
            betaP[t,c] ~ dnorm(betaPCohort[c],1)
        }
    }
    ##    
    # back-transform for examining output
    for (c in 1:nCohorts){
        betaPhiCohortOut[c] <- 1/(1 + exp(-betaPhiCohort[c]))
        betaPCohortOut[c] <- 1/(1 + exp(-betaPCohort[c]))
        for (t in 1:(T-1)){ 
            betaPhiOut[t,c] <- 1/(1 + exp(-betaPhi[t,c]))
            betaPOut[t,c] <- 1/(1 + exp(-betaP[t,c])) 
        }
    }
    ##    
    for (s in 1:nSeasons){
        for (c in 1:nCohorts){
            betaFlow[1,s,c] ~ dnorm(betaFlowCohort[1,c],1)
            betaFlow[2,s,c] ~ dnorm(betaFlowCohort[2,c],1)
        }   
    }
    ##    
    # likelihood
    for (i in 1:N){
        ## DT changes:
        ##z[i,first[i]] ~ dcat(delta[1:2])
        ##for (j in (first[i]+1):(last[i])){
        ##    z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1, i])
        ##    y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1, i])
        ##}
        y[i,first[i]:last[i]] ~ dDHMMo(init = delta[1:2],
                                       probTrans = gamma[1:2, 1:2, first[i]:last[i], i],
                                       probObs = omega[1:2, 1:2, first[i]:last[i], i],
                                       len = length[i],
                                       checkRowSums = 1)
    }
})









set.seed(0)

## you'll get warnings that the data 'yCJS' is not used, and the 'z' initial
## values are not in the model.  Those don't cause any problems,
## and let us use the same myData and initialValue() for both models.
system.time(
    Rmodel <- nimbleModel(
        code = hmm.phiT_pT_cohort_flowCohortHierDHMM,
        constants = myConstants,
        data = myData,              
        inits = initialValues(),
        calculate = FALSE
    )
)





Rmodel$calculate()
## latent state: -29141.62
## latent state (200 ind): -2847.746
## dDHMMo (200 ind): -1199.098 (same as CJS)
## dDHMMo (all but last observation): -1217.127 (same as CJS)


parametersToSave <- c("betaInt", 
                      "betaPhi", "betaP", "betaPhiCohort", "betaPCohort",
                      "betaPhiOut", "betaPOut", "betaPhiCohortOut", "betaPCohortOut", 
                      "betaFlow",
                      "betaFlowCohort", "betaFlowTop")  




system.time(
    conf <- configureMCMC(
        Rmodel,
        monitors = parametersToSave
    )
)
