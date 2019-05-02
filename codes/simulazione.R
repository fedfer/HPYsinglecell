rm(list = ls())


# number of runs
Runs=3


### independent proposals
### different alphas and ds for the j populations


########### set parameters ##########

# total number of species among all populations
N=3000

# number of arms
J=8
#J=15

# number of species present in the J populations
NN=c(2500,2500,2500,2500,2500,2500,2500,2500)
#NN=c(2500,2500,2500,2500,2500,2500,2500,2500,2500,2500,2500,2500,2500,2500,2500)

# set the parameters for Zipf laws in the J populations
Zipfpar=c(1.3,1.3,2,2,2,2,2,2)
#Zipfpar=c(1.2,1.2,1.4,1.4,1.1,1.5,1.5,1.6,1.6,1.4,1.7,1.7,1.8,1.1,1.8)

# Initial sample sizes
n_init=c(30,30,30,30,30,30,30,30)
#n_init=c(100,100,100,100,100,100,100,100,100,100,100,100,100,100,100)

# additional sample
addsample=300

# iterations in MCMC for number of tables and hyperp HPY given the initial sample
iter=100
burnin=5

# iterations in MCMC for hyper HPY after every observations
iter1=4
burnin1=2

# Good Turing tuning parameter
C=(1+sqrt(2))*sqrt(3)

# function to generate random samples from a (finite) discrete distribution 
gendiscr <- function(supp,masses) {
k=length(supp)
Pcum=c(0, cumsum(masses))
U=runif(1,0,1)
for (i in 1:k){
        if(Pcum[i]<U & U<=Pcum[i+1]){
        answer=supp[i];
}
}
return(answer)
}

# function for log of the generalized factorial coefficients
loggfc <- function(alpha,n) {
if (n<=0) {
return(vector(mode="numeric",length=0))
}
logalpha <- log(alpha)
c <- vector(mode="numeric",length=1)
c[1] <- logalpha
if (n==1) {
return(c)
}
for (m in seq(1,len=n-1)) {
x <- vector(mode="numeric",length=m+1)
a <- log(m-alpha*seq(2,len=m-1))+c[seq(2,len=m-1)]
b <- logalpha+c[seq(1,len=m-1)]
themax <- pmax(a,b)
x[1] <- log(m-alpha*1)+c[1]
x[seq(2,len=m-1)] <- themax + log(exp(a-themax) + exp(b-themax))
x[m+1] <- logalpha + c[m]
c <- x
}
return(c)
}




DATAfinal=list()
DATAfinal[[paste('HPY', 1)]] <- matrix(rep(0,addsample*Runs),nrow=Runs) 
DATAfinal[[paste('uniform', 2)]] <- matrix(rep(0,addsample*Runs),nrow=Runs) 
DATAfinal[[paste('Oracle', 3)]] <- matrix(rep(0,addsample*Runs),nrow=Runs) 
DATAfinal[[paste('Good-Turing', 4)]] <- matrix(rep(0,addsample*Runs),nrow=Runs) 



# weights
WEIGTHS=list()
WEIGTHS[[paste('weights.HPY', 1)]] <- matrix(rep(0,J*Runs),nrow=Runs) 
WEIGTHS[[paste('weight.uniform', 2)]] <- matrix(rep(0,J*Runs),nrow=Runs) 
WEIGTHS[[paste('weight.Oracle', 3)]] <- matrix(rep(0,J*Runs),nrow=Runs) 
WEIGTHS[[paste('weight.Good-Turing', 4)]] <- matrix(rep(0,J*Runs),nrow=Runs) 



# progress bar
pbRUN <- winProgressBar(title="RUN", label="0% done", min=0, max=100, initial=0)

for (III in 1:Runs){



########## simulate true distributions in the J populations ########

# labels of the species
labels=1:N

# frequencies in the populations
freq <- list() 
for (i in 1:J) freq[[paste('frequecies', i)]] <- rep(0,NN[i]) 

for (i in 1:J){
freq1=1:NN[i]
freq2=sum(freq1^(-Zipfpar[i]))
for (j in 1:NN[i]){
freq[[i]][j]=((1/j)^Zipfpar[i])/freq2
}
}

# labebls in the populations
pop <- list() 
for (i in 1:J) pop[[paste('population', i)]] <- sample(labels, NN[i], replace = FALSE, prob = NULL) 
## with this choice, the sharing of species in different populations will be very low. Alternatevely, we can use
# for (i in 1:J) pop[[paste('population', i)]] <- 1:NN[i] 

#### pop[[i]] are the labels of arm i, with frequencies freq[[i]]




####### generate initial data for every population #########

data = list()
for (i in 1:J) data[[paste('data population', i)]] <- rep(0,n_init[i]) 
for (j in 1:J){
for (i in 1:n_init[j]){
data[[j]][i]=gendiscr(pop[[j]],freq[[j]])
}
}

# distinct values in the initial sample (among all J pop)
Kini = unique(unlist(data))
# total number of distincts in the initial sample
tot.dist=length(Kini)




####### MCMC - inference on the number of tables and hyperparameters HPY #######

# space to storage results (hyperpar and tables mjk)
alpha=matrix(rep(0,iter*J),ncol=J)
d=matrix(rep(0,iter*J),ncol=J)
gamma=rep(0,iter)
nu=rep(0,iter)

mjk_est <- list() 
for (i in 1:J) mjk_est[[paste('population', i)]] <- matrix(rep(0,tot.dist*iter), ncol=tot.dist)

# set initial values for the hyperpar of the HPY
alpha[1,]=1
d[1,]=0.5
gamma[1]=1
nu[1]=0.5

# find an initial value for the table allocations (running a chinese franchise conditional to obs)
tables = list()
for (i in 1:J) {tables[[paste('tables_assignemnts', i)]] <- rep(0,(n_init[i]))}

for (j in 1:J) {
# first customer sits at table 1
tables[[j]][1]=1

for (i in 2:(n_init[j])) {
# if the i-th observation is new, then the customer is sitting at a new table
if( data[[j]][i] %in% data[[j]][1:(i-1)] == FALSE){
tables[[j]][i]=max(tables[[j]][1:(i-1)])+1 
}
else {
# if the i-th obs is old, then the customer can sit either at a new table or one of the other tables serving 
# the same dish (having the same observed value). We find these tables and their numbers of customers njt. .
equalcust=which(data[[j]][1:(i-1)] %in% data[[j]][i])
possibletables=unique(tables[[j]][equalcust])
lengpos=length(possibletables)
njt.=rep(0,lengpos)
for (u in 1:lengpos) {
njt.[u]=length(which(tables[[j]][1:(i-1)] %in% possibletables[u]))
}

freqtables=c(njt.-d[1],alpha[1,j]+(length(possibletables))*d[1,j])
freqtables=freqtables/sum(freqtables)

possibletables=c(possibletables,max(tables[[j]][1:(i-1)])+1)

tables[[j]][i]=gendiscr(possibletables,freqtables)
}
}
}
# having resampled all tji, compute all mjk
for (w in 1:tot.dist){
for (u in 1:J){
equalcust=which(data[[u]] %in% Kini[w])
mjk_est[[u]][1,w]=length(unique(tables[[u]][equalcust]))
}
}



### MCMC table allocations tji and hyperpar

# progress bar
pbtab <- winProgressBar(title="Tables progress bar", label="0% done", min=0, max=100, initial=0)

for (i in 2:iter){

## update table allocations
for (j in 1:J){
for (n in 1:n_init[j]){

# update only if there are other obs with the same value (otherwise it means the customer is 
# sitting alone at his table and there are no other tables serving the same dish in that restourant)
if (data[[j]][n] %in% data[[j]][-n] == TRUE) {

# find the possible tables where tjn can sit and their frequencies njt.
equalcust=which(data[[j]] %in% data[[j]][n])
possibletables=unique(tables[[j]][equalcust])

# compute njt.
lengpos=length(possibletables)
njt.=rep(0,lengpos)
for (u in 1:lengpos) {
njt.[u]=length(which(tables[[j]][-n] %in% possibletables[u]))
}

# compute m.. , total number of tables
m.. =0
for (u in 1:J){
m.. = m.. + length(unique(tables[[u]]))
}

# compute m.k , total number of tables eating the same dish of n-th customer in restourant j
mjk=rep(0,J)
for (u in 1:J){
equalcust=which(data[[u]] %in% data[[j]][n])
samedish=tables[[u]][equalcust]
mjk[u]=length(unique(samedish))
}
m.k=sum(mjk)

# if the individual we are resampling is sitting alone at his tables (but there are other tables 
# serving that dish) remove the table and subtract one from tables counts
if (tables[[j]][n] %in% tables[[j]][-n] == FALSE){
njt.=njt.[!possibletables %in% tables[[j]][n]]
possibletables=possibletables[!possibletables %in% tables[[j]][n]]
m.k=m.k-1
m..=m..-1
}

# compute mj. , total number of tables in rest j
mj. = length(possibletables)

# compute the probabilities of sitting at each of the possible tables
possibletables=c(possibletables, max(tables[[j]]+1))
proballoc=c(njt. - d[i-1,j], (alpha[i-1,j] + mj.*d[i-1,j])*((m.k - nu[i-1])/(gamma[i-1] + m..)))
proballoc= proballoc/(sum(proballoc))

# generate the new allocations value tjn
tables[[j]][n]=gendiscr(possibletables,proballoc)   
}
}
}

# store the mjk for all restourants and distinct values
for (w in 1:tot.dist){
for (u in 1:J){
equalcust=which(data[[u]] %in% Kini[w])
mjk_est[[u]][i,w]=length(unique(tables[[u]][equalcust]))
}
}


## update hyperparameters using metropolis hastings with independent proposals

for (w in 1:J) {

# update d
dprop=runif(1,0,1)
# compute acceptance propabilities
accept=lgamma((alpha[i-1,w]/dprop)+sum(mjk_est[[w]][i,]))+loggfc(dprop,n_init[w])[sum(mjk_est[[w]][i,])]-lgamma((alpha[i-1,w]/d[i-1,w])+sum(mjk_est[[w]][i,]))-loggfc(d[i-1,w],n_init[w])[sum(mjk_est[[w]][i,])]+lgamma(alpha[i-1,w]/d[i-1,w])-lgamma(alpha[i-1,w]/dprop)+log((1-dprop)/(1-d[i-1,w]))
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) d[i,w]=dprop  else d[i,w]=d[i-1,w]

# update alpha 
alphaprop=rexp(1,1)
# compute acceptance propabilities
accept=lgamma((alphaprop/d[i,w])+sum(mjk_est[[w]][i,]))-lgamma(alphaprop+n_init[w])-lgamma((alpha[i-1,w]/d[i,w])+sum(mjk_est[[w]][i,]))+lgamma(alpha[i-1,w]+n_init[w])+lgamma(alpha[i-1,w]/d[i,w])-lgamma(alphaprop/d[i,w])+lgamma(alphaprop)-lgamma(alpha[i-1,w])
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) alpha[i,w]=alphaprop else alpha[i,w]=alpha[i-1,w]

}

# compute m.. , total number of tables
m.. = rep(0,J)
for (q in 1:J){
m..[q]=sum(mjk_est[[q]][i,])
}
m..=sum(m..)

# update nu
nuprop=runif(1,0,1)
# compute acceptance propabilities
accept=lgamma((gamma[i-1]/nuprop)+tot.dist)+loggfc(nuprop,m..)[tot.dist]-lgamma(gamma[i-1]/nuprop)+log(1-nuprop)-lgamma((gamma[i-1]/nu[i-1])+tot.dist)-loggfc(nu[i-1],m..)[tot.dist]+lgamma(gamma[i-1]/nu[i-1])+log(1-nu[i-1])
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) nu[i]=nuprop  else nu[i]=nu[i-1]

# update gamma
gammaprop=rexp(1,1)
# compute acceptance propabilities
accept=lgamma((gammaprop/nu[i])+tot.dist)+lgamma(gammaprop)-lgamma(gammaprop/nu[i])-lgamma(gammaprop+m..)-lgamma((gamma[i-1]/nu[i])+tot.dist)-lgamma(gamma[i-1])+lgamma(gamma[i-1]/nu[i])+lgamma(gamma[i-1]+m..)
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) gamma[i]=gammaprop else gamma[i]=gamma[i-1]


# update progress bar
info <- sprintf("%d%% done", round((i/iter)*100))
	setWinProgressBar(pbtab, i/(iter)*100, label=info)
}
close(pbtab)



## results MCMC

# estimate of mjk
mjk=matrix(rep(0,J*tot.dist),nrow=J)
for (j in 1:J) {
for (i in 1:tot.dist){
mjk[j,i]=mean(mjk_est[[j]][burnin:iter,i])
}
}
mjk=round(mjk)

# estimate of m.k
m.k=rep(0,tot.dist)
for (i in 1:tot.dist){
m.k[i]=sum(mjk[,i])
}

# estimates of mj.
mj. = rep(0,J)
for (i in 1:J){
mj.[i]=sum(mjk[i,])
}

# estimate of m..
m..=sum(mj.)


# keep the results just to check acceptances/rejections of proposals
alpha.ini=alpha
d.ini=d
gamma.ini=gamma
nu.ini=nu

# estimates of the hyperparameters
alpha=rep(0,J)
d=rep(0,J)
for (j in 1:J){
alpha[j]=mean(alpha.ini[burnin:iter,j])
d[j]=mean(d.ini[burnin:iter,j])
}
gamma=mean(gamma.ini[burnin:iter])
nu=mean(nu.ini[burnin:iter])




######### algorithms - unif, HPY, Oracle, GT - observe new data ########

# uniform sampling
supp=1:J
masses=rep(1/J,J)
KuniUni=Kini

# HPY
nn=n_init
KuniHPY=Kini
bigK=length(Kini)
# compute n.k
nj.k=matrix(rep(0,tot.dist*J),nrow=J)
for (j in 1:J){
for (w in 1:tot.dist){
nj.k[j,w]=length(which(data[[j]] %in% KuniHPY[w]))
}
}

# oracle
KuniOracle=Kini
missingmass=rep(0,J)
for (j in 1:J) {
labelsseen=which(pop[[j]] %in% KuniOracle)
missingmass[j]=1-sum(freq[[j]][labelsseen])
}

# Good Turing
KuniGT=Kini
nGT=n_init
# find the number of obs with frequency 1 in the joint sample
unGT_tot=as.numeric(names(which(table(unlist(data))==1)))
# find the number of obs with frequency 1 in each population
unGT=list()
for (i in 1:J) {unGT[[paste('frequencies 1 in group', i)]] <- 0}
for (j in 1:J) {
uniGTp=as.numeric(names(which(table(data[[j]])==1)))
unGT[[j]]=unGT_tot[which(unGT_tot %in% uniGTp)]
}

# indicators of discoveries
newobsindHPY=rep(0,addsample)
newobsindOra=rep(0,addsample)
newobsindUni=rep(0,addsample)
newobsindGT=rep(0,addsample)

# indicators for weights
w.HPY=rep(0,J)
w.Ora=rep(0,J)
w.Uni=rep(0,J)
w.GT=rep(0,J)

# space to storage results for estimation of hyperp (to check acc/rej)
d.stor <- list() 
for (i in 1:J) {d.stor[[paste('d.stor', i)]] <- rep(0,addsample*iter1)}
alpha.stor <- list()
for (i in 1:J) {alpha.stor[[paste('alpha.stor', i)]] <- rep(0,addsample*iter1)}
nu.stor=rep(0,addsample*iter1)
gamma.stor=rep(0,addsample*iter1)


# progress bar
pb <- winProgressBar(title="Algorithm progress bar", label="0% done", min=0, max=100, initial=0)


### start the algorithms

for (i in 1:addsample) {

# choose arm - unif
unif=gendiscr(supp,masses)

w.Uni[unif]=w.Uni[unif]+1


# choose arm - HPY
betadraws=rep(0,J)
betazero=rbeta(1, gamma+nu*bigK, m..-bigK*nu)
for (j in 1:J) {
betadraws[j]=rbeta(1, (betazero*(alpha[j]+(mj.[j])*d[j])),((1-betazero)*(alpha[j]+(mj.[j])*d[j])+ nn[j]- mj.[j]*d[j]))
}
armchosen=which(betadraws==max(betadraws))

w.HPY[armchosen]=w.HPY[armchosen]+1


# choose arm - oracle
armchosenOrac=which(missingmass==max(missingmass))
if (length(armchosenOrac)>1){
Oracprob=rep((1/length(armchosenOrac)),length(armchosenOrac))
armchosenOrac=gendiscr(armchosenOrac,Oracprob)
}

w.Ora[armchosenOrac]=w.Ora[armchosenOrac]+1


# choose arm - GoodTuring
R_GT=rep(0,J)
for (j in 1:J) {
R_GT[j]=(length(unGT[[j]])/nGT[j])+C*sqrt(log(4*(sum(nGT))/nGT[j]))
}
armchosenGT=which(R_GT==max(R_GT))
if (length(armchosenGT)>1){
GTprob=rep((1/length(armchosenGT)),length(armchosenGT))
armchosenGT=gendiscr(armchosenGT,GTprob)
}

w.GT[armchosenGT]=w.GT[armchosenGT]+1


# sample a new value for each arm
newobservations=rep(0,J)
for (j in 1:J) {
newobservations[j]=gendiscr(pop[[j]],freq[[j]])
}

# unif
newunif=newobservations[unif]
# hpy
newobsHPY=newobservations[armchosen]
# oracle
newobsOrac=newobservations[armchosenOrac]
# GT
newobsGT=newobservations[armchosenGT]


# update uniform
if(newunif %in% KuniUni == FALSE){
newobsindUni[i]=1
KuniUni=c(KuniUni,newunif)
}

# update HPY tables and distincts
if(newobsHPY %in% KuniHPY == FALSE){
bigK=bigK+1
mj.[armchosen]=mj.[armchosen]+1
m..=m..+1
KuniHPY=c(KuniHPY,newobsHPY)
newobsindHPY[i]=1
m.k=c(m.k,1)

zeror=matrix(rep(0,J),ncol=1)
nj.k=cbind(nj.k,zeror)
mjk=cbind(mjk,zeror)

nj.k[armchosen,bigK]=1
mjk[armchosen,bigK]=1
}
else {

olddistinct=which(newobsHPY %in% KuniHPY)
# compute the probability of sitting at a new table conditionally that the observation is old
probnewold=c(nj.k[armchosen,olddistinct]-d[armchosen]*mjk[armchosen,olddistinct],(alpha[armchosen]+mj.[armchosen]*d[armchosen])*((m.k[olddistinct]-nu)/(gamma+m..)))
probnewold=probnewold/sum(probnewold)
bern=rbinom(1,1,probnewold[1])

if(bern == 0){
# if the old observation forms a new table, update tables counts
mj.[armchosen]=mj.[armchosen]+1
m.k[olddistinct]=m.k[olddistinct]+1
mjk[armchosen,olddistinct]=mjk[armchosen,olddistinct]+1
m..=m..+1
}

nj.k[armchosen,olddistinct]=nj.k[armchosen,olddistinct]+1

}

nn[armchosen]=(nn[armchosen])+1


## update hyperparameters using metropolis hastings with independent proposals

# set initial values
d1=matrix(rep(0,iter1*J),ncol=J)
alpha1=matrix(rep(0,iter1*J),ncol=J)
nu1=rep(0,iter1)
gamma1=rep(0,iter1)

# set initial values
d1[1,]=d
alpha1[1,]=alpha
nu1[1]=nu
gamma1[1]=gamma

for (y in 2:iter1){

for (j in 1:J) {
# update d
dprop=runif(1,0,1)
# compute acceptance propabilities
accept=lgamma((alpha1[y-1,j]/dprop)+mj.[j])+loggfc(dprop,nn[j])[mj.[j]]-lgamma((alpha1[y-1,j]/d1[y-1,j])+mj.[j])-loggfc(d1[y-1,j],nn[j])[mj.[j]]+lgamma(alpha1[y-1,j]/d1[y-1,j])-lgamma(alpha1[y-1,j]/dprop)+log((1-dprop)/(1-d1[y-1,j]))
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) d1[y,j]=dprop  else d1[y,j]=d1[y-1,j]

# update alpha 
alphaprop=rexp(1,1)
# compute acceptance propabilities
accept=lgamma((alphaprop/d1[y,j])+mj.[j])-lgamma(alphaprop+nn[j])-lgamma((alpha1[y-1,j]/d1[y,j])+mj.[j])+lgamma(alpha1[y-1,j]+nn[j])+lgamma(alpha1[y-1,j]/d1[y,j])-lgamma(alphaprop/d1[y,j])+lgamma(alphaprop)-lgamma(alpha1[y-1,j])
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) alpha1[y,j]=alphaprop else alpha1[y,j]=alpha1[y-1,j]

}

# update nu
nuprop=runif(1,0,1)
# compute acceptance propabilities
accept=lgamma((gamma1[y-1]/nuprop)+bigK)+loggfc(nuprop,m..)[bigK]-lgamma(gamma1[y-1]/nuprop)+log(1-nuprop)-lgamma((gamma1[y-1]/nu1[y-1])+bigK)-loggfc(nu1[y-1],m..)[bigK]+lgamma(gamma1[y-1]/nu1[y-1])+log(1-nu1[y-1])
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) nu1[y]=nuprop  else nu1[y]=nu1[y-1]

# update gamma
gammaprop=rexp(1,1)
# compute acceptance propabilities
accept=lgamma((gammaprop/nu1[y])+bigK)+lgamma(gammaprop)-lgamma(gammaprop/nu1[y])-lgamma(gammaprop+m..)-lgamma((gamma1[y-1]/nu1[y])+bigK)-lgamma(gamma1[y-1])+lgamma(gamma1[y-1]/nu1[y])+lgamma(gamma1[y-1]+m..)
accepttt=exp(accept)
acceptprob=min(c(accepttt,1))
bern_accept=rbinom(1,1,acceptprob)
if (bern_accept==1) gamma1[y]=gammaprop else gamma1[y]=gamma1[y-1]

}

# storage all the steps of the chain to check acceptances/rejections
for (j in 1:J){
d.stor[[j]][((i-1)*iter1+1):(i*iter1)]=d1[,j]
alpha.stor[[j]][((i-1)*iter1+1):(i*iter1)]=alpha1[,j]
}
nu.stor[((i-1)*iter1+1):(i*iter1)]=nu1
gamma.stor[((i-1)*iter1+1):(i*iter1)]=gamma1

d=rep(0,J)
alpha=rep(0,J)
for (j in 1:J){
d[j]=mean(d1[burnin1:iter1,j])
alpha[j]=mean(alpha1[burnin1:iter1,j])
}
nu=mean(nu1[burnin1:iter1])
gamma=mean(gamma1[burnin1:iter1])


# update Oracle parameters
if (newobsOrac %in% KuniOracle == FALSE) {
KuniOracle=c(KuniOracle,newobsOrac)
newobsindOra[i]=1
# update the missing masses
for (j in 1:J) {
newlabel=which(pop[[j]] %in% newobsOrac)
# the following 'if' is needed only in case of distinct supports of the true distributions, 
# since some populations could not have newobsOrac in their support, so newlabel=NaN (alternatevely, we can use is.na(newlabel)==FALSE)
if (length(newlabel)==1){
missingmass[j]=missingmass[j]-(freq[[j]][newlabel])
}
}
}

# update Good Turing
if (newobsGT %in% KuniGT == FALSE) {
unGT[[armchosenGT]]=c(unGT[[armchosenGT]],newobsGT)
unGT_tot=c(unGT_tot,newobsGT)
KuniGT=c(KuniGT,newobsGT)
newobsindGT[i]=1
}
if (newobsGT %in% unGT == TRUE) {
for (j in 1:J){
if (newobsGT %in% unGT[[j]] == TRUE) {
unGT[[j]]=unGT[[j]][-which(unGT[[j]]==newobsGT)]
}
}
unGT_tot=unGT_tot[-which(unGT_tot==newobsGT)]
}
nGT[armchosenGT]=nGT[armchosen]+1

# update progress bar
info <- sprintf("%d%% done", round((i/addsample)*100))
	setWinProgressBar(pb, i/(addsample)*100, label=info)

}
close(pb)


##### results

distcumHPY=cumsum(newobsindHPY)
distcumUni=cumsum(newobsindUni)
distcumOra=cumsum(newobsindOra)
distcumGT=cumsum(newobsindGT)


DATAfinal[[1]][III,]=distcumHPY 
DATAfinal[[2]][III,]=distcumUni 
DATAfinal[[3]][III,]=distcumOra 
DATAfinal[[4]][III,]=distcumGT 


WEIGTHS[[1]][III,]=w.HPY 
WEIGTHS[[2]][III,]=w.Uni 
WEIGTHS[[3]][III,]=w.Ora 
WEIGTHS[[4]][III,]=w.GT



# update progress bar
info <- sprintf("%d%% done", round((III/Runs)*100))
	setWinProgressBar(pbRUN, (III/Runs)*100, label=info)
}
close(pbRUN)


