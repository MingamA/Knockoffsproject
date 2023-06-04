library(doParallel)
library(knockoff)
library(matlib)
library(pracma)

set.seed(1234)

# Problem parameters
n = 400          # number of observations
p = 120           # number of variables
k = 60            # number of variables with nonzero coefficients
amplitude = 3.5   # signal amplitude (for noise level = 1)
q = 0.2        # niveau de controle du FDR


# Generate the variables from a multivariate normal distribution
unscaled_X = matrix(rnorm(n*p),n,p)
X = scale(unscaled_X)

# Generate the response from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero) / sqrt(n)
y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)

invsig = Inverse(t(X) %*% X)

zvalues = invsig %*% t(X) %*% y
for ( i in 1 : p) {
  zvalues[i] = zvalues[i]/sqrt(invsig[i,i])}

precision = 100 # taille de la liste de thresholds testés
tlist = linspace(0,max(zvalues),precision) #liste des potentiels thresholds

P = function(t) {2*(1-pnorm(t))} # Proba que  |N(0,1)|>t

nbselectBH = function(list,z) 
{ {for (i in 1:length(list) ) (list[i] <- sum(abs(z)>=list[i]))} 
  list }     # cardinal des coordonnées sélectionnées par BHq pour les t dans tlist

threshold = function(z,list) { min(which((((p*P(list))/nbselectBH(list,z))<=q))) }#le threshold

selectBHq = function(z,list) which(abs(z)>= tlist[threshold(z,list)]) #les hypothèse sélectionnées

result1 = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=q)
print(result1)  

result2 = selectBHq(zvalues,tlist)
print(result2)

fdp = function(selected) {sum(beta[selected] == 0) / max(1, length(selected))} #calcul du FDP

fdpknockoff = fdp(result1$selected)
print(fdpknockoff)

fdpBHq = fdp(result2)
print(fdpBHq)



