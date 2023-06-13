library(knockoff)
library(doParallel)
library(matlib)
library(pracma)

set.seed(1234)

# Paramètres du problème
n = 3000          # nombres d'observations
p = 1000           # nombres de variables
k = 30            # nombre de variables d'hypothèse non-nulles
amplitude = 3.5   # amplitude du signal 
q = 0.2           # niveau de controle du FDR
precision = 100 # taille de la liste des thresholds testés

P = function(t) {2*(1-pnorm(t))} # Proba que  |N(0,1)|>t

nbselectBH = function(list,z) 
{ {for (i in 1:length(list) ) (list[i] <- sum(abs(z)>=list[i]))} 
  list }     # cardinal des coordonnées sélectionnées par BHq (|Z_j |>=t)pour chaque t dans tlist

threshold = function(z,list) { min(which((((p*P(list))/nbselectBH(list,z))<=q))) }#le threshold
# Plus petit t permettant un controle d'uune estimation du FDR
selectBHq = function(z,list) which(abs(z)>= tlist[threshold(z,list)]) #les hypothèse sélectionnées

fdp = function(selected) {sum(beta[selected] == 0) / max(1, length(selected))} #calcul du FDP

power = function(selected) {sum(beta[selected] != 0) / max(1,k)} # calcul de la puissance
#Les N simulations
N = 10
listfdpknockoffplus = c()
listfdpknockoff = c()
listfdpBHq = c()
listpowerknockoffplus = c()
listpowerknockoff = c()
listpowerBHq = c()
startime = Sys.time()
for (i in 1:N) {
  # création de la matrice de design X
  unscaled_X = matrix(rnorm(n*p),n,p)
  
  #Renormalisation
  X = scale(unscaled_X)
  
  # création de la réponse Y
  nonzero = sample(p, k) # Choix des variables non nulles
  beta = amplitude * (1:p %in% nonzero) / sqrt(n) # Construction de beta
  y.sample = function(X) X %*% beta + rnorm(n)
  y = y.sample(X) # Création de Y
  
  invsig = solve(t(X) %*% X) # Calcul de l'inverse de la matrice de Gram
  
  zvalues = invsig %*% t(X) %*% y
  for ( i in 1 : p) {zvalues[i] = zvalues[i]/sqrt(invsig[i,i])}# Calcul des  Z-valeurs
  
  result1 = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=q)
  result1bis = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=q,offset=0)
  listfdpknockoffplus = c(listfdpknockoffplus,fdp(result1$selected))
  listfdpknockoff = c(listfdpknockoff,fdp(result1bis$selected))
  tlist = linspace(0,max(zvalues),precision) #liste des potentiels thresholds
  result2 = selectBHq(zvalues,tlist)
  listfdpBHq = c(listfdpBHq,fdp(result2))
  
  listpowerknockoffplus = c(listpowerknockoffplus,power(result1$selected))
  listpowerknockoff = c(listpowerknockoff,power(result1bis$selected))
  listpowerBHq = c(listpowerBHq, power(result2))}

print(mean(listfdpknockoffplus))
print(mean(listfdpknockoff))
print(mean(listfdpBHq))
print(mean(listpowerknockoffplus))
print(mean(listpowerknockoff))
print(mean(listpowerBHq))
print(endtime-startime)