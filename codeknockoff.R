library(doParallel)
library(knockoff)
library(matlib)
library(pracma)

set.seed(1234)

# Paramètres du problème
n = 400           # nombres d'observations
p = 120           # nombres de variables
k = 60            # nombre de variables d'hypothèses non-nulles
amplitude = 3.5   # amplitude du signal 
q = 0.2           # niveau de controle du FDR


# création de la matrice de design X
unscaled_X = matrix(rnorm(n*p),n,p)
X = scale(unscaled_X)              #Renormalisation

# création de la réponse Y
nonzero = sample(p, k)                           # Choix des variables non nulles
beta = amplitude * (1:p %in% nonzero) / sqrt(n)  # Construction de beta (+ renormalisation ??)
y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)                                  # Création de Y

invsig = Inverse(t(X) %*% X)                     # Calcul de l'inverse de la matrice de Gram

zvalues = invsig %*% t(X) %*% y                  # Calcul des  Z-valeurs
for ( i in 1 : p) {
  zvalues[i] = zvalues[i]/sqrt(invsig[i,i])}

precision = 100 # taille de la liste des thresholds testés
tlist = linspace(0,max(zvalues),precision) #liste des potentiels thresholds

P = function(t) {2*(1-pnorm(t))} # Proba que  |N(0,1)|>t

nbselectBH = function(list,z) 
{ {for (i in 1:length(list) ) (list[i] <- sum(abs(z)>=list[i]))} 
  list }     # cardinal des coordonnées sélectionnées par BHq (|Z_j |>=t)pour chaque t dans tlist

threshold = function(z,list) { min(which((((p*P(list))/nbselectBH(list,z))<=q))) }#le threshold
# Plus petit t permettant un controle d'uune estimation du FDR
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
