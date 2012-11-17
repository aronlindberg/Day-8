# Loading libraries and data
library(TraMineR)
library(cluster)
library(WeightedCluster)
data(biofam)

# Setting up the dataset
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
                     labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"), right=FALSE)

bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child", "Left/Child", "Left/Married/Child", "Divorced")

bf.shortlab <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")

biofam.seq <-seqdef(biofam[;10:25], states=bf.shortlb, labels=bf.states, weights=biofam$wp00tbgs)

## distance matrix based on state properties
properties <- matrix(c(# left, married, child, divorced
  0, 0, 0, 0,  # parent
  1, 0, 0, 0,  # left
  0, 1, .5, 0, # marr
  1, 1, 0, 0,  # left+marr
  0, 0, 1, 0,  # child
  1, 0, 1, 0,  # left+child
  1, 1, 1, 0,  # left+marr+child
  .5, 1, .5, 1 # divorced
), 8, 4, byrow=TRUE)
sm <- as.matrix(dist(properties))
indel <- .5*max(sm)
dOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=sm, full.matrix=FALSE)


## 2. Hierchical clustering using weights: Ward and WPGMA
weight <- attr(biofam.seq, "weight")
cl.ward <- hclust(dOM, method="ward", members=weight)
cl.mcquitty <- hclust(dOM, method="mcquitty", members=weight)

## 3. Dendrogram side-by-side
par(mfrow=c(1,2))
plot(cl.ward, labels=FALSE, main="Ward")
plot(cl.mcquitty, labels=FALSE, main="WPGMA")
dev.off()

par(mfrow=c(1,2))
plot(cl.ward, labels=FALSE, main="Ward")
plot(cl.mcquitty, labels=FALSE, main="WPGMA")
dev.off()

## 4. 3-cluster solution
cl.hw3 <- cutree(cl.ward, k=3)

## Dtermine labels by looking at I-plots
seqIplot(biofam.seq, group=cl.hw3, sort="from.end")

## Cluster 1: Parental, Cluster 2: Left Home - No Child, 3: Staying with parents

## 5. Silhouette
plot(silhouette(cl.hw3, dmatrix=as.matrix(dOM)))
## Last group is clearly defined.
## The second one looks badly defined with negative silhouette
##   values for a lot of sequences.

## 6. Comparing outcome accounting for weights with
##     ASW without accounting for weights.
wcClusterQuality(dOM, cl.hw3, weight = weight)
summary(silhouette(cl.hw3, dmatrix=as.matrix))


wcClusterQuality(dOM, cl.hw3, weight = weight)
summary(silhouette(cl.hw3, dmatrix=as.matrix(dOM)))


## 7. Logistic regression

cl1 <- as.numeric(cl.hw3==1)
logreg1 <- glm(cl1 ~ sex + birthyr + plingu02, family=binomial, data=biofam)
summary(logreg1)

cl2 <- as.numeric(cl.hw3==2)
logreg2 <- glm(cl2 ~ sex + birthyr + plingu02, family=binomial, data=biofam)
summary(logreg2)

cl3 <- as.numeric(cl.hw3==3)
logreg3 <- glm(cl3 ~ sex + birthyr + plingu02, family=binomial, data=biofam)
summary(logreg3)


## 8. Quality measure for a range of PAM solutions
(pamRange <- wcKMedRange(dOM, 2:10, weight = weight))
plot(pamRange,stat=c("PBC", "HG", "ASW"))
## k=6 seems to be a good compromise solution


## 9. Label the cluster of the 6-group solution
set.seed(13425)
cluster.pam6 <- wcKMedoids(dOM, k = 6, weight = weight)
cl.pam6 <- factor(cluster.pam6$clustering)
seqIplot(biofam.seq, group=cl.pam6, sort="from.end")
pam.labels <- c("Staying with parents","Late Parenthood","Married without child","Solo","Staying married with parents","Early Parenthood")
cl.pam6.fact <- factor(cl.pam6, labels=pam.labels)
seqmtplot(biofam.seq, group=cl.pam6.fact)