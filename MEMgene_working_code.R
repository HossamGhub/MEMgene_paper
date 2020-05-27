#########################################################################
# Spatial outlier locus detection with Forester et al. 2015 simulations #
#########################################################################
# Set working directory!
# Code assumes that the following files and folder are in the working directory:
# - Folder 'Simulations' from Supplementary Materials of Forester et al. 2016
# - Folder 'Surfaces_Sample', dito
# - File 'sample500.csv'
# - File 'subsample100.csv'

# Load R packages that are needed:
# --------------------------------
require(maptools)
require(SDMTools)
require(spdep)
require(adespatial)

#####################
# Prepare data import
#####################

filenames.Gen <- paste("./Simulations/",
                       rep("Gsims/G", 240),
                       rep(c("S01", "S05", "S10", "S50"), each=60),
                       rep(rep(c("D03", "D05", "D10", "D15", "D25","D50"), each=10),4),
                       paste("/R", 0:9, sep=""), ".csv", sep="")

filenames.Gen <- c(filenames.Gen, 
                   paste("./Simulations/",
                         rep("H1sims/H1", 240),
                         rep(c("S01", "S05", "S10", "S50"), each=60),
                         rep(rep(c("D03", "D05", "D10", "D15", "D25","D50"), each=10),4),
                         paste("/R", 0:9, sep=""), ".csv", sep=""))
filenames.Gen <- c(filenames.Gen, 
                   paste("./Simulations/",
                         rep("H5sims/H5", 240),
                         rep(c("S01", "S05", "S10", "S50"), each=60),
                         rep(rep(c("D03", "D05", "D10", "D15", "D25","D50"), each=10),4),
                         paste("/R", 0:9, sep=""), ".csv", sep=""))
filenames.Gen <- c(filenames.Gen,
                   paste("./Simulations/",
                         rep("H9sims/H9", 240),
                         rep(c("S01", "S05", "S10", "S50"), each=60),
                         rep(rep(c("D03", "D05", "D10", "D15", "D25","D50"), each=10),4),
                         paste("/R", 0:9, sep=""), ".csv", sep=""))

filenames.Env <- paste("L10", rep(c("H1", "H5", "H9"), each=10), 
                       paste("R", 1:10, sep=""), "_aa.asc", sep="")

filenames.Env <- paste("./Surfaces_Sample/", filenames.Env, sep="")

Selection <- list()
Selection$G <- readAsciiGrid(
  paste("./Surfaces_Sample/gradient1024x1024_N0_S1.txt", sep=""))
Selection$H1 <- list()
for(i in 1:10) Selection$H1[[i]] <- readAsciiGrid(filenames.Env[i])
Selection$H5 <- list()
for(i in 11:20) Selection$H5[[i-10]] <- readAsciiGrid(filenames.Env[i])
Selection$H9 <- list()
for(i in 21:30) Selection$H9[[i-20]] <- readAsciiGrid(filenames.Env[i])


######################
# Spatial coordinates 
######################

# Coordinates 
# -----------
Sample <- unlist(read.csv("sample500.csv", header=FALSE))
coord <- data.matrix(read.csv(filenames.Gen[1])[Sample, c(2,3)])

# Subsample of n=100:
Sample.Rand <- unlist(read.csv("subsample100.csv", header=FALSE))
coord.Rand <- data.matrix(read.csv(filenames.Gen[1]))[Sample.Rand, c(2,3)]

# Extract values of env variable at sampling locations
# ----------------------------------------------------
Env <- list()
Env$G <- extract.data(coordinates(coord), Selection$G)
Env$G <- data.frame(matrix(Env$G, length(Sample), 10))
Env$H1 <- as.data.frame(Reduce(cbind,lapply(Selection$H1, function(ls) 
  extract.data(coordinates(coord), ls))))
Env$H5 <- as.data.frame(Reduce(cbind,lapply(Selection$H5, function(ls) 
  extract.data(coordinates(coord), ls))))
Env$H9 <- as.data.frame(Reduce(cbind,lapply(Selection$H9, function(ls) 
  extract.data(coordinates(coord), ls))))
names(Env$H1) <-names(Env$H5) <-names(Env$H9) <- names(Env$G) <- paste("R", 0:9, sep="")

Env.Rand <- list()
Env.Rand$G <- extract.data(coordinates(coord.Rand), Selection$G)
Env.Rand$G <- data.frame(matrix(Env.Rand$G, length(Sample.Rand), 10))
Env.Rand$H1 <- as.data.frame(Reduce(cbind,lapply(Selection$H1, function(ls) 
  extract.data(coordinates(coord.Rand), ls))))
Env.Rand$H5 <- as.data.frame(Reduce(cbind,lapply(Selection$H5, function(ls) 
  extract.data(coordinates(coord.Rand), ls))))
Env.Rand$H9 <- as.data.frame(Reduce(cbind,lapply(Selection$H9, function(ls) 
  extract.data(coordinates(coord.Rand), ls))))
names(Env.Rand$H1) <-names(Env.Rand$H5) <-names(Env.Rand$H9) <- names(Env.Rand$G) <- 
  paste("R", 0:9, sep="")

remove(Selection)  # Remove large data set from memory (optional)

#########################
# MEM: vectors and values
#########################
# Currently using Gabriel graph with inverse distance weights

# MEM for full sample (n=500):

nb<-graph2nb(gabrielneigh(coord), sym=TRUE)      # Gabriel graph: neighbor definition
listW <- nb2listw(nb,style="W")                  # Spatial weights matrix
# Use inverse distance weights:
disttri<-nbdists(nb,coord)
fdist<-lapply(disttri,function(x) x^(-1))        # Use inverse distance weights
listW <-nb2listw(nb,glist=fdist,style="W")       # Revised spatial weights matrix
tmp <- scores.listw(listW, MEM.autocor = "all")  # Eigen analysis
mem <- list(vectors = as.matrix(tmp), values = attr(tmp, "values"))
mem$values <- mem$values / abs(sum(mem$values))  # Rescale eigenvalues to Moran's I

# MEM for subsample (n=100)

nb<-graph2nb(gabrielneigh(coord.Rand), sym=TRUE)      # Gabriel graph: neighbor definition
listW <- nb2listw(nb,style="W")                  # Spatial weights matrix
# Use inverse distance weights:
disttri<-nbdists(nb,coord.Rand)
fdist<-lapply(disttri,function(x) x^(-1))        # Use inverse distance weights
listW <-nb2listw(nb,glist=fdist,style="W")       # Revised spatial weights matrix
tmp <- scores.listw(listW, MEM.autocor = "all")  # Eigen analysis
mem.Rand <- list(vectors = as.matrix(tmp), values = attr(tmp, "values"))
mem.Rand$values <- mem.Rand$values / abs(sum(mem.Rand$values))






##############################
# Begin of simulations
##############################

# Global settings:
cutoffs <- abs(qnorm(c(0.025, 0.005, 0.0005)))   # Cutoffs for step 1
nPerm <- 199                                     # Number of replicates for MSR test

# Prepare storage of results:
File <- TPR <- FPR <- MC <- TPR.msr <- FPR.msr <- TPR.both <- FPR.both <- 
  Trend.max.neutral <- 
  array(NA, dim=c(10, 6, 4, 4, 2, length(cutoffs)),
        dimnames=list(paste("R", 0:9, sep=""),
                      c("D03", "D05", "D10", "D15", "D25","D50"),
                      c("S01", "S05", "S10", "S50"),
                      c("G", "H1", "H5", "H9"),
                      c("Full", "Rand"),
                      c("0.05", "0.01", "0.001")))

# Input files:
File[]  <- filenames.Gen

# Loop through simulated data
# ---------------------------

for(i in 1:4)      # Habitat spatial distribution
{
  cat("\n")
  cat("i=", i)
  for(s in 1:4)    # Selection strength
  {
    cat("\n")
    cat(" s =", s, ": ")
    for(d in 1:6)  # Dispersal limitation
    {
      cat(d)
      for(r in 1:10)  # Replicate landscapes
      {
        # Read in genetic data:
        # ---------------------
        Loci <- data.matrix(read.csv(File[r,d,s,i,1,1])
                            [Sample, seq(from=8, by=2, length.out=100)])
        Loci.Rand <- data.matrix(read.csv(File[r,d,s,i,1,1])
                                 [Sample.Rand, seq(from=8, by=2, length.out=100)])
        
        # Maximum trend in neutral loci: (for Supplementary Material)
        # ------------------------------
        Trend.max.neutral[r,d,s,i,1,1] <- 
          max(unlist(lapply(summary(lm(Loci[,-1] ~ coord)), 
                            function(ls) ls$r.squared)))
        Trend.max.neutral[r,d,s,i,2,1] <- 
          max(unlist(lapply(summary(lm(Loci.Rand[,-1] ~ coord.Rand)), 
                            function(ls) ls$r.squared)))
        
        # Correlations with MEM axes:
        # ---------------------------
        R.YV <- cor(Loci, mem$vectors, use="pairwise.complete.obs")
        R.YV.Rand <- cor(Loci.Rand, mem.Rand$vectors, use="pairwise.complete.obs")
        
        S <- apply(R.YV^2, 2, mean)          # S=Average power spectrum
        S.Rand <- apply(R.YV.Rand^2, 2, mean)
        
        MC[r,d,s,i,1,1] <- S %*% mem$values             # Moran's I of S (n=500)
        MC[r,d,s,i,2,1] <- S.Rand %*% mem.Rand$values   # Moran's I of S (n=100)
        
        
        # Calculate z-scores for power spectra:
        # -------------------------------------
        # For n=500:
        Dev <- sweep(R.YV^2, 2, S, "/") - 1    # Variance ratio minus 1
        Dev[Dev > 0] <- 0                      # Set positive deviations to 0
        Dev <- apply(Dev, 1, sum)              # Sum of negative deviations
        a <- scale(Dev)                        # Standardize (z-scores)
        
        # For n=100:
        Dev.Rand <- sweep(R.YV.Rand^2, 2, S.Rand, "/") - 1 
        Dev.Rand[Dev.Rand > 0] <- 0        
        Dev.Rand <- apply(Dev.Rand, 1, sum)
        a.Rand <- scale(Dev.Rand)                   
        
        
        # Step 1: identify outlier loci (3 cutoffs):
        # ------------------------------------------
        
        Candidates <- Candidates.Rand <- list()
        
        for(h in 1:length(cutoffs))
        {
          Candidates[[h]] <- c(1:100)[abs(a)>cutoffs[h]]                 # Outlier loci 
          TPR[r,d,s,i,1,h]  <-as.numeric(is.element(1, Candidates[[h]])) # True pos. rate
          FPR[r,d,s,i,1,h]  <-
            (length(Candidates[[h]]) - TPR[r,d,s,i,1,h])/99              # False pos. rate
          
          Candidates.Rand[[h]] <- c(1:100)[abs(a.Rand)>cutoffs[h]]       # For n=100
          TPR[r,d,s,i,2,h] <-
            as.numeric(is.element(1, Candidates.Rand[[h]]))     
          FPR[r,d,s,i,2,h]  <-
            (length(Candidates.Rand[[h]]) - TPR[r,d,s,i,2,h])/99 
        }
        
        
        # Preparation for Step 2: Spectral randomization test 
        # ---------------------------------------------------
        
        # MEM for Env and coordinates (as spurious predictors):
        #------------------------------------------------------
        R.XV.Env <- cor(Env[[i]][,r], mem$vectors)
        R.XV.xcoord <- cor(coord[,1], mem$vectors)
        R.XV.ycoord <- cor(coord[,2], mem$vectors)
        
        R.XV.Rand.Env <- cor(Env.Rand[[i]][,r], mem.Rand$vectors)
        R.XV.Rand.xcoord <- cor(coord.Rand[,1], mem.Rand$vectors)
        R.XV.Rand.ycoord <- cor(coord.Rand[,2], mem.Rand$vectors)
        
        # function to perform MSR test:
        #------------------------------
        get.pvalue.msr <- function(r.XV=R.XV, r.YV=R.YV, nPerm=199)
        {
          R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
          R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
          Cor.obs <- abs(as.vector(r.YV %*% t(r.XV)))
          Cor.rand <- abs(r.YV %*% t(R.XV.rand))
          P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
          P.values.MSR
        }
        
        # MSR test for ALL loci (i.e., Step 2 only, without Step 1):
        # ----------------------------------------------------------
        b <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV, nPerm=nPerm)
        b.Rand <- get.pvalue.msr(r.XV=R.XV.Rand.Env, r.YV=R.YV.Rand, nPerm=nPerm)
        TPR.msr[r,d,s,i,1,1]  <- as.numeric(b[1] < 0.05)       # True positive rate
        FPR.msr[r,d,s,i,1,1]  <- sum(b[-1] < 0.05)/99          # False positive rate
        TPR.msr[r,d,s,i,2,1]  <- as.numeric(b.Rand[1] < 0.05)  # True positive rate
        FPR.msr[r,d,s,i,2,1]  <- sum(b.Rand[-1] < 0.05)/99     # False positive rate
        
        # MSR test for outlier loci (i.e., Step 2 after Step 1):
        # ------------------------------------------------------
        for(h in 1:length(cutoffs))
        {
          b2 <- b2.Rand <- c()
          if(length(Candidates[[h]]) > 0)
          {
            b2 <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV[Candidates[[h]],], nPerm=nPerm)
            if(is.element(1, Candidates[[h]]))
            {
              TPR.both[r,d,s,i,1,h]  <- as.numeric(b2[1] < 0.05)      # True positive rate
              FPR.both[r,d,s,i,1,h]  <- sum(b2[-1] < 0.05)/99         # False positive rate
            } else {
              TPR.both[r,d,s,i,1,h]  <- 0                             # True positive rate
              FPR.both[r,d,s,i,1,h]  <- sum(b2 < 0.05)/99         # False positive rate
            }
            
          } else 
          {
            TPR.both[r,d,s,i,1,h] <- FPR.both[r,d,s,i,1,h] <- 0
          }
          if(length(Candidates.Rand[[h]]) > 0)
          {
            b2.Rand <- get.pvalue.msr(r.XV=R.XV.Rand.Env, 
                                      r.YV=R.YV.Rand[Candidates.Rand[[h]],], nPerm=nPerm)
            if(is.element(1, Candidates.Rand[[h]]))
            {
              TPR.both[r,d,s,i,2,h]  <- as.numeric(b2.Rand[1] < 0.05)     
              FPR.both[r,d,s,i,2,h]  <- sum(b2.Rand[-1] < 0.05)/99        
            } else {
              TPR.both[r,d,s,i,2,h]  <- 0     
              FPR.both[r,d,s,i,2,h]  <- sum(b2.Rand < 0.05)/99        
            }
          } else 
          {
            TPR.both[r,d,s,i,2,h] <- FPR.both[r,d,s,i,2,h] <- 0
          }
        }
      }
    }
  }
}

####################
# END OF SIMULATIONS
####################

