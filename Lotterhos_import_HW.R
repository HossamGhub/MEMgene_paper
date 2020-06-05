############################################################################# 
# Spatial outlier locus detection with Lotterhos and Whitlock (2015)'s Data #
#############################################################################
# Set working directory!
# Code assumes that the folder 'dryad/SimFilesLFMM' from the Supplementary
# Materials of Lotterhos and Whitlock (2015) on Dryad is in the working directory. 
# In addition, the following files from the Dryad repository of the present 
# paper should also be in the working directory:
#  '1351142954_453EnviMat.txt'                      # Environmental factor
#  '1351142970_988EnviMat.txt'
#  '1351142986_950EnviMat.txt' 
#  1351142954_453EnviMatPAIRS_ED.txt                # Paired sampling locations
#  1351142970_988EnviMatPAIRS_ED.txt
#  1351142986_950EnviMatPAIRS_ED.txt
#  1351142954_453EnviMatTRANSECTS_ED_Design.txt     # Transect sampling locations
#  1351142970_988EnviMatTRANSECTS_ED_Design.txt
#  1351142986_950EnviMatTRANSECTS_ED_Design.txt
#  SchemeRandom1.txt                                # Random sampling locations
 
# Load packages 
# -------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# if (!requireNamespace("qvalue", quietly = TRUE))
#    BiocManager::install(c("qvalue"))


 require(spdep)
 require(adespatial)
# require(moments)
# require(qvalue)

######################
# Spatial coordinates: 
######################   
   
   Coords <- list()
   #Coords$Pairs <- list()
   # Coords$Pairs$E453 <- read.table("1351142954_453EnviMatPAIRS_ED.txt")
   # Coords$Pairs$E988 <- read.table("1351142970_988EnviMatPAIRS_ED.txt")
   # Coords$Pairs$E950 <- read.table("1351142986_950EnviMatPAIRS_ED.txt")
   # Coords$Transect <- list()
   # Coords$Transect$E453 <- read.table("1351142954_453EnviMatTRANSECTS_ED_Design.txt")
   # Coords$Transect$E988 <- read.table("1351142970_988EnviMatTRANSECTS_ED_Design.txt")
   # Coords$Transect$E950 <- read.table("1351142986_950EnviMatTRANSECTS_ED_Design.txt")
   Coords$Random <- list()
   Coords$Random$E453 <- Coords$Random$E988 <- Coords$Random$E950 <- 
     read.table("SchemeRandom1.txt") 
   
   # for(k in 1:length(Coords$Transect))
   # {
   #   names(Coords$Transect[[k]])[2:3] <- names(Coords$Pairs[[1]])[2:3]
   #   Coords$Transect[[k]][,7:9][is.na(Coords$Transect[[k]][,7:9])] <- FALSE
   # }
   for(i in 1:length(Coords))
   {
     for(k in 1:length(Coords[[i]]))
     {
       b <- order(Coords[[i]][[k]][,3], Coords[[i]][[k]][,2]) # Sort by y, then x
       Coords[[i]][[k]] <- Coords[[i]][[k]][b,]  # correct order!!
     }
   }

###################
# Extract env data:
###################
   
   # Environment data: find all file names in the folder "SimFilesLFMM" that contain 'env':
   Filenames.lfmm <- list.files(paste0(here::here(), "/dryad/SimFilesLFMM/"), pattern="env") ##
   
   # Env <- list()
   # for(i in 1:length(Filenames.env))
   # {
   #   tmp <- read.table(paste("./SimFilesLFMM/", Filenames.env[[i]], sep=""))
   #   Env[[i]] <- unlist(tmp)
   # }
   # 
######################################################################   
# Read in file names, create table of design information for each file
######################################################################   

   # Genetic data: find all file names in the folder "SimFilesLFMM" that contain 'lfmm':
   Filenames.lfmm <- list.files(paste0(here::here(), "/dryad//SimFilesLFMM/"), pattern="lfmm") ##

   
   test <- Reduce(rbind, strsplit(Filenames.lfmm, split="[_=]"))
   test <- cbind(test, Reduce(rbind, strsplit(test[,ncol(test)], split="[.]")))
   tmp <- strsplit(test[,2], split="[.xs]")
   test2 <- matrix(NA, nrow(test), 3)
   for(i in 1:nrow(test)) test2[i,1:length(tmp[[i]])] <- tmp[[i]]
   test2 <- gsub("[A-Z, a-z]","", test2)
   test <- cbind(test, test2)
   dimnames(test) <- list(NULL, c("Demography", "Design", "ID", "Env", "ID2", 
                                  "V6", "NumPops", "V8", "V9", "NumInd", "Design2",
                                  "NumPops2", "NumTrans", "NumInd2"))
   
   # Design matrix (parameter space): each row is one combination of parameter settings:
   Design <- as.data.frame(test[,c(1,2,4,7,10,13,14)])
   Design$Env <- ordered(Design$Env, levels=c(453,988,950))
   Design$Type <- ordered(substr(Design$Design, 1, 1), levels=c("P", "T", "R"))
   Design$Design <- as.character(Design$Design)
   Design$Design[Design$Design == "T30.T3x10"] <- "T30.3x10s"
   Design$Design[Design$Design == "T30.T6x5s"] <- "T30.6x5s"
   Design$NumPops <- as.numeric(as.character(Design$NumPops))

   head(Design)
   
   # Demography: single refugium (1R), two refugia (2R), isolation by distane (IBD), island model (IM)
   # Design: sampling design (R, P, T; see below) and number of pops sampled
   # Env: which of the three replicate landscapes '453', '950', '988'
   # NumPops: how many populations are sampled
   # NumInd: how many individuals sampled per pop
   # Type: random (R), pairs (P), transects (T)
   
####################################################################
# BEGIN OF SIMULATIONS: DEFINE VARIABLES TO GENERATE
####################################################################
# Note: high number of MSR replicates (nPerm) for calculation of empirical p-values
# For normal applications, nPerm=199 should be sufficient for alpha = 0.05

   
   # Global settings:
   # ----------------
   # cutoffs = abs(qnorm(c(0.025, 0.005, 0.0005)))      # z-score cutoff levels (one-sided)
   #                                                    # two-sided are used (i.e., double)
   # nPerm = 199                                        # Number of MSR permutations   
   
   Sites.90 <- c(1:nrow(Design))[Design$NumPops==90]  # All samples with 90 populations:
   
   # Prepare collection of results
   # -----------------------------
   NumLoci <- Rsq <- MoransI <- MoransI.w <- Fst.total <- Fst.sample <- rep(NA, nrow(Design))

   
   # Loop through simulated datasets:
   # --------------------------------
   j=1
   for(j in 1:length(Sites.90))
   {
     # Select the sites that need to be sampled for this run:
     i=Sites.90[j]
     cat("j:", j, ", i:", i, "\n")
     
     # Sample size: 
     # NumPops the number of populations (sites) sampled
     # NumInd is the number of individuals sampled per population
     Site <- rep(1:Design$NumPops[i],each=as.numeric(as.vector(Design$NumInd)[i]))
     
     # Each file has NumPops x NumInd rows (sampled individuals) and 9996 columns (loci:)
     tmp <- read.table(paste0(here::here(),"/dryad/SimFilesLFMM/", Filenames.lfmm[[i]]))
     dim(tmp)
     NumLoci[i] <- ncol(tmp) 
     
     # Split the genetic data by site (pop) and take the mean for each locus = column
     # Note: each locus is one column
     Allele.freqs <- t(sapply(split(tmp, Site), colMeans)/2)
     
     # Check out a few loci for all individuals sampled from the first site:
     split(tmp, Site)[[1]][,1:20]
     
     # Extract the grid coordinates of the sampled sites
     coord <- data.matrix(Coords[[as.numeric(Design$Type[i])]]
                          [[as.numeric(Design$Env[i])]][,2:3])
   
     # Check the first few rows (x and y grid coordinates of the first few sites)
     head(coord)
     
     # Maximum trend in neutral / selected loci:
     # -----------------------------------------
     # Trend <- unlist(lapply(summary(lm(Allele.freqs ~ coord)), 
     #                        function(ls) ls$r.squared))
     # Trend.max.neutral <- max(Trend[1:9900])
     # Trend.max.selected <- max(Trend[-c(1:9900)])
     
     # MEM
     # -----
     nb<-graph2nb(gabrielneigh(coord), sym=TRUE)      # Gabriel graph: neighbor definition
     listW <- nb2listw(nb,style="W")                  # Spatial weights matrix
     disttri<-nbdists(nb,coord)                       # Determine distances
     fdist<-lapply(disttri,function(x) x^(-1))        # Use inverse distance weights
     listW <-nb2listw(nb,glist=fdist,style="W")       # Revised spatial weights matrix
     tmp <- scores.listw(listW, MEM.autocor = "all")  # Eigen analysis of W 
     mem <- list(vectors = as.matrix(tmp), values = attr(tmp, "values"))
     mem$values <- mem$values / abs(sum(mem$values))  # Rescale eigenvalues to Moran's I
     
     # Plot Gabriel graph:
     plot(nb, coord)
     
     # Object 'mem' has two attributes:
     # mem$vectors: mem spatial eigenvectors (n = 90 rows, n - 1 = 89 columns)
     # mem$values: eigen values rescaled to return Moran's I for each mem spatial eigenvector (n - 1= 89 values)
     
     # Range of rescaled eigenvalues determines maximum range of Moran's I for dataset:
     range(mem$values)
     
     # Mean of rescaled eigenvalues is the expected value of Moran's I if there is zero spatial structure.
     # E[Moran's I] = -1/(n-1)
     n=90
     -1/(n-1)
     mean(mem$values) ##-------------------------------------------------------------------
     
     # Step 1: Outlier detection from power spectrum
     # ---------------------------------------------
     # R.XV.Env <- cor(sapply(split(Env[[i]], Site), unique), mem$vectors)
     # R.YV <- cor(Allele.freqs, mem$vectors, use="pairwise.complete.obs")
     # S <- apply(R.YV^2, 2, mean)            # S=Average power spectrum
     # MC[i] <- S %*% mem$values              # Moran's I of S
     # Dev <- sweep(R.YV^2, 2, S, "/") - 1    # Variance ratios minus one
     # Dev[Dev > 0] <- 0                      # Set positive deviations to zero
     # Dev <- apply(Dev, 1, sum)              # Sum of negative deviations
     # a <- scale(Dev)                        # Standardize
     # 
     # Step 1: identify outlier loci
     # -----------------------------
     # Candidates <- list()
     # for(h in 1:length(cutoffs))
     # {
     #   Candidates[[h]] <- c(1:length(a))[abs(a)>cutoffs[h]]  
     #   TP[[h]][i] <- sum(Candidates[[h]] > 9900)
     #   TPR[[h]][i] <- sum(Candidates[[h]] > 9900)/(length(a)-9900)
     #   FPR[[h]][i]  <-sum(Candidates[[h]] <= 9900)/9900
     # }

     # Empirical p-values (compatibility with Lotterhos and Whitlock 2015):
     # -------------------------------------------------------------------
     # p.values <- rep(NA, 100)
     # for(loc in 9901:NumLoci[i])
     # {
     #   p.values[loc-9900] <- (sum(abs(a[1:9900]) > abs(a[loc]))+1)/9901
     # }
     # P.values[[i]] <- p.values
     # Q.values <- qvalue(p.values, lambda=0)
     # TP.empirical[i] <- sum(Q.values$qvalues < 0.01, na.rm=TRUE)
     # TPR.empirical[i] <- mean(Q.values$qvalues < 0.01, na.rm=TRUE)
     
     # Evaluate skewness and kurtosis for neutral loci:
     # ------------------------------------------------
   #   S.neutral <- apply(R.YV[1:9900,]^2, 2, mean)          # S=Average power spectrum
   #   Dev.neutral <- sweep(R.YV[1:9900,]^2, 2, S.neutral, "/") - 1
   #   Dev.neutral[Dev.neutral > 0] <- 0        # TRYING THIS OUT!
   #   Dev.neutral <- apply(Dev.neutral, 1, sum)
   #   Skewness[i] <- skewness(Dev.neutral)        # Skewness: expected = 0
   #   Kurtosis[i] <- kurtosis(Dev.neutral) - 3    # Excess kurtosis: expected = 0
   #   S.neutral.sd[i] <- sd(Dev.neutral)
   #   
   #   # Step 2 only: MSR test
   #   # ---------------------
   #   get.pvalue.msr <- function(r.XV=R.XV, r.YV=R.YV, nPerm=199)
   #   {
   #     R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
   #     R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
   #     Cor.obs <- abs(as.vector(r.YV %*% t(r.XV)))
   #     Cor.rand <- abs(r.YV %*% t(R.XV.rand))
   #     P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
   #     P.values.MSR
   #   }
   #   
   #   R.XV.Env <- cor(sapply(split(Env[[i]], Site), unique), mem$vectors)
   #   tmp.e <- sapply(split(Env[[i]], Site), unique)
   #   Env.R2[i] <- summary(lm(tmp.e ~ coord))$r.squared                 # Linear trend in env
   #   
   #   b <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV, nPerm=nPerm)        # MSR p-values
   #   TPR.msr[i]  <- mean(b[-c(1:9900)] < 0.05)                         # True positive rate
   #   for(loc in 9901:NumLoci[i])
   #   {
   #     p.values[loc-9900] <- (sum(b[1:9900] < b[loc])+1)/9901
   #   }
   #   Q.values <- qvalue(p.values, lambda=0)
   #   TPR.empirical.msr[i] <- mean(Q.values$qvalues < 0.01, na.rm=TRUE)
   #   FPR.msr[i]  <- sum(b[1:9900] < 0.05)/9900                         # False positive rate
   #   
   #   # Steps 1 an 2: MSR test
   #   # ----------------------     
   #   b2 <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV[Candidates[[1]],], nPerm=nPerm)
   #   b3 <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV[1:9900,], nPerm=nPerm)
   #   TPR.both[i]  <- sum(b2[Candidates[[1]] > 9900] < 0.05)/(length(a)-9900)       
   #   n.neutral = sum(Candidates[[1]] <= 9900)
   #   
   #   p.values <- rep(NA, sum(Candidates[[1]] > 9900))
   #   for(loc in (n.neutral+1):length(Candidates[[1]]))
   #   {
   #     p.values[loc - n.neutral] <- (sum(b3 < b2[loc])+1)/9901
   #   }
   #   Q.values <- qvalue(p.values, lambda=0)
   #   TPR.empirical.both[i] <- sum(Q.values$qvalues < 0.01, na.rm=TRUE)/(length(a)-9900)
   #   FPR.both[i]  <- sum(b2[Candidates[[1]] <= 9900] < 0.05)/9900 
   # 
   #   P.values.both[[j]] <- b2[Candidates[[1]] > 9900]
   #   P.empirical.both[[j]] <- p.values
   # }
   # 
   }
   
####################
# End of simulations
####################
   
####################################   
# Prepare result data set: Design.90
####################################

   # Design.90 <- data.frame(Design, NumLoci=NumLoci,
   #                         TPR.1=TPR[[1]], TPR.2=TPR[[2]], TPR.3=TPR[[3]], 
   #                         TPR.both=TPR.both, TPR.msr=TPR.msr, 
   #                         FPR.1=FPR[[1]], FPR.2=FPR[[2]], FPR.3=FPR[[3]], 
   #                         FPR.both=FPR.both, FPR.msr=FPR.msr,
   #                         Skewness=Skewness, Kurtosis=Kurtosis, Env.R2=Env.R2, 
   #                         MC=MC, S.neutral.sd=S.neutral.sd,
   #                         TPR.empirical.both=TPR.empirical.both,
   #                         TPR.empirical=TPR.empirical,
   #                         TPR.empirical.msr=TPR.empirical.msr)[Sites.90,]
   # Design.90$NumInd <- factor(as.character(Design.90$NumInd), levels=c("20", "6"))
   # Design.90$Type <- factor(Design.90$Type, levels=c("R", "T", "P"))
   # Design.90$Demography <- factor(Design.90$Demography, levels=c("2R", "1R","IBD", "IM"))

   Design.90 <- data.frame(Design, NumLoci=NumLoci)[Sites.90,]
   Design.90$NumInd <- factor(as.character(Design.90$NumInd), levels=c("20", "6"))
   Design.90$Type <- factor(Design.90$Type, levels=c("R", "T", "P"))
   Design.90$Demography <- factor(Design.90$Demography, levels=c("2R", "1R","IBD", "IM"))
   
