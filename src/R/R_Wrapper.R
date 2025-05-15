rm(list=ls())

# Note: Install these packages using install.packages("pkgname") if not already available
library(Biostrings)
library(vcfR)
library(expm)
library(reticulate)

# Set the Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define Directory Paths
input.data.dir <- "../../additional_files/"
out.data.dir <- "../../out/"
vcf.data.dir <- "../../vcf_files/"
pedigree_dir <- "../../pedigree_files/"

# Load R Scripts
source("makePHYLO.R")
source("makeVCFpedigreeTEMPv2.R")
source("bootSOMA.R")

#Set Up Python Virtual Environment
virtualenv_create("mutSOMA_env")
use_virtualenv("mutSOMA_env", required = TRUE)
py_install(c("numpy", "scipy", "matplotlib", "pandas"))

#Load Python mutSOMA Implementation
source_python("../Python/mutSOMA_new.py")

#Compute Nucleotide Base Frequencies from Reference Genome
poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
freqBASE<-colSums(freqBASE[1:19,])
probBASE<-freqBASE/sum(freqBASE)

#Define Input VCF Files
tree13.vec<-file.path(vcf.data.dir, c("tree13_relaxed.vcf",
                                      "tree13_moderate.vcf",
                                      "tree13_strict1.vcf",
                                      "tree13_strict2.vcf",
                                      "tree13_strict3.vcf",
                                      "tree13_strict4.vcf",
                                      "tree13_strict5.vcf",
                                      "tree13_strict6.vcf"))


tree14.vec<-file.path(vcf.data.dir, c("tree14_relaxed.vcf",
                                      "tree14_moderate.vcf",
                                      "tree14_strict1.vcf",
                                      "tree14_strict2.vcf",
                                      "tree14_strict3.vcf",
                                      "tree14_strict4.vcf",
                                      "tree14_strict5.vcf",
                                      "tree14_strict6.vcf"))

#Define Corresponding Effective Genome Sizes 
gs13 <- c(
  389369818,   #relaxed
  387491779,   #moderate
  382632366,   #strict1
  387491779,   #strict2
  323689759,   #strict3
  323189259,   #strict4
  259329982,   #strict5
  259329982    #strict6
)

gs14 <- c(
  389848942,  # relaxed
  388684187,  # moderate
  384857231,  # strict1
  388684187,  # strict2
  332172596,  # strict3
  332672096,  #strict4
  283905454,  #strict5
  283905454   #strict6
)

gamma<-NULL
gamma_flagged<-NULL
LSE<-NULL
LSE_flagged<-NULL

#Run mutSOMA Analysis for Each VCF Pair
for (a in 1:length(tree13.vec)){
  
  tree13.name<-tree13.vec[a]
  tree14.name<-tree14.vec[a]
  gs.in<-mean(gs13[a],gs14[a])

  pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in,
                                  input.dir = input.data.dir,
                                  tree13 = tree13.name,
                                  tree14 = tree14.name)
  pedigree<-pedigree[[1]]
  
  tree13.basename <- basename(tree13.name)
  tree14.basename <- basename(tree14.name)
  
  ped_path <- file.path(pedigree_dir, paste0("pedigree_", tree13.basename, "_", tree14.basename, ".txt"))
  
  if (!dir.exists(out.data.dir)) {
    dir.create(out.data.dir, recursive = TRUE)
  }
  
  if (!dir.exists(pedigree_dir)) dir.create(pedigree_dir, recursive = TRUE)
  
  write.table(pedigree,file = ped_path, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)
  
  out_name <- paste("EST", tree13.basename, tree14.basename, sep="_")
  
  out<-mutSoma(pedigree_path = ped_path,
               base_probs = c(probBASE[1], probBASE[2], probBASE[3], probBASE[4]),
               prop_het = 0.1,
               num_starts = 10,
               out_dir = out.data.dir,
               out_name = out_name)
  
              gamma[a] <- out$estimates[[1]]$gamma
              LSE[a] <- out$estimates[[1]]$LSE
              
              if (length(out$estimates_flagged) > 0) {
                gamma_flagged[a] <- out$estimates_flagged[[1]]$gamma
                LSE_flagged[a] <- out$estimates_flagged[[1]]$LSE
              } else {
                gamma_flagged[a] <- NA
                LSE_flagged[a] <- NA
              }
              
    pedigree <- as.data.frame(out$input)
    colnames(pedigree) <- c("time0", "time1", "time2", "D.value")
    
    plot_df <- as.data.frame(do.call(rbind, out$plot))
    predicted_divergence <- plot_df[, 4]
    observed_divergence <- pedigree$D.value
    
    pedigree$div.obs <- observed_divergence
    pedigree$div.pred <- predicted_divergence
    pedigree$residual <- observed_divergence - predicted_divergence
    
    out$pedigree <- pedigree
    
    # Run Bootstrap for File 
    bout <- bootSOMA(
      pedigree.data = out,
      Nboot = 10,
      out.dir = out.data.dir,
      out.name = paste("BOOT", tree13.name, tree14.name, sep = "_")
    )
    
 }

#Save Estimates in a Summary File 
outcollect<-cbind(gamma, LSE, gamma_flagged, LSE_flagged, gs13, gs14, tree13.vec, tree14.vec)
write.table(outcollect,
            file = file.path(out.data.dir, paste0("complete_EST_summary_table.txt")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


