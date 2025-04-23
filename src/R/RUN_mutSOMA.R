

######################################### 
#### Clean calls (Illumina) Paul + Sujan
#### COMPLETE TREE
#########################################
rm(list=ls())
library(Biostrings)
library(vcfR)
library(expm)
library(optimx)

input.data.dir <- "/Users/adaakinci/Desktop/"
out.data.dir <- "/Users/adaakinci/Desktop/"
vcf.data.dir <- "/Users/adaakinci/Dropbox/RAW/"
#func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"

## Loading source code
source(paste0("~/Desktop/mutsomaR/makePHYLO.R"))
source(paste0("~/Desktop/mutsomaR/makeVCFpedigreeTEMPv2.R"))
source(paste0("~/Desktop/mutSOMA.R"))
source(paste0("~/Desktop/mutsomaR/bootSOMA.R"))

## Reading in the poplar fasta file and determining the nucleotide frequency
poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
freqBASE<-colSums(freqBASE[1:19,])
probBASE<-freqBASE/sum(freqBASE)

## Reading in the filter info
basesin<-read.table(paste(input.data.dir, "Ptricocarpa_bases_at_each_step.txt", sep=""), header=T)
basesin[,5]<-sub(",", "", basesin[,5])
basesin[,5]<-sub(",", "", basesin[,5])
basesin[,5]<-as.numeric(basesin[,5])
#print(basesin[,5])


      #tree13.name<-"tree13.filtered_somatic_SNPs.v1.vcf"
      #tree14.name<-"tree14.filtered_somatic_SNPs.v1.vcf"
      
      tree13.vec<-file.path(vcf.data.dir, c("tree13.somatic.down_20.rep_1.positions.vcf",
                    "tree13.somatic.down_20.rep_2.positions.vcf",
                    "tree13.somatic.down_20.rep_3.positions.vcf",
                    "tree13.somatic.down_25.rep_1.positions.vcf",
                    "tree13.somatic.down_25.rep_2.positions.vcf",
                    "tree13.somatic.down_25.rep_3.positions.vcf",
                    "tree13.somatic.down_30.rep_1.positions.vcf",
                    "tree13.somatic.down_30.rep_2.positions.vcf",
                    "tree13.somatic.down_30.rep_3.positions.vcf",
                    "tree13.somatic.down_35.rep_1.positions.vcf",
                    "tree13.somatic.down_35.rep_2.positions.vcf",
                    "tree13.somatic.down_35.rep_3.positions.vcf",
                    "tree13.somatic.down_40.rep_1.positions.vcf",
                    "tree13.somatic.down_40.rep_2.positions.vcf",
                    "tree13.somatic.down_40.rep_3.positions.vcf",
                    "tree13.somatic.down_45.rep_1.positions.vcf",
                    "tree13.somatic.down_45.rep_2.positions.vcf",
                    "tree13.somatic.down_45.rep_3.positions.vcf"))
      
      
      tree14.vec<-file.path(vcf.data.dir, c("tree14.somatic.down_20.rep_1.positions.vcf",
                    "tree14.somatic.down_20.rep_2.positions.vcf",
                    "tree14.somatic.down_20.rep_3.positions.vcf",
                    "tree14.somatic.down_25.rep_1.positions.vcf",
                    "tree14.somatic.down_25.rep_2.positions.vcf",
                    "tree14.somatic.down_25.rep_3.positions.vcf",
                    "tree14.somatic.down_30.rep_1.positions.vcf",
                    "tree14.somatic.down_30.rep_2.positions.vcf",
                    "tree14.somatic.down_30.rep_3.positions.vcf",
                    "tree14.somatic.down_35.rep_1.positions.vcf",
                    "tree14.somatic.down_35.rep_2.positions.vcf",
                    "tree14.somatic.down_35.rep_3.positions.vcf",
                    "tree14.somatic.down_40.rep_1.positions.vcf",
                    "tree14.somatic.down_40.rep_2.positions.vcf",
                    "tree14.somatic.down_40.rep_3.positions.vcf",
                    "tree14.somatic.down_45.rep_1.positions.vcf",
                    "tree14.somatic.down_45.rep_2.positions.vcf",
                    "tree14.somatic.down_45.rep_3.positions.vcf"))
                    
                    
      gs13<-c(rep(40388007, 3),
               rep(44900782, 3),
               rep(38229529, 3),
               rep(32546285, 3),
               rep(24989388, 3),
               rep(18264725, 3))
               
      gs14<-c(rep(54998919, 3),
               rep(58223515, 3),
               rep(49555495, 3),
               rep(40561041, 3),
               rep(28972798, 3),
               rep(19351957, 3))

      rate.collect<-NULL
      rate.collect2<-NULL
      lsq.collect<-NULL
      lsq.collect2<-NULL
      
        for (a in 1:length(tree13.vec))
        {
            tree13.name<-tree13.vec[a]
            tree14.name<-tree14.vec[a]
            gs.in<-mean(gs13[a],gs14[a])
              
            pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in, 
                                              input.dir = input.data.dir, 
                                              tree13 = tree13.name,
                                              tree14 = tree14.name)
            pedigree<-pedigree[[1]]
            
            out<-mutSOMA(pedigree.data = pedigree, 
                           p0aa= probBASE[1],
                           p0cc= probBASE[2],
                           p0tt= probBASE[3],
                           p0gg= probBASE[4],
                           Nstarts=10,
                           prop.het = 0.1,
                           out.dir = out.data.dir,
                           out.name = paste("complete_EST", tree13.name, tree14.name, sep="_"))
              
            rate.collect[a]<-out$estimates[1,1]
            lsq.collect[a]<-out$estimates[1,3]
            rate.collect2[a]<-out$estimates.flagged[1,1]
            lsq.collect2[a]<-out$estimates.flagged[1,3]
              
        }
  
  outcollect<-cbind(rate.collect, lsq.collect, rate.collect2, lsq.collect2, gs13, gs14, tree13.vec, tree14.vec)
  write.table(outcollect, paste(out.data.dir, "complete_EST_summary_table.text", sep=""))
  
  bout<-bootSOMA(pedigree.data =out,    
                 Nboot=100, 
                 out.dir=out.data.dir, 
                 out.name = paste("complete_BOOT", tree13.name, tree14.name, sep="_"))
  

# 
#   #rm(list=ls())
#   input.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
#   out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/FIGURES/"
#   func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/epiFUNC/"
#   
#   ### Loading source code
#   source(paste(func.dir, "plotFITS.r", sep=""))
#   source(paste(func.dir, "plotfitsFINDdims.R", sep=""))
#   
#   ### Finding plotting dimensions
#   pedigree.list<-c(
#     "complete_EST_tree13.somatic.down_35.rep_2.positions.vcf_tree14.somatic.down_35.rep_2.positions.vcf.Rdata")
#   
#   plotFITS(pedigree.names = pedigree.list, 
#            input.dir=input.data.dir,
#            output.dir=out.data.dir, 
#            plot.type="both",
#            f.w=1.25,
#            f.h=1.75,
#            f.res=800,
#            x.dim=c(0,650),
#            y.dim=c(0, 15e-07),
#            col.vec.data=c("red"),
#            col.vec.line=c("red"),
#            out.name=paste("PLOT", pedigree.list, sep="_"), 
#            header= "",
#            lsq.line="theory",
#            plot.null=FALSE)
# 




  
  ######################################### 
  #### Clean calls (Illumina) Paul + Sujan
  #### COMPLETE TREE
  #### Bootstrapping
  #########################################
  rm(list=ls())
  library(Biostrings)
  library(vcfR)
  library(expm)
  library(optimx)
  
  in.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
  out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
  func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"
  
  ## Loading source code
  source(paste(func.dir, "makePHYLO.R", sep=""))
  source(paste(func.dir, "makeVCFpedigreeTEMPv2.R", sep=""))
  source(paste(func.dir, "mutSOMA.R", sep=""))
  source(paste(func.dir, "bootSOMA.R", sep=""))
  
  
  out<-dget(paste(in.data.dir, "complete_EST_tree13.somatic.down_20.rep_1.positions.vcf_tree14.somatic.down_20.rep_1.positions.vcf.Rdata", sep=""))
  
  bout<-bootSOMA(pedigree.data =out,    
                 Nboot=500, 
                 out.dir=out.data.dir, 
                 out.name = "complete_BOOT_complete_EST_tree13.somatic.down_20.rep_1.positions.vcf_tree14.somatic.down_20.rep_1.positions.vcf.Rdata")
  
 check<-bout$boot.results
 check<-check[which(check[,1] > 0),]
 quantile(check[,1], probs=c(0.025, 0.975))



 #1.325858e-10
 #1.182789e-11 4.088842e-10



     
######################################### 
#### depthX_hetAt1_step2 
#### COMPLETE TREE
#########################################
     rm(list=ls())
     library(Biostrings)
     library(vcfR)
     library(expm)
     library(optimx)
     
     input.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/RAW/"
     out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
     func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"
     
     ## Loading source code
     source(paste(func.dir, "makePHYLO.R", sep=""))
     source(paste(func.dir, "makeVCFpedigreeTEMPv2.R", sep=""))
     source(paste(func.dir, "mutSOMA.R", sep=""))
     source(paste(func.dir, "bootSOMA.R", sep=""))
     
     ## Reading in the poplar fasta file and determining the nucleotide frequency
     poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
     freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
     freqBASE<-colSums(freqBASE[1:19,])
     probBASE<-freqBASE/sum(freqBASE)
     
     ## Reading in the filter info
     basesin<-read.table(paste(input.data.dir, "Ptricocarpa_bases_at_each_step.txt", sep=""), header=T)
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-as.numeric(basesin[,5])
     
     tree13bin<-basesin[which(basesin[,1] == "tree13" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     tree14bin<-basesin[which(basesin[,1] == "tree14" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     
     depth.vec<-tree13bin[,2]
     
     
    for (a in 1:length(depth.vec))
    {
      
      cat("Progress: ----------------------------------------------", a/length(depth.vec), "\n")
      
        gs.13<-tree13bin[which(tree13bin[,2] == depth.vec[a]), 5]
        gs.14<-tree14bin[which(tree14bin[,2] == depth.vec[a]), 5]
        gs.in<-mean(gs.13, gs.14)
        
        tree13.name<-paste("tree13.depth", depth.vec[a], "_hetAt1_step2", sep="")
        tree14.name<-paste("tree14.depth", depth.vec[a], "_hetAt1_step2", sep="")
        
         pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in, 
                                         input.dir = input.data.dir, 
                                         tree13 = tree13.name,
                                         tree14 = tree14.name)
         pedigree<-pedigree[[1]]
        
         ## Running models
         out<-mutSOMA(pedigree.data = pedigree, 
                      p0aa= probBASE[1],
                      p0cc= probBASE[2],
                      p0tt= probBASE[3],
                      p0gg= probBASE[4],
                      Nstarts=20,
                      prop.het = 0.4,
                      out.dir = out.data.dir,
                      out.name = paste("complete_EST", tree13.name, tree14.name, sep="_"))
         
         bout<-bootSOMA(pedigree.data =out,    
                        Nboot=100, 
                        out.dir=out.data.dir, 
                        out.name = paste("complete_BOOT", tree13.name, tree14.name, sep="_"))
         
    }    
     
     
     
     

     
######################################### 
#### depthX_hetAt1_step3
#### COMPLETE TREE
#########################################
     rm(list=ls())
     library(Biostrings)
     library(vcfR)
     library(expm)
     library(optimx)
     
     input.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/RAW/"
     out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
     func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"
     
     ## Loading source code
     source(paste(func.dir, "makePHYLO.R", sep=""))
     source(paste(func.dir, "makeVCFpedigreeTEMPv2.R", sep=""))
     source(paste(func.dir, "mutSOMA.R", sep=""))
     source(paste(func.dir, "bootSOMA.R", sep=""))
     
     ## Reading in the poplar fasta file and determining the nucleotide frequency
     poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
     freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
     freqBASE<-colSums(freqBASE[1:19,])
     probBASE<-freqBASE/sum(freqBASE)
     
     ## Reading in the filter info
     basesin<-read.table(paste(input.data.dir, "Ptricocarpa_bases_at_each_step.txt", sep=""), header=T)
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-as.numeric(basesin[,5])
     
     tree13bin<-basesin[which(basesin[,1] == "tree13" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     tree14bin<-basesin[which(basesin[,1] == "tree14" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     
     depth.vec<-tree13bin[,2]
     
     
     for (a in 9:length(depth.vec))
     {
       
       cat("Progress: ----------------------------------------------", a/length(depth.vec), "\n")
       
       gs.13<-tree13bin[which(tree13bin[,2] == depth.vec[a]), 5]
       gs.14<-tree14bin[which(tree14bin[,2] == depth.vec[a]), 5]
       gs.in<-mean(gs.13, gs.14)
       
       tree13.name<-paste("tree13.depth", depth.vec[a], "_hetAt1_step3", sep="")
       tree14.name<-paste("tree14.depth", depth.vec[a], "_hetAt1_step3", sep="")
       
       pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in, 
                                       input.dir = input.data.dir, 
                                       tree13 = tree13.name,
                                       tree14 = tree14.name)
       pedigree<-pedigree[[1]]
       
       ## Running models
       out<-mutSOMA(pedigree.data = pedigree, 
                    p0aa= probBASE[1],
                    p0cc= probBASE[2],
                    p0tt= probBASE[3],
                    p0gg= probBASE[4],
                    Nstarts=20,
                    prop.het = 0.4,
                    out.dir = out.data.dir,
                    out.name = paste("complete_EST", tree13.name, tree14.name, sep="_"))
       
       bout<-bootSOMA(pedigree.data =out,    
                      Nboot=100, 
                      out.dir=out.data.dir, 
                      out.name = paste("complete_BOOT", tree13.name, tree14.name, sep="_"))
       
     }    
     
   
     
     
     
     
######################################### 
#### depthX_hetAt1_step2 
#### WITHIN TREE
#########################################
     rm(list=ls())
     library(Biostrings)
     library(vcfR)
     library(expm)
     library(optimx)
     
     input.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/RAW/"
     out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
     func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"
     
     ## Loading source code
     source(paste(func.dir, "makePHYLO.R", sep=""))
     source(paste(func.dir, "makeVCFpedigreeTEMPv2.R", sep=""))
     source(paste(func.dir, "mutSOMA.R", sep=""))
     source(paste(func.dir, "bootSOMA.R", sep=""))
     
     ## Reading in the poplar fasta file and determining the nucleotide frequency
     poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
     freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
     freqBASE<-colSums(freqBASE[1:19,])
     probBASE<-freqBASE/sum(freqBASE)
     
     ## Reading in the filter info
     basesin<-read.table(paste(input.data.dir, "Ptricocarpa_bases_at_each_step.txt", sep=""), header=T)
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-as.numeric(basesin[,5])
     
     tree13bin<-basesin[which(basesin[,1] == "tree13" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     tree14bin<-basesin[which(basesin[,1] == "tree14" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     
     depth.vec<-tree13bin[,2]
     
     
     for (a in 1:length(depth.vec))
     {
       
       cat("Progress: ----------------------------------------------", a/length(depth.vec), "\n")
       
       gs.13<-tree13bin[which(tree13bin[,2] == depth.vec[a]), 5]
       gs.14<-tree14bin[which(tree14bin[,2] == depth.vec[a]), 5]
       gs.in<-mean(gs.13, gs.14)
       
       tree13.name<-paste("tree13.depth", depth.vec[a], "_hetAt1_step2", sep="")
       tree14.name<-paste("tree14.depth", depth.vec[a], "_hetAt1_step2", sep="")
       
       pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in, 
                                       input.dir = input.data.dir, 
                                       tree13 = tree13.name,
                                       tree14 = tree14.name)
       pedigree<-pedigree[[1]]
       pedigree<-pedigree[which(pedigree[,1] > 0),]
       
       ## Running models
       out<-mutSOMA(pedigree.data = pedigree, 
                    p0aa= probBASE[1],
                    p0cc= probBASE[2],
                    p0tt= probBASE[3],
                    p0gg= probBASE[4],
                    Nstarts=20,
                    prop.het = 0.4,
                    out.dir = out.data.dir,
                    out.name = paste("within_EST", tree13.name, tree14.name, sep="_"))
       
       bout<-bootSOMA(pedigree.data =out,    
                      Nboot=100, 
                      out.dir=out.data.dir, 
                      out.name= paste("within_BOOT", tree13.name, tree14.name, sep="_"))
       
     }    
     
     
     
     
     
     
######################################### 
#### depthX_hetAt1_step3
#### WITHIN TREE
#########################################
     rm(list=ls())
     library(Biostrings)
     library(vcfR)
     library(expm)
     library(optimx)
     
     input.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/RAW/"
     out.data.dir<-"/home/user/Dropbox/jlab/Projects/mutSOMA/DATA/PRODUCED/"
     func.dir<-"/home/user/Dropbox/jlab/Projects/FUNCTIONS/mutSOMA/"
     
     ## Loading source code
     source(paste(func.dir, "makePHYLO.R", sep=""))
     source(paste(func.dir, "makeVCFpedigreeTEMPv2.R", sep=""))
     source(paste(func.dir, "mutSOMA.R", sep=""))
     source(paste(func.dir, "bootSOMA.R", sep=""))
     
     ## Reading in the poplar fasta file and determining the nucleotide frequency
     poplar<-readDNAStringSet(paste(input.data.dir, "PtrichocarpaStettler14_532_v1.0.fa", sep=""))
     freqBASE<-letterFrequency(poplar, letters = c("A", "C", "T", "G"), as.prob=FALSE)
     freqBASE<-colSums(freqBASE[1:19,])
     probBASE<-freqBASE/sum(freqBASE)
     
     ## Reading in the filter info
     basesin<-read.table(paste(input.data.dir, "Ptricocarpa_bases_at_each_step.txt", sep=""), header=T)
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-sub(",", "", basesin[,5])
     basesin[,5]<-as.numeric(basesin[,5])
     
     tree13bin<-basesin[which(basesin[,1] == "tree13" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     tree14bin<-basesin[which(basesin[,1] == "tree14" & basesin[,3] == "Step0" & basesin[,4] == 1),]
     
     depth.vec<-tree13bin[,2]
     
     
     for (a in 1:length(depth.vec))
     {
       
       cat("Progress: ----------------------------------------------", a/length(depth.vec), "\n")
       
       gs.13<-tree13bin[which(tree13bin[,2] == depth.vec[a]), 5]
       gs.14<-tree14bin[which(tree14bin[,2] == depth.vec[a]), 5]
       gs.in<-mean(gs.13, gs.14)
       
       tree13.name<-paste("tree13.depth", depth.vec[a], "_hetAt1_step3", sep="")
       tree14.name<-paste("tree14.depth", depth.vec[a], "_hetAt1_step3", sep="")
       
       pedigree<-makeVCFpedigreeTEMPv2(genome.size=gs.in, 
                                       input.dir = input.data.dir, 
                                       tree13 = tree13.name,
                                       tree14 = tree14.name)
       pedigree<-pedigree[[1]]
       pedigree<-pedigree[which(pedigree[,1] > 0),]
       
       
       ## Running models
       out<-mutSOMA(pedigree.data = pedigree, 
                    p0aa= probBASE[1],
                    p0cc= probBASE[2],
                    p0tt= probBASE[3],
                    p0gg= probBASE[4],
                    Nstarts=20,
                    prop.het = 0.4,
                    out.dir = out.data.dir,
                    out.name = paste("within_EST", tree13.name, tree14.name, sep="_"))
       
       bout<-bootSOMA(pedigree.data =out,    
                      Nboot=100, 
                      out.dir=out.data.dir, 
                      out.name= paste("within_BOOT", tree13.name, tree14.name, sep="_"))
       
     }    
     
     
     