

bootSOMA<-function(pedigree.data, Nboot, out.dir, out.name)
{

  pedigree.data<-out
  
  ## Reading the dataset for bootstrapping and extracting the parameter settings
  settings<-pedigree.data$settings
  est<-pedigree.data$estimates
  optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
  pedigree<-pedigree.data$pedigree
  

  ##### Defining the divergence function
  divergence <- function(pedigree, p0aa, p0cc, p0tt, p0gg, param)
  {
    
    ## Initializing parameters
    # c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CA", "CT", "CG", "TA", "TC", "TG", "GA", "GC", "GT")
    PrAA <- p0aa; PrCC <- p0cc; PrTT <- p0tt; PrGG <- p0gg; PrAC <- p0ac; PrAT <- p0at; PrAG <- p0ag; PrCA <- p0ca; PrCT <- p0ct
    PrCG <- p0cg; PrTA <- p0ta; PrTC <- p0tc; PrTG <- p0tg; PrGA <- p0ga; PrGC <- p0gc; PrGT <- p0gt
    
    g <- param[1]
    
    
    ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    ## Keeping this in this format because I may want to introduce weights for the "unobserved" initial heterozygotes
    svGzero    <- c(PrAA,
                    PrCC,
                    PrTT,
                    PrGG,
                    PrAC,
                    PrAT,
                    PrAG,
                    PrCA,
                    PrCT,
                    PrCG,
                    PrTA,
                    PrTC,
                    PrTG,
                    PrGA,
                    PrGC,
                    PrGT)
    
    ## Defining the generation (or transition) matrix for the mitotic case
    Tmat<-matrix(NA, nrow=16, ncol=16)
    
    Ta<-matrix(c(
      (1-g)^2, 1/9*g^2, 1/9*g^2, 1/9*g^2, 
      1/9*g^2, (1-g)^2, 1/9*g^2, 1/9*g^2, 
      1/9*g^2, 1/9*g^2, (1-g)^2, 1/9*g^2, 
      1/9*g^2, 1/9*g^2, 1/9*g^2, (1-g)^2), nrow=4, byrow=TRUE)
    
    
    #              1                2            3          4            5           6              7           8            9            10           11           12
    Tb<-matrix(c(
      1/3*(1-g)*g,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
      1/3*(1-g)*g,   1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
      1/9*g^2,       1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
      1/9*g^2,       1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,     1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g), 
      nrow=4, byrow=TRUE)
    
    
    Tc<-matrix(c(
      1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
      1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
      1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
      1/3*(1-g)*g,  1/3*(1-g)*g, 1/9*g^2,     1/9*g^2, 
      1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2, 
      1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
      1/3*(1-g)*g,  1/9*g^2,     1/3*(1-g)*g, 1/9*g^2,
      1/9*g^2,      1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*g^2,
      1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
      1/3*(1-g)*g,  1/9*g^2,     1/9*g^2,     1/3*(1-g)*g,
      1/9*g^2,      1/3*(1-g)*g, 1/9*g^2,     1/3*(1-g)*g,
      1/9*g^2,      1/9*g^2,     1/3*(1-g)*g, 1/3*(1-g)*g), nrow=12, byrow=TRUE)
    
    #              1                2            3           4            5           6            7            8           9            10           11           12
    Td<-matrix(c(
      (1-g)^2,     1/3*(1-g)*g,  1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
      1/3*(1-g)*g, (1-g)^2,      1/3*(1-g)*g,  1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
      1/3*(1-g)*g, 1/3*(1-g)*g,  (1-g)^2,      1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
      1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
      1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,
      1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,
      1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,
      1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,
      1/9*(g)^2,   1/9*(g)^2,    1/3*(1-g)*g,  1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2,     1/9*(g)^2,   1/9*(g)^2, 1/9*(g)^2,
      1/9*(g)^2,   1/9*(g)^2,    1/9*(g)^2,    1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   (1-g)^2,     1/3*(1-g)*g, 1/3*(1-g)*g,
      1/3*(1-g)*g, 1/9*(g)^2,    1/9*(g)^2,    1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,  1/3*(1-g)*g, (1-g)^2,     1/3*(1-g)*g, 
      1/9*(g)^2,   1/3*(1-g)*g,  1/9*(g)^2,    1/9*(g)^2,   1/3*(1-g)*g, 1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/9*(g)^2,   1/3*(1-g)*g, 1/3*(1-g)*g, (1-g)^2), nrow=12, byrow=TRUE)
    
    
    ## Transition matrix
    Tmat[1:4, 1:4]<-Ta
    Tmat[1:4, 5:16]<-Tb
    Tmat[5:16, 1:4]<-Tc
    Tmat[5:16, 5:16]<-Td
    Genmatrix <- Tmat
    
    ## Coding the divergence between samples i and j
    Deffects<-matrix(as.numeric(c("0", "1", "1", "1", "0.5", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", 
                                  "0.5", "0.5", "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "0.5", "0.5", 
                                  "1", "1", "0.5", "1", "1", "1", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "1", 
                                  "1", "0", "0.5", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "0.5", "0", "0.5", "1", "0.5", 
                                  "1", "1", "1", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                  "1", "0.5", "0.5", "1", "1", "1", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", 
                                  "0.5", "1", "0.5", "0", "0.5", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "1", "0.5", "1", "1", "0.5", "0.5", "0.5", "0", "1", 
                                  "1", "0.5", "1", "1", "1", "0.5", "1", "0.5", "1", "1", "1", "1", "0.5", "1", "1", "0", "0.5", "0.5", "0.5", "1", "1", "1", 
                                  "0.5", "0.5", "1", "0.5", "1", "1", "1", "1", "1", "0.5", "0", "0.5", "1", "0.5", "1", "1", "1", "0.5", "0.5", "1", "1", "0.5", 
                                  "1", "1", "0.5", "0.5", "0.5", "0", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "0.5", "1", "1", "0.5", "1", "1", 
                                  "0", "0.5", "0.5", "1", "0.5", "1", "0.5", "0.5", "1", "1", "1", "1", "1", "1", "0.5", "1", "0.5", "0", "0.5", "1", "1", "0.5", 
                                  "0.5", "1", "0.5", "1", "1", "0.5", "1", "1", "1", "1", "0.5", "0.5", "0")), nrow=16, ncol=16, byrow=TRUE)
    
    
    
    ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL
    
    ## Matrix of all possible initial states
    tvec.mat<-diag(16)
    
    
    for (p in 1:nrow(pedigree))
    {
      
      svt0<-NULL
      svt1<-list()
      svt2<-list()
      
      svt0<-t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      
      for (i in 1:nrow(tvec.mat))
      {
        
        svt1[i]<-list(t(tvec.mat[1,])%*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1])))
        svt2[i]<-list(t(tvec.mat[1,]) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1])))
      }
      
      
      
      DivProbt1t2<-NULL
      
      for (j in 1:nrow(tvec.mat))
      {
        
        t1in<-svt1[[j]]
        t2in<-svt2[[j]]
        
        jointPROB     <- expand.grid(a = t1in, b = t2in)
        jointPROB     <- jointPROB[,1] * jointPROB[,2] 
        jointPROBt1t2 <- matrix(jointPROB, 16, 16, byrow=F)
        DivProbt1t2[j]   <- svt0[j]*sum(Deffects * jointPROBt1t2)
      }
      
      
      
      ## Total (weighted) divergence 
      Dt1t2[p]<- sum(DivProbt1t2)
      
      
    }
    
    divout<-list(Dt1t2)
    
    return(divout)
    
  }
  
  
  ###### Defining the Least Square function to be minimized
  ###### Note the equilibrium constraint, which can be made as small as desired.
  
  LSE_intercept<-function(param_int) 
  {
    sum((pedigree[,4] - param_int[2] - divergence(pedigree, p0aa, p0cc, p0tt, p0gg, param_int[1])[[1]])^2)	
  }
  
  
  
  ###### Calculating the initial proportions 
  p0aa<-as.numeric(as.character(settings[which(as.character(settings[,1]) == "p0aa"),2]))
  p0cc<-as.numeric(as.character(settings[which(as.character(settings[,1]) == "p0cc"),2]))
  p0tt<-as.numeric(as.character(settings[which(as.character(settings[,1]) == "p0tt"),2]))
  p0gg<-as.numeric(as.character(settings[which(as.character(settings[,1]) == "p0gg"),2]))
  #p0aa<-p0aain/sum(p0aain, p0ccin, p0ttin, p0ggin)
  #p0cc<-p0ccin/sum(p0aain, p0ccin, p0ttin, p0ggin)
  #p0tt<-p0ttin/sum(p0aain, p0ccin, p0ttin, p0ggin)
  #p0gg<-p0ggin/sum(p0aain, p0ccin, p0ttin, p0ggin)
  prop.het<-as.numeric(as.character(settings[which(as.character(settings[,1]) == "prop.het"),2]))
  p0ac <- (prop.het)*1/12
  p0at <- (prop.het)*1/12
  p0ag <- (prop.het)*1/12
  p0ca <- (prop.het)*1/12
  p0ct <- (prop.het)*1/12
  p0cg <- (prop.het)*1/12
  p0ta <- (prop.het)*1/12
  p0tc <- (prop.het)*1/12
  p0tg <- (prop.het)*1/12
  p0ga <- (prop.het)*1/12
  p0gc <- (prop.het)*1/12
  p0gt <- (prop.het)*1/12
  
  if(as.numeric(as.character(sum(c(p0aa, p0cc, p0tt, p0gg, p0ac, p0at, p0ag, p0ca, p0ct, p0cg, p0ta, p0tc, p0tg, p0ga, p0gc, p0gt), na.rm =T))) != 1) 
  {stop("The initial state probabilities don't sum to 1")}
  
  
  ##### Initializing
  optim.method<-as.character(settings[which(as.character(settings[,1]) == "optim.method"),2])
  final<-NULL
  counter<-0
  opt.out<-NULL
  
  ## Defining starting values
  gamma.start  <-est[1,1]
  intercept.start <-est[1,2]
  param_int0 = c(gamma.start, intercept.start)
  
  
  ## Start of boostrap loops
  for (booting in 1:Nboot)
  {
    
    opt.out<-NULL
    pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=T)
    
    counter<-counter+1
    
    cat("Bootstrap interation: ", counter/Nboot, "\n")
    
    opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
    opt.out <-cbind(opt.out, gamma.start, intercept.start)
    final<-rbind(final, opt.out)
    
    
  } # End of Bootstrap loops
  
  colnames(final)[1:2]<-c("gamma", "intercept")
  
  SE.gamma<-sd(final[,1],na.rm=TRUE)
  SE.inter<-sd(final[,2],na.rm=TRUE)
  CI.gamma<-quantile(final[,1],probs=c(0.025, 0.975))
  CI.inter<-quantile(final[,2],probs=c(0.025, 0.975))
 
  
  SE<-c(SE.gamma, SE.inter) 
  CI<-rbind(CI.gamma, CI.inter) 
  
  SE.out<-cbind(SE, CI)
  colnames(SE.out)[1]<-"SE"
  rownames(SE.out)<-c("gamma", "intercept")
  
  final<-data.frame(final)
  good.boots<-length(which(is.na(final[,"gamma"]) == FALSE))
  
  SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final)
  names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results")
  
  ## Ouputting result datasets
  dput(SE.out, paste(out.dir, out.name, ".Rdata", sep=""))
  
  return(SE.out)
  

    

  
} #End of function
	
		
		
		
		
		
		
		