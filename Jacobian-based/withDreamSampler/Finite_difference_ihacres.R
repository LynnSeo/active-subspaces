#RS of Sacramento
library('zoo')
library('hydromad')
library('xts')
library('stringr')
library(sensitivity)
#set up working directory
wd='C:/UserData/seol/Sensitivity Analyses/IHACRES/Gradient using finite difference/withDreamSampler/'
setwd(wd)

#read input
#read input
####################################
t_year='80s'           #target year for calibration
t_catchment='Murrindindi'  #target catchment 
####################################
name_tspqe=str_c('C:/UserData/seol/Sensitivity Analyses/PQE input/',t_catchment,'/',t_catchment,'.csv')
tsPQE=read.zoo(name_tspqe,sep=',',header=TRUE)
tsPQE=as.xts(tsPQE)
#sub-decades
ts80s <- tsPQE["1980-01-01::1989-12-31"]
ts90s <- tsPQE["1990-01-01::1999-12-31"]
ts00s <- tsPQE["2000-01-01::2009-12-31"]

years=c('80s','90s','00s')

which_year = years[1]
no_seed = 2026
  
  
  #set up obj function
  # hydromad.options(objective = ~hmadstat("r.squared")(Q, X)/(2-hmadstat("r.squared")(Q, X)))
  hydromad.options(objective = hmadstat('RMSE'))
  
  #load model property including parameter ranges
  ts_which_year=get(str_c('ts',which_year))
  # obj1  = hydromad(ts_which_year,sma='sacramento',routing=NULL,uztwm=c(15,150))
  obj1 = hydromad(ts_which_year,sma='cmd',routing='expuh',
                  f=range(0.5,1.5),e=range(0.99,1.01),d=range(50,550),tau_q=range(0,10), tau_s=range(15,1000),  v_s=range(0,1))
  
  obj1$parlist$shape=NULL
  #set up same seed number
  
  set.seed(no_seed)
  
  #read original Dream samples
  sampled.dream = read.csv('sampled.by.dream.csv')
  head(sampled.dream)
  sampled.dream$X =NULL
  
  #separate cooridinates and output
  coordinate.dream.sample = sampled.dream[,1:6]
  
  #caculate RMSE
  hydromad.getOption('objective')
  hydromad.options(objective = hmadstat('RMSE'))
  y.dream = evalPars(coordinate.dream.sample,obj1)
  
  #check values
  summary(y.dream)
  
  
  pts.dream = cbind(coordinate.dream.sample,y.dream)
  
  write.csv(pts.dream,'origianl_dream_sample.csv')
  # summary(selected.sampled.dream)
  first.pt = 1
  last.pt = 2000
  
  plot(pts.dream[first.pt:last.pt,1],pts.dream[first.pt:last.pt,7],ylim = c(0.3,1.8))
  
  tiff('dreamSamples.tiff')
  plot(pts.dream[first.pt:last.pt,1],pts.dream[first.pt:last.pt,7],ylim = c(0.3,1.8))
  dev.off()
  
  new.set = pts.dream[first.pt:last.pt,]
  
  #linear fitting
  a = lm(new.set[,7] ~ new.set[,1])
  
  #include fitted values and residuals
  new.set = cbind(new.set,a$fitted.values,a$residuals)
  
  colnames(new.set)[7] = 'y'
  colnames(new.set)[8] = 'fitted.value'
  colnames(new.set)[9] = 'residuals'
  
  new.set2 = new.set[order(new.set$residuals),]
  head(new.set2)
  
  new.set3 = new.set2[1:1000,]
  head(new.set3)
  
  tiff('DreamSample-afterScreen.tiff')
  plot(new.set3[,1],new.set3[,7])
  dev.off()
  
  write.csv(new.set3,'screened.dream.sample.csv')
  
  
  #generate random numbers on uniform distribution
  x0=matrix(NA,ncol=Kpar,nrow=Nsamp)
  x0= new.set3[,1:6]
  
  Kpar = ncol(x0)           #number of parameter
  Nsamp = nrow(new.set3)          #number of samples = 1,000
  
  #perturb random numbers to 1e-7
  # perts=c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
  perts=c(1e-5,1e-6,1e-7)
  scale_pert = perts[2]
  
    
    x1=matrix(NA,ncol=Kpar,nrow=Nsamp)
    #generate input matrix to be run on Sacramento
    
    for (rsamp in 1:Nsamp){  
      for (jpar in 1:Kpar){
    x1[rsamp,jpar] = x0[rsamp,jpar] +(obj1$parlist[[jpar]][2]-obj1$parlist[[jpar]][1])*scale_pert    
        
      }
    }
    
    X2 = do.call(rbind, lapply(1:ncol(x0), function(i) {
      X2i = x0
      X2i[, i] = x1[, i] 
      X2i
    }))
    
    
    #final input matrix (sample coordinates) to be run on a model
    X=rbind(x0,X2)
    
    colnames(X)=names(obj1$parlist)#attach parlist header
    head(X)

    #run model and return model performance measures.
    y=evalPars(X,obj1)
    
    xy=cbind(X,y)
    
    
    #all the derivatives
    out <- as.numeric(y)
    deriv=matrix(NA, ncol = Kpar, nrow = Nsamp)
    
    for (rsamp in 1:Nsamp){
      for (jpar in 1:Kpar){
        idx.pert=Nsamp*jpar+rsamp
        deriv[rsamp, jpar] = (out[idx.pert] - out[rsamp])/(X[idx.pert, jpar] - X[rsamp, jpar])    
      }
    }
    
    
    #attach header which indicates parameter names on the result matrix
    colnames(deriv)=names(obj1$parlist)
    
    summary(deriv)
    summary(xy)
    
    deriv=as.data.frame(deriv)
    xy=as.data.frame(xy)
    
    deriv$shape = NULL
    xy$shape = NULL
    
    dir.create(str_c(which_year,'/','seed-',no_seed,'/',scale_pert),recursive = TRUE)
    
    
    #write derivates which is 1000X13 (Nsamp X Kpar) bvmatrix
    write.csv(deriv,str_c(wd,which_year,'/','seed-',no_seed,'/',scale_pert,'/',scale_pert,'.csv'))
    #coordinates and RMSE
    write.csv(xy,str_c(wd,which_year,'/','seed-',no_seed,'/',scale_pert,'/','xy.csv'))
    print(no_seed)
    print(which_year)
    print(str_c(which_year,'/','seed-',no_seed,'/',scale_pert))
  


