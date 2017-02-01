#RS of Sacramento
library('zoo')
library('hydromad')
library('xts')
library('stringr')

#compute Jacobian using finite differencing
est_Jac = function (t_year='70s', no_seed=2025, t_catchment = 'Gingera',sub_catchment='Constrained range',
                    Nsamp=1000){
  name_tspqe=str_c('C:/UserData/seol/Sensitivity Analyses/PQE input/',t_catchment,'/',t_catchment,'.csv')
  tsPQE=read.zoo(name_tspqe,sep=',',header=TRUE)
  tsPQE=as.xts(tsPQE)
  head(tsPQE)
  #sub-decades
  ts70s <- tsPQE["1970-01-01::1979-12-31"]
  ts80s <- tsPQE["1980-01-01::1989-12-31"]
  ts90s <- tsPQE["1990-01-01::1999-12-31"]
  ts00s <- tsPQE["2000-01-01::2009-12-31"]
  
  
  #set up obj function
  # hydromad.options(objective = ~hmadstat("r.squared")(Q, X)/(2-hmadstat("r.squared")(Q, X)))
  # hydromad.options(objective = ~hmadstat("r.squared")(Q, X))
  hydromad.options(objective = hmadstat('RMSE'))
  
  
  #load model property including parameter ranges
  ts_t_year=get(str_c('ts',t_year))
  
  # obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',
  #                 f=range(0.5,1.5),e=range(0.99,1.01),d=range(50,550),tau_q=range(0,10), tau_s=range(15,1000),  v_s=range(0,1))
  # obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',
  #                 f=range(0.5,1.5),e=range(0.99,1.01),d=range(50,550),tau_q=range(2,10), tau_s=range(200,1000),  v_s=range(0.2,1))
  # obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',
  #                 f=range(0.5,1.5),e=range(0.99,1.01),d=range(200,550),tau_q=range(2,10), tau_s=range(200,1000),  v_s=range(0.6,1))
  obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',
                  f=range(0.5,1.5),e=range(0.99,1.01),d=range(200,550),tau_q=range(3,10), tau_s=range(200,1000),  v_s=range(0.7,1))
  
  obj1$parlist$shape=NULL
  #set up same seed number
  set.seed(no_seed)
  
  #generate random numbers on uniform distribution
  r0=matrix(data=runif(Nsamp*length(obj1$parlist),0,1),nrow=Nsamp,ncol=length(obj1$parlist))
  Kpar = ncol(r0)           #number of parameter
  Nsamp = nrow(r0)          #number of samples = 1,000
  r1=matrix(data=NA,nrow=Nsamp,ncol=Kpar)
  
  #perturb random numbers to 1e-7
  # perts=c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
  perts=c(1e-5,1e-6,1e-7)
  
  for (scale_pert in perts[2]){
    print(  scale_pert)
    
    #r0=x, r1=x+h
    for (j in 1:Kpar){
      r1[r0[,j]>0.5,j]=r0[r0[,j]>0.5,j]-scale_pert
      r1[r0[,j]<=0.5,j]=r0[r0[,j]<=0.5,j]+scale_pert
    }
    
    
    x0=matrix(NA,ncol=Kpar,nrow=Nsamp)
    x1=matrix(NA,ncol=Kpar,nrow=Nsamp)
    #generate input matrix to be run on Sacramento
    for (rsamp in 1:Nsamp){
      for (jpar in 1:Kpar){
        x0[rsamp,jpar]=qunif(r0[rsamp,jpar],obj1$parlist[[jpar]][1],obj1$parlist[[jpar]][2])
        x1[rsamp,jpar]=qunif(r1[rsamp,jpar],obj1$parlist[[jpar]][1],obj1$parlist[[jpar]][2])
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
    
    #bind cooridnates and output array
    xy=cbind(X,y)
    
    
    #caculate the derivatives
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
    
    deriv=as.data.frame(deriv)
    xy=as.data.frame(xy)
    
    deriv$shape = NULL
    xy$shape = NULL
    
    #create output folder
    name_dir_path = str_c(t_catchment,'/',sub_catchment,'/',t_year,'/','seed-',no_seed,'/',scale_pert,'/')
    dir.create(name_dir_path,recursive = TRUE)
    #print parameter range
    sink(str_c(t_catchment,'/',sub_catchment,'/',t_year,'/parameter_range.txt'))
    print(obj1)
    sink()
    
    #write derivates which is (Nsamp X Kpar) matrix
    write.csv(deriv,str_c(name_dir_path,scale_pert,'.csv'),row.names = F)
    #coordinates and RMSE
    write.csv(xy,str_c(name_dir_path,'xy.csv'),row.names = F)
    
    
    print(str_c(wd,name_dir_path))
  }
}

#set up working directory
wd='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/Jacobian-based/'
setwd(wd)

years = c('70s','80s','90s','00s')
seeds = c(2025,2026)


for (i in years){
  for (j in seeds){
    time_taken = proc.time()
    est_Jac(sub_catchment='constrained-range-1',Nsamp = 1000, t_year = i, no_seed = j)    
    time_taken2 = proc.time()- time_taken
    print (time_taken2)
    }
}
  


