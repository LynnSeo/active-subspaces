library('hydromad')
library('stringr')
library('xts')
load_tsPQE <- function(t_catchment='Gingera', t_year='70s'){
  name_tspqe=str_c('C:/UserData/seol/Sensitivity Analyses/PQE input/',t_catchment,'/',t_catchment,'.csv')
  tsPQE=read.zoo(name_tspqe,sep=',',header=TRUE);  tsPQE=as.xts(tsPQE)
  if(t_year=='70s'){
    assign('ts_t_year', tsPQE["1970-01-01::1979-12-31"],envir = .GlobalEnv)
  }else if(t_year=='80s'){
    assign('ts_t_year', tsPQE["1980-01-01::1989-12-31"],envir = .GlobalEnv)
  }else if(t_year=='90s'){
    assign('ts_t_year', tsPQE["1990-01-01::1999-12-31"],envir = .GlobalEnv)
  }else if(t_year=='00s'){
    assign('ts_t_year', tsPQE["2000-01-01::2009-12-31"],envir = .GlobalEnv)
  }
  return(ts_t_year)
}


#compute Jacobian and Hessian using finite differencing
est_J_H <- function (t_catchment = 'Gingera', t_year = '70s',
                     sub_catchment = 'test1', no_seed = 2025, perts = '1e-04', Nsamp=10){

  ts_t_year = load_tsPQE(t_catchment = t_catchment, t_year=t_year)
  #set up obj function
  # hydromad.options(objective = ~hmadstat("r.squared")(Q, X)/(2-hmadstat("r.squared")(Q, X)))
  hydromad.options(objective = hmadstat('RMSE'))
  
  # obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',
  #                 f=range(0.5,1.5),e=range(0.99,1.01),d=range(50,550),tau_q=range(0,10), tau_s=range(15,1000),  v_s=range(0,1))
  obj1 = hydromad(ts_t_year,sma='cmd',routing='expuh',f=range(0.5,1.5),e=range(0.99,1.01),d=range(200,550),tau_q=range(3,10), tau_s=range(200,1000),  v_s=range(0.7,1))
  # obj1_ts70s <- fitBySCE(obj1, control = list(trace = 1, ncomplex = 20))
  obj1$parlist['shape']=NULL
  
  print(obj1$parlist)
  
  #set up same seed number
  scale_pert = as.numeric(perts)
  set.seed(no_seed)
  print(no_seed)
  
  list_pert = lapply(obj1$parlist, function(x) (x[2]-x[1])*scale_pert  )
  df_pert = do.call('cbind',list_pert)
  
  # getFreeParsRanges(obj1)*scale_perts
  #generate random numbers on uniform distribution
  x0 = parameterSets(obj1$parlist,Nsamp,method='random'); set.seed(no_seed)
  x1 = as.data.frame(t(apply(x0, 1, function (x) x +  df_pert )))
  x2 = as.data.frame(t(apply(x0, 1, function (x) x -  df_pert )))
  Npar = length(names(x0))
  par_list = names(x0)
  colnames(x1) <- par_list
  colnames(x2) <- par_list
  
  #return jacobian and hessian
  j_h_deriv <- function (Nsamp, x0=x0, x1=x1, x2=x2){
    #x0 = x, x1 = x+h, x2 = x-h
    Jaco_mat = as.data.frame(matrix(NaN, Nsamp, Npar))
    time_taken = proc.time()
    #iterate as Nsamp
    for (k in 1:Nsamp){
      #temp_hess_mat = hessian at one point
      temp_h = as.data.frame(matrix(NaN,Npar,Npar))
      for (i in 1:Npar){  #row in Hessian
        for (j in i:Npar){  #col in Hessian
          p = as.data.frame(matrix(NaN, nrow=7, Npar ))
          colnames(p) = par_list;        
          p[1:7,] = x0[k,]
          p[1,i] =  x1[k,i];   p[1,j] =  x1[k,j]
          p[2,i] =  x1[k,i]
          p[3,i] =  x1[k,i];   p[3,j] =  x2[k,j]
          p[5,i] =  x2[k,i];   p[5,j] =  x1[k,j]
          p[6,i] =  x2[k,i]
          p[7,i] =  x2[k,i];   p[7,j] =  x2[k,j]
          u = evalPars(p, obj1)
          #u1 = u i+1, j+1        at p1
          #u2 = u i+1,j           at p2
          #u3 = u i+1,j-1         at p3
          #u4 = u i, j            at p4
          #u5 = u i-1, j+1        at p5
          #u6 = u i-1, j          at p6
          #u7 = u i-1, j-1        at p7
          # p : parameter values on grid
          # u : model output by parameter values defined above
          # i : 1st factor
          # j : 2nd factor 
          # Jacobian 
          # hessian returned has elements for half of it. coz it is symmetric. 
          if (i == j){
            # Jaco_mat[k, i] = (u[2]-u[6])/ (2*df_pert[i])
            Jaco_mat[k, i]= (u[2]-u[6])/ (2*df_pert[i])
            temp_h[i,j]  =  (u[2]-2*u[4]+u[6])/ ( 2*df_pert[i]^2 )
          }else{temp_h[i,j] =  (u[1]-u[3]-u[5] +u[7])/ (4* df_pert[i] * df_pert[j] )
          print (paste(k,i,j))}
        }
      }
      
      if (k==1){Hess_mat=temp_h
      }else{Hess_mat = cbind(Hess_mat, temp_h)}
      
    }
    time_taken2 = proc.time()- time_taken
    print (time_taken2)
    derivs = list(jacobian=Jaco_mat, hessian=Hess_mat)
    return(derivs)
  }
  
  
  ls_jh = j_h_deriv(Nsamp = 2, x0=x0, x1=x1, x2=x2)
  
  #generate dir : Gingera/70s/seed-2025/1e-06
  name_dir_path = str_c('Jacobian_Hessian/',t_catchment,'/',sub_catchment,'/',t_year,'/','seed-',no_seed,'/',perts,'/')
  dir.create(name_dir_path,recursive = TRUE, showWarnings = F)
  
  #print parameter range
  sink(str_c('Jacobian_Hessian/',t_catchment,'/',sub_catchment,'/parameter_range.txt'))
  print(paste('Nsamp =',Nsamp))
  print(obj1)
  sink()
  
  u0 = evalPars(x0, obj1)
  xy=as.data.frame(cbind(x0,u0))
  xy$shape = NULL
  
  write.csv(xy,str_c(name_dir_path,'xy.csv'),row.names = FALSE)
  write.csv(ls_jh$jacobian, str_c(name_dir_path,perts,'.csv'),row.names = FALSE)
  write.csv(ls_jh$hessian, str_c(name_dir_path,perts,'_Hessian.csv'),row.names = FALSE)
  print(str_c(name_dir_path))
  
  
}

#set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

perts = c('1e-02','1e-04','1e-06')
years = c('70s')
seeds = c(2025,2026)

for (p in perts){
  for (y in years){
    for (s in seeds){
      est_J_H( t_year = y,sub_catchment = 'constrained-range-1', no_seed = s, perts = p, Nsamp=1000)        
    }
  
  }
}


