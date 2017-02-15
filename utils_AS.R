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
                     sub_catchment = '1sample', no_seed = 2025, perts = '1e-04', Nsamp=1, method='random'){
  
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
  x0 = parameterSets(obj1$parlist,Nsamp,method=method); set.seed(no_seed)
  x1 = as.data.frame(t(apply(x0, 1, function (x) x +  df_pert )))
  x2 = as.data.frame(t(apply(x0, 1, function (x) x -  df_pert )))
  Npar = length(names(x0))
  par_list = names(x0)
  colnames(x1) <- par_list
  colnames(x2) <- par_list
  
  #return jacobian and hessian
  j_h_deriv <- function (Nsamp, Npar, x0=x0, x1=x1, x2=x2){
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
          }else{temp_h[i,j] =  (u[1]-u[3]-u[5] +u[7])/ (4* df_pert[i] * df_pert[j] );temp_h[j,i]=temp_h[i,j]
          print (paste(k,i,j))}
        }
      }
      
      if (k==1){Hess_mat=temp_h
      u0=u[4]
      }else{Hess_mat = rbind(Hess_mat, temp_h)
      u0=rbind(u0,u[4])
      }
      
    }
    time_taken2 = proc.time()- time_taken
    print (time_taken2)
    
    derivs = list(jacobian=Jaco_mat, hessian=Hess_mat,u0=u0)
    return(derivs)
  }
  
  ls_jh = j_h_deriv(Nsamp = Nsamp, Npar=Npar, x0=x0, x1=x1, x2=x2)
  
  #generate dir : Gingera/70s/seed-2025/1e-06
  name_dir_path = str_c('Jacobian-Hessian/',t_catchment,'/',sub_catchment,'/',t_year,'/','seed-',no_seed,'/',perts,'/')
  dir.create(name_dir_path,recursive = TRUE, showWarnings = F)
  
  #print parameter range
  sink(str_c('Jacobian-Hessian/',t_catchment,'/',sub_catchment,'/parameter_range.txt'))
  print(paste('Nsamp =',Nsamp))
  print(obj1)
  sink()
  
  xy=as.data.frame(cbind(x0,ls_jh$u0, row.names=NULL))
  
  xy$shape = NULL
  
  write.csv(xy,str_c(name_dir_path,'xy.csv'),row.names = FALSE)
  write.csv(ls_jh$jacobian, str_c(name_dir_path,perts,'.csv'),row.names = FALSE)
  write.csv(ls_jh$hessian, str_c(name_dir_path,perts,'_Hessian.csv'),row.names = FALSE)
  print(str_c(name_dir_path))
  
  
}

plot_eigen = function (t_catchment = 'Jacobiand-Hessian/Gingera', sub_catchment='test1',
                                     t_year = '70s',no_seed = 2025,perts = '1e-04' ,additional_dir='' ){
  # prim_dir='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/'
  # sub_dir = 'Jacobian-based'
  # t_catchment = 'Gingera'
  # t_year = '70s'
  # no_seed = 2025
  # perts = 1e-06
  
  ndir = str_c(t_catchment,'/',sub_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/',additional_dir)
  
  #read eigenvectors
  eigen_values = read.csv(str_c(ndir,'eigenvalues.csv'),header=F)
  eigen_vectors = read.csv(str_c(ndir,'eigenvectors.csv'),header=F) #column wise
  
  Npar = nrow(eigen_vectors)
  ncol(eigen_vectors)
  df_values = as.data.frame(matrix(NaN,Npar,Npar))
  df_values = as.data.frame(eigen_vectors)
  colnames(df_values) = c('vector1','vector2','vector3','vector4','vector5','vector6')
  #index for orinal parameters
  df_values$Index = c('1','2','3','4','5','6')
  #reorder levels
  df_values$Index = factor(df_values$Index, levels = df_values$Index[1:Npar])
  
  head(df_values)
  p=ggplot(data = df_values)
  size.font = 15
  size.font.xaxis = 15
  
  p1 = p +geom_point(aes(x=Index, y=vector1),size=5)+ylab('Active Variable Weight')+xlab('Input Parameter')+
    theme(text = element_text(size=size.font),axis.text.x = element_text(size=size.font.xaxis)) 
  p2 = p +geom_point(aes(x=Index, y=vector2),size=5)+ylab('Active Variable Weight')+xlab('Input Parameter')+
    theme(text = element_text(size=size.font),axis.text.x = element_text(size=size.font.xaxis)) 
  p3 = p +geom_point(aes(x=Index, y=vector3),size=5)+ylab('Active Variable Weight')+xlab('Input Parameter')+
    theme(text = element_text(size=size.font),axis.text.x = element_text(size=size.font.xaxis)) 
  p4 = p +geom_point(aes(x=Index, y=vector4),size=5)+ylab('Active Variable Weight')+xlab('Input Parameter')+
    theme(text = element_text(size=size.font),axis.text.x = element_text(size=size.font.xaxis)) 
  
  dir.create(str_c(ndir,'/figs'),recursive = TRUE, showWarnings = F)
  
  tiff(str_c(ndir,'/figs/Eigenvectors.tiff'),width=800)
  grid.arrange(p1,p2,p3,p4, ncol=2, nrow =2)
  dev.off()
  
  #plot eigenvalues
  df_ev = eigen_values
  df_ev$index = seq(Npar)
  head(df_ev)
  colnames(df_ev) = c('Eigenvalues','Index')
  
  # p5 =p5 + geom_point(aes(x=index, y=V1),size=5)+ylab('Eigenvalues')+xlab('Active variable index')+
  #   theme(text = element_text(size=size.font),axis.text.x = element_text(size=size.font.xaxis))+scale_y_log10() 
  p5 = ggplot(df_ev, aes(x=Index, y=Eigenvalues))+geom_line(size=2) +geom_point(size=6) +
    # scale_x_discrete(breaks=seq(6))+
    scale_y_log10(breaks = 10^seq(10,-5)) + 
    theme(text = element_text(size=20),axis.text.x = element_text(size=20))+xlab('Index of linear combinations of parameters') 
  
  tiff(str_c(ndir,'figs/Eigenvalues plot.tiff'),width=800)
  print(p5)
  dev.off()
  
}

# construc RSM based on function g which consists of active variabnles and plot RSM
# high-dimensioanl linear fitting used
# g(y) =a W1T X
construct_RSM = function (t_catchment = 'Jacobian-Hessian/Gingera', sub_catchment='test1',
                          t_year = '70s', no_seed = 2025, perts = '1e-04', seperation.n = 3, plot=FALSE,
                          additional_dir=''){
  ndir = str_c(t_catchment,'/',sub_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/',additional_dir)
  #read coordinates
  xy = read.csv(str_c(ndir,'xy.csv'))
  #delete the list of ids 
  xy$X=NULL
  head(xy)
  #original model output
  q = xy[,ncol(xy)]
  xy[,ncol(xy)] = NULL
  #original paramaeter set
  # m is number of dimension
  m = ncol(xy)
  Nsamp = nrow(xy)
  x= xy # Nsamp by m 
  x=t(x) # m by Nsamp
  
  #read eigenvectors
  ei_vec = read.csv(str_c(ndir,'eigenvectors.csv'),header=F)
  
  
  #separate eigenvectors to active subspaces and inactive subspaces
  # seperation.n = 3 # n 
  # Nsamp = 6
  
  w1t = t(ei_vec[,1:seperation.n])  # n by m
  w2t =t(ei_vec[,(seperation.n+1):m]) # m-n by m
  
  # new parametger sets = theta
  # n by Nsamp = n by m * m by Nsamp
  w1t = data.matrix(w1t)
  # z = y, active variables
  z = w1t %*% x
  
  write.csv(t(z),str_c(ndir,'y.csv'),row.names = F)
  
  mydata = as.data.frame(cbind(t(z),q))
  head(mydata)
  #mydata2 = df with sqaured elements. 
  mydata2= mydata[,1:seperation.n]^2
  head(mydata2)
  
  colnames(mydata)[1:seperation.n]=c('theta1','theta2','theta3','theta4','theta5','theta6')[1:seperation.n]
  head(mydata)
  
  if (seperation.n==3){
    fit <- lm(q ~ theta1 +theta2 +theta3, data=mydata)
  }else if(seperation.n==1){fit <- lm(q ~ theta1, data=mydata)    
  }else if(seperation.n==2){fit <- lm(q ~ theta1 +theta2, data=mydata)    
  }else if(seperation.n==4){fit <- lm(q ~ theta1 +theta2 +theta3+theta4, data=mydata)    
  }
  
  
  #multi-dimensional regression
  a = fit$coefficients[(1+1):(seperation.n+1)]
  write.csv(a,str_c(ndir,'A.csv'),row.names=F)
  #function g values
  predicts = a %*% z  + fit$coefficients[1]
  const_g = fit$coefficients[1]
  result =list(z, const_g)
  
  #r = errors, residuals
  r = q - predicts
  write.csv(t(r),str_c(ndir,'r.csv'),row.names = F)
  write.csv(t(predicts),str_c(ndir,'g.csv'),row.names = F)
  
  #generate matrix for plot
  matrix_plot = cbind(t(x),q[1:Nsamp],predicts[1:Nsamp],r[1:Nsamp],t(z))
  
  colnames(matrix_plot)[(m+1):(m+3)] = c('original','fitted','residual')
  colnames(matrix_plot)[(m+4):(m+3+seperation.n)] = c('y1','y2','y3','y4')[1:seperation.n]
  head(matrix_plot)
  
  matrix_plot=as.data.frame(matrix_plot)
  
  
  #plotting options
  range.ylim = c(-0.5,2)
  size.font = 20
  size.font.xaxis = 20
  
  #calculate Std of scatters
  var.original = round(sd(matrix_plot$original)^2, digits=4)
  var.fitted = round(sd(matrix_plot$fitted)^2, digits=4)
  var.resi = round(sd(matrix_plot$residual)^2, digits=4)
  
  #plotting function
  plot_original_para <- function (name_para){
    if (is.numeric(name_para)==TRUE){
      name_para = colnames(matrix_plot)[name_para]
    }
    x.location.text = (max(matrix_plot[name_para])-min(matrix_plot[name_para]))*0.5+min(matrix_plot[name_para])
    val_alpha = 1/10
    p1 <- ggplot(matrix_plot,aes_string(name_para,'original'))+geom_point(alpha = val_alpha)+ylab('output')+ggtitle('f(x)')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
      annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.original)) )+
      # annotate('text', size=7, x=x.location.text, y=-0.5, label= var.original ) + 
      theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis,angle=90),
            axis.title.x = element_text(size=(size.font+10)))
    p2 <- ggplot(matrix_plot,aes_string(name_para,'fitted'))+geom_point(alpha = val_alpha)+ggtitle('g(y)') +ylab('output')+xlab('1')+  xlab(name_para)+ylim(range.ylim)+
      annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.fitted)) )+
      theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
            axis.title.x = element_text(size=(size.font+10)))
    p3 <- ggplot(matrix_plot,aes_string(name_para,'residual'))+geom_point(alpha = val_alpha) +ggtitle('f - g')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
      annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.resi)) )+
      theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
            axis.title.x = element_text(size=(size.font+10)))
    
    tiff(str_c(ndir,'figs/',name_para,'_RS.png'), width=800)
    grid.arrange(p1,p2,p3, ncol=3, nrow =1)
    dev.off()
    print(  grid.arrange(p1,p2,p3, ncol=3, nrow =1))
  }
  
  if (isTRUE(plot)){
    
    dir.create(str_c(ndir,'figs'),recursive = TRUE, showWarnings = F)
    para.array = c('f','e','d','tau_q','tau_s','v_s','y1','y2','y3')[1:(m+seperation.n)]
    for (i in para.array){
      plot_original_para(i)
    }
    
  }
  
}

