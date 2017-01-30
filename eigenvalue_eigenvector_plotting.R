#plotting eigenvectors
library(stringr)
library(ggplot2)
library(gridExtra)

plot_eigenvalues_vectors = function (wd='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/Hessian-based/',
                                     t_catchment = 'Gingera',sub_catchment,
                                     t_year = '70s',no_seed = 2025,perts = 1e-06  ){
  # prim_dir='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/'
  # sub_dir = 'Jacobian-based'
  # t_catchment = 'Gingera'
  # t_year = '70s'
  # no_seed = 2025
  # perts = 1e-06
  setwd(wd)
  
  ndir = str_c(t_catchment,'/',sub_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/')
  setwd(ndir)
  #read eigenvectors
  eigen_values = read.csv('eigenvalues.csv',header=F)
  eigen_vectors = read.csv('eigenvectors.csv',header=F) #column wise
  
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
  
  dir.create('figs',recursive = TRUE, showWarnings = F)
  
  tiff('figs/Eigenvectors.tiff',width=800)
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
  
  tiff('figs/Eigenvalues plot.tiff',width=800)
  print(p5)
  dev.off()

}


years = c('70s','80s','90s','00s')

for (year in years){
  
  plot_eigenvalues_vectors(wd='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/Hessian-based/',
                           t_catchment='Gingera',sub_catchment = 'original-range-1',
                           no_seed = 2025, t_year=year)

    
}

