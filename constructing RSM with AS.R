library(hydromad)
library(xts)
library(stringr)
library(ggplot2)
source('C:/UserData/seol/Sensitivity Analyses/Lynx custom functions/multiplot.R')

#construc RSM based on function g which consists of active variabnles and plot RSM
construct_RSM = function (prim_dir='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/',
                         sub_dir = 'Jacobian-based', t_catchment = 'Gingera/Constrained range',
                         t_year = '00s', no_seed = 2025, perts = 1e-06, seperation.n = 3){

# prim_dir='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/'
# sub_dir = 'Jacobian-based'
# t_catchment = 'Gingera/Constrained range'
# t_year = '70s'
# no_seed = 2025
# perts = 1e-06

ndir = str_c(prim_dir,sub_dir,'/',t_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/')
print (ndir)
setwd(ndir)
print (getwd())
#read coordinates
xy = read.csv('xy.csv')
#delete the list of ids 
xy$X=NULL
head(xy)

#original model output
q = xy[,ncol(xy)]
xy[,ncol(xy)] = NULL

#original paramaeter set
# m is number of dimension
m = ncol(xy)
Nsamp = nrow(xy)/(m+1)


x= xy # Nsamp by m 
x=t(x) # m by Nsamp

#read eigenvectors
ei_vec = read.csv('eigenvectors.csv',header=F)


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

write.csv(z,'y.csv',row.names = F)

mydata = as.data.frame(cbind(t(z),q))
head(mydata)
#mydata2 = df with sqaured elements. 
mydata2= mydata[,1:seperation.n]^2
head(mydata2)

colnames(mydata)[1:seperation.n]=c('theta1','theta2','theta3','theta4','theta5','theta6')[1:seperation.n]
head(mydata)

if (seperation.n==3){
  fit <- lm(q ~ theta1 +theta2 +theta3, data=mydata)
}else{fit <- lm(q ~ theta1 +theta2 +theta3 +theta4, data=mydata)    
}

#multi-dimensional regression
a = fit$coefficients[(1+1):(seperation.n+1)]
write.csv(a,'A.csv',row.names=F)
#function g values
predicts = a %*% z  + fit$coefficients[1]

const_g = fit$coefficients[1]

result =list(z, const_g)


#r = errors, residuals
r = q - predicts
write.csv(r,'r.csv',row.names = F)
write.csv(predicts,'g.csv',row.names = F)

#generate matrix for plot
matrix_plot = cbind(t(x),q[1:1000],predicts[1:1000],r[1:1000],t(z))

colnames(matrix_plot)[(m+1):(m+3)] = c('original','fitted','residual')
colnames(matrix_plot)[(m+4):(m+3+seperation.n)] = c('y1','y2','y3','y4')[1:seperation.n]
head(matrix_plot)

matrix_plot=as.data.frame(matrix_plot)
library(gridExtra)

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
  p1 <- ggplot(matrix_plot,aes_string(name_para,'original'))+geom_point(alpha = val_alpha)+ylab('RMSE')+ggtitle('f(x)')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
    annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.original)) )+
    # annotate('text', size=7, x=x.location.text, y=-0.5, label= var.original ) + 
    theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis,angle=90),
          axis.title.x = element_text(size=(size.font+10)))
  p2 <- ggplot(matrix_plot,aes_string(name_para,'fitted'))+geom_point(alpha = val_alpha)+ggtitle('g(y)') +ylab('RMSE')+xlab('1')+  xlab(name_para)+ylim(range.ylim)+
    annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.fitted)) )+
    theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
          axis.title.x = element_text(size=(size.font+10)))
  p3 <- ggplot(matrix_plot,aes_string(name_para,'residual'))+geom_point(alpha = val_alpha) +ggtitle('f - g')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
    annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.resi)) )+
    theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
          axis.title.x = element_text(size=(size.font+10)))
  
  tiff(str_c('figs/',name_para,'_RS.png'), width=800)
  grid.arrange(p1,p2,p3, ncol=3, nrow =1)
  dev.off()
  print(  grid.arrange(p1,p2,p3, ncol=3, nrow =1))
}

dir.create('figs',recursive = TRUE, showWarnings = F)

para.array = c('f','e','d','tau_q','tau_s','v_s','y1','y2','y3')
for (i in para.array){
  plot_original_para(i)
}

}


# construct_RSM(sub_dir = 'Hessian-based', t_catchment = 'Gingera/constrained-range-1',
#               t_year = '70s', no_seed = 2025,
#               seperation.n =3)


# years=c('70s','80s','90s')
# seeds = c(2025, 2026)
# 
# for (year in years){
#   for (seed in seeds){
#     construct_RSM(sub_dir = 'Hessian-based', t_catchment = 'Gingera/constrained-range-1',
#                   t_year = year, no_seed = seed,
#                   seperation.n =3)
# 
#   }
# }

construct_RSM(sub_dir = 'Jacobian-based', t_catchment = 'Gingera/constrained-range-1',
              t_year = '70s', no_seed = 2025,
              seperation.n =3)



