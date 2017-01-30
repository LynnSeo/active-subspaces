library(hydromad)
library(xts)
library(stringr)
wd1='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/RSM with AS/Constrained-e/'
sub_dir ='80s/'
setwd(str_c(wd1,sub_dir))
getwd()
#read coordinates
xy = read.csv(str_c('seed-2025/1e-06/','xy.csv'))
#delete the list of ids 
xy$X=NULL
head(xy)

#original model output
q = xy$y
xy$y = NULL

#original paramaeter set
# m is number of dimension
Nsamp = nrow(xy)

x= xy # Nsamp by m 
x=t(x) # m by Nsamp

#read eigenvectors
ei_vec = read.table(str_c('figs/','eigenvectors.txt'))


#separate eigenvectors to active subspaces and inactive subspaces
seperation.n = 4 # n 
Nsamp = 6

w1t = ei_vec[1:seperation.n,]  # n by m
w2t =ei_vec[seperation.n+1:6,] # m-n by m

# new parametger sets = theta
# n by Nsamp = n by m * m by Nsamp
w1t = data.matrix(w1t)
z = w1t %*% x

mydata1 = as.data.frame(cbind(t(z),q))
head(mydata)

mydata2 = mydata1[,1:seperation.n]^2
head(mydata2)

mydata = cbind(mydata2,mydata1)
head(mydata)

colnames(mydata)[1:8]=c('sq1','sq2','sq3','sq4','theta1','theta2','theta3','theta4')
head(mydata)

fit <- lm(q ~ sq1+sq2+sq3+sq4+theta1 +theta2 +theta3 +theta4, data=mydata)


summary(fit)

#multi-dimensional regression
a = fit$coefficients[2:9]

predicts = a %*% t(mydata[,1:8])    + fit$coefficients[1]

#r = errors, residuals
r = q - predicts

source('C:/UserData/seol/Sensitivity Analyses/Lynx custom functions/multiplot.R')

library(ggplot2)
matrix_plot = cbind(t(x),q[1:1000],predicts[1:1000],r[1:1000],mydata[,1:8])


colnames(matrix_plot)[(Nsamp+1):(Nsamp+3)] = c('original','fitted','residual')
colnames(matrix_plot)[(Nsamp+4):ncol(matrix_plot)] = c('sq_y1','sq_y2','sq_y3','sq_y4','y1','y2','y3','y4')
head(matrix_plot)

matrix_plot=as.data.frame(matrix_plot)
library(gridExtra)

#plotting options
range.ylim = c(-0.5,1.5)
size.font = 20
size.font.xaxis = 20

#calculate Std of scatters
var.original = round(sd(matrix_plot$original)^2, digits=4)
var.fitted = round(sd(matrix_plot$fitted)^2, digits=4)
var.resi = round(sd(matrix_plot$residual)^2, digits=4)

#plotting custom function
plot_original_para <- function (name_para){
          if (is.numeric(name_para)==TRUE){
            name_para = colnames(matrix_plot)[name_para]
          }
  x.location.text = (max(matrix_plot[name_para])-min(matrix_plot[name_para]))*0.5+min(matrix_plot[name_para])
  
  p1 <- ggplot(matrix_plot,aes_string(name_para,'original'))+geom_point(alpha = 1/50)+ylab('RMSE')+ggtitle('f(x)')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
  annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.original)) )+
  # annotate('text', size=7, x=x.location.text, y=-0.5, label= var.original ) + 
  theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis,angle=90),
        axis.title.x = element_text(size=(size.font+10)))
    p2 <- ggplot(matrix_plot,aes_string(name_para,'fitted'))+geom_point(alpha = 1/50)+ggtitle('g(y)') +ylab('RMSE')+xlab('1')+  xlab(name_para)+ylim(range.ylim)+
      annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.fitted)) )+
  theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
        axis.title.x = element_text(size=(size.font+10)))
  p3 <- ggplot(matrix_plot,aes_string(name_para,'residual'))+geom_point(alpha = 1/50) +ggtitle('f - g')+xlab('1')+xlab(name_para)+ylim(range.ylim)+
    annotate('text', size=6, x=x.location.text, y=-0.4  , parse=T, label= paste0('sigma^2==',as.character(var.resi)) )+
  theme(text = element_text(size=size.font),axis.text = element_text(size=size.font.xaxis, angle=90),
        axis.title.x = element_text(size=(size.font+10)))

    tiff(str_c(name_para,'_RS.png'),width=800)
  grid.arrange(p1,p2,p3, ncol=3, nrow =1)
  dev.off()
  print(  grid.arrange(p1,p2,p3, ncol=3, nrow =1))
}

# plot_original_para('f')
# plot_original_para('e')

para.array = c('f','e','d','tau_q','tau_s','v_s','y1','y2','y3','y4','sq_y1','sq_y2','sq_y3','sq_y4')
for (i in para.array){
  plot_original_para(i)
}


