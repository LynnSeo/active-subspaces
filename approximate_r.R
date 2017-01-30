library(hydromad)
library(xts)
library(stringr)
library(ggplot2)
source('C:/UserData/seol/Sensitivity Analyses/Lynx custom functions/multiplot.R')
#####
#compute r = f - g with Jacobain-based AS and Hessian-based AS
#####
 
prim_dir='C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/'
sub_dir1 = 'Jacobian-based'
sub_dir2 = 'Hessian-based'

t_catchment = 'Gingera/constrained-range-1'
t_year = '70s'
no_seed = 2025
perts = 1e-06

ndir = str_c(prim_dir,sub_dir1,'/',t_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/')
ndir2 = str_c(prim_dir,sub_dir2,'/',t_catchment,'/',t_year,'/seed-',no_seed,'/',perts,'/')

j_df = read.csv(str_c(ndir,perts,'.csv'))
h_df = read.csv(str_c(ndir2,perts,'_Hessian.csv'))

# number of variables and samples
Npar = ncol(j_df)
Nsamp = nrow(j_df)

#read xy which are coordinates and output scalar
xy = read.csv(str_c(ndir,'xy.csv'))
xy=xy[seq(Nsamp),]

#active variable y = W1T dot X
y = read.csv(str_c(ndir,'y.csv'))
y = y[,seq(Nsamp)]

#residual = f - g
r = read.csv(str_c(ndir,'r.csv'))
r = as.numeric(r[seq(Nsamp)])
#output of function g which is truncated model
g = read.csv(str_c(ndir,'g.csv'))
g = as.numeric(g[seq(Nsamp)])
#output of function f which is original model
q = as.numeric(xy[,(Npar+1)])


g0 = rowMeans(g)


