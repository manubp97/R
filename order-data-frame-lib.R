 

#if (!require("RWeka"))
#install.packages("RWeka") ## arff format
#library(RWeka)
if (!require("cluster"))
install.packages("cluster") ## pam cluster
library(cluster)
if (!require("dplyr"))
install.packages("dplyr")
library(dplyr)
if (!require("GenSA"))
install.packages("GenSA")
library(GenSA)
if (!require("mlr"))
install.packages("mlr")
library(mlr)
if (!require("infotheo"))
install.packages("infotheo")
library(infotheo)

###
discretize <- function( df.X, levels=2) {
print(names(df.X))
print(dim(df.X))
for( i in 1:(dim(df.X)[2]-1)) {
  brks <- quantile( df.X[ ,i], probs = seq(0.0, 1.0, 0.25))
  print(brks)
  for( b in 2:5)
    ##df.X[ df.X[ ,i] < brks[b],i] <- ( brks[b]-brks[b-1])/2 
    df.X[ df.X[ ,i] < brks[b],i] <- b 
}
print(dim(df.X))
  return( df.X)
}

## class items count ~ fitness
##  class variable is at last column 
fitness <- function( df.X) {  
it <- 1
class <- df.X[1,dim(df.X)[2]] ## first row class
for( i in 2:dim(df.X)[1])
  if ( class != df.X[i,dim(df.X)[2]]) { ## change test
    it <- it+1
    class <- df.X[i,dim(df.X)[2]] ## change!
  }
  cat("it num: ", it,"\n")
  return( it)
}

## class items count ~ fitness
##  class variable is at last column 
FitnessGenSA <- function( par, df.X) {  
fn.call <<- fn.call + 1
par.ix <- order( par[-dim(df.X)[2]]) ## == permut 1:dim(df.X)[2]
df.X <- df.X[ , c(par.ix,dim( df.X)[2])] ## base exchange, permutacion
df.X <- order.index.data.frame( df=df.X)  ## ~ KBMR
it <- 1
class <- df.X[1,dim(df.X)[2]] ## first row class
for( i in 2:dim(df.X)[1])
  if ( class != df.X[i,dim(df.X)[2]]) { ## change test
    it <- it+1
    class <- df.X[i,dim(df.X)[2]] ## change!
  }
  cat("it num*: ", it,"\n")
  print(par)
  print(par.ix)
  if ( IT > it) {
    df.OPT <<- df16 
    IT <<- it
    } else if ( IT. < it) IT. <<- it
  
  return( it)
}

##
remove.const.vars <- function( df.X) {
vcte <- c() ## remove const variables
for( j in 1:dim( df.X)[2]) {
  df.X[ , j] <- as.numeric( df.X[ , j])-1
  if ( min( df.X[ , j]) == max( df.X[ , j])) vcte <- c( vcte, j)
  }
if ( length( vcte) > 0) {
  if ( length( vcte) < dim(df.X)[2])  
    df.X <- df.X[,-vcte]  
    else stop("No vars!")
    }
return( df.X)
}

## permutacion colunms with invariant class variable vector order
index.permutation <- function( df.X) {
  df.X <- df.X[ , c( sample( 1:( dim( df.X)[2]-1)), dim( df.X)[2])] ## index permutation
  return( df.X)
}

##  order by sorted vars from left to rigth and proper class variable permutacion
order.index.data.frame <- function( df.X) {

if (0) { ## ?
cols.list <- list( df.X[ , 1])
for( i in 2:( dim(df.X)[2]-1))
  cols.list <- c( cols.list, list( df.X[ , i]))
#print( str( cols.list))
df.X <- df.X[ do.call( order, ( cols.list)),] ## all columns
return( df.X)
} else {
  df.X <- df.X[ order( df.X[ , 1]),] ## first column
  cat("i: \n")
  for( i in 2:( dim(df.X)[2]-1)) { ## columns
    cat(" ",i)
    ###   gc(); mask <- df.X[,c(1:i-1)]; ## ordered variables 
    j<-1;
    while( j < dim( df.X)[1]) { ## blocks
      for( z in ( j+1):dim( df.X)[1]) {
        df.jz.rows <- as.data.frame( df.X[ c(j,z), c(1:i-1)]); 
        if ( sum( df.jz.rows[ 1,] != df.jz.rows[ 2,]) >= 1) {
          k<-z-1;break; } else if ( z == dim( df.X)[1]) k<-z;
        }
      if ( (k-j) > 0 & ( min( df.X[ j:k,i]) < max( df.X[ j:k,i]))) {
        df.X[ j:k, i:dim( df.X)[2]] <- df.X[ order( df.X[ j:k,i]) + j - 1, i:dim( df.X)[2]]; ## block
        cat(".")
        }
      j<-z;
      }    
    }
  cat("\n")
  return( df.X)
  }
}

## feature selection from mlr
feature.selection <- function(df, thr){
fv = generateFilterValuesData(df, method = "FSelectorRcpp_information.gain")
filtered.task = filterFeatures(df, fval = fv, threshold = thr)
df <- filtered.task
return(df)
}

## conditional mutual information computation
cond.info <- function(df){
for (i in 1:(dim(df)[2]-2))
  for (j in (i+1):(df(dim[2]-1)))
    
matrix(dim(df)[2]-1:dim(df[2]-1)) <- condinformation(df[,i], df[,j], df[,dim(df)[2]], method="emp")
}