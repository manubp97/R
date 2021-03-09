

source("order-data-frame-lib.R")

#################################################################
FILE <- "genbase" 
FILE <- "medical" 
FILE <- "yeast"

LOOP.iteractions <- 50
NVARS <- 8
CLASS <- 4
SEED <- 2028
set.seed( SEED)
df <- read.table( paste( FILE, ".csv", sep=""), sep=",", header=TRUE, dec=".")
names( df);#summary( df);
dim(df)

range.class.medical <- 1450:dim(df)[2]
range.preds.medical <- 1:1449
range.class.yeast <- 104:dim(df)[2]
range.preds.yeast <- 1:103
range.class.genbase <- 1187:dim(df)[2]
range.preds.genbase <- 2:1186

range.class <- range.class.genbase 
range.class <- range.class.medical 
range.class <- range.class.yeast 

range.preds <- range.preds.genbase 
range.preds <- range.preds.medical 
range.preds <- range.preds.yeast 

cat("order-data-frame-opt2\n",
  "FILE", FILE," NVARS: ",NVARS," CLASS",CLASS," LOOP: ",LOOP.iteractions," SEED: ",SEED,"\n")

#################################################################
df16r <- remove.const.vars( df[,range.preds])  ##
print( dim(df16r))
df16 <- df16r[,sample(1:dim(df16r)[2],NVARS)]
print( names( df16));#summary( df16);
print( dim(df16))
print(summary(df16))
CLASS <- min(CLASS,dim(df16)[1])
CLUSTERING <- pam( df[,range.class], k=CLASS)$clustering

df16 <- discretize( df16)
df16 <- cbind( df16, class=CLUSTERING)
print(summary(df16))
df16 <- order.index.data.frame( df=df16)  ## ~ KBMR
## set unKB items to fill KBMR (tiny index ~ 20 binary vars; up to 20 ~ KBMR!)
fitness( df=df16)
print( names( df16));#summary( df16);
print( dim(df16))

## optimization by chance permutations
df.OPT <- NULL
IT <- dim(df16)[1] ## list lenght
IT. <- 0
if (0) {
###
for( p in 1:LOOP.iteractions) {
  cat("LOOP: ",p,"\n")
  df16 <- index.permutation( df=df16)
  df16 <- order.index.data.frame( df=df16)  ## ~ KBMR
  it <- fitness( df=df16)
  if ( IT > it) {
    df.OPT <- df16 
    IT <- it
    } else if ( IT. < it) IT. <- it
}

write.table(df.OPT,file=paste( FILE, NVARS, "-opt", IT, ".csv", sep=""),sep=";",row.names=F)

cat("opt IT: ",IT,"   bad IT.:",IT.,"\norder-data-frame-opt2\n")

} else {
  cat("LOOP: ",p,"\n")

set.seed(12364)
fn.call <<- 0
out <- GenSA(par=NULL, lower = rep(1,NVARS), upper = rep(NVARS,NVARS), fn = FitnessGenSA, df16,
             control=list(threshold.stop=NULL, verbose=TRUE, max.call=3000, smooth=FALSE, max.time=180))
fn.call.GenSA = fn.call             
out[c("value","par","counts")]

par.ix <- order( out$par)
print( par.ix)
df16 <- df16[ , c(par.ix,dim( df16)[2])] ## base exchange, permutacion
df.OPT <- order.index.data.frame( df=df16)  ## ~ KBMR

write.table(df.OPT,file=paste( FILE, NVARS, "-opt", IT, ".csv", sep=""),sep=";",row.names=F)

cat("opt IT: ",IT,"   bad IT.:",IT.,"\norder-data-frame-opt2\n")

}

###
#out <- GenSA(lower = lower, upper = upper, fn = scoring,
#             control=list(threshold.stop=global.min+tol,verbose=TRUE))
#out[c("value","par","counts")]




