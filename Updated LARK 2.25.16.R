
#An example: the relationship between intelligence and first five mood items.
#dataset to test
setwd("C:\\Users\\Lara Mercurio\\Dropbox\\R\\LARK")
setwd("~/Dropbox/R/LARK")
setwd("C:/Users/Mark/Desktop/LARK/data")
load("Students.RData")
load("test.RData")


# LARK STEP 1 - Creating DATAFRAME and SORTED CORRELATION TABLE
# LARK takes in a dataframe where the first column is the DV and all other columns are the IVs

library(boot)
library(psych)

lark.prep <- function(data) {
  data<-na.omit(data)
  
  #dataframe of names and correlations
  df.cor<-data.frame(c(names(data[2:ncol(data)])),
                     c(abs(cor(data[1],data[2:ncol(data)]))),
                     c(cor(data[1],data[2:ncol(data)]))) 
  
  #TO RETURN: table of the sorted correlations by abs value
  df.sort.cor<-df.cor[order(df.cor[2],decreasing=T),] #includes raw and abs
  sortcor<-data.frame(df.sort.cor[1],round(df.sort.cor[3],2)) #only raw
  names(sortcor)<-c("Variable","Correlation")
  
  #Creating the final dataframe (ff) to use for bootstrapping
  ff<-data.frame(data[1])#final dataframe to use
  for (i in 1:nrow(sortcor)){
    ff[i+1] = scale(data[names(data)==as.character(sortcor[i,1])],scale=F)
    colnames(ff)[i+1] <- as.character(sortcor[i,1]) ## use original names
  }
  
  result <- list(ff=ff,sortcor=sortcor)
  return(result)
  
}

# all linear model functions can use this to get the AIC on the bootstrapped dataset
lm.aic <- function(formula, data, w) {
  d <- data[w,] #allows boot to select sample
  fit <- lm(formula, data = d)
  return(extractAIC(fit)[2])
}


                 
#Function to output the previous variable and type to enter to future steps.
#Note, we would need to change this to output "prior.var" for each
#previous predictor.

var.fit<-function(data.in,y,q=q, besttype=besttype) {
  
  #to create an initial formula to update
  if(q == 3){ 
    pre.form <- lm(data.in[,y] ~ 1)
  }
  
  if(is.null(besttype)){
  prior.var <- paste("1")
  }
    else
    if(besttype=="linear"){
    prior.var <- paste("data.in[,q-1]")
    } 
      else
      if(besttype=="quadratic"){
      prior.var <- paste("data.in[,q-1] + I(data.in[,q-1]^2)")
      }
        else
        if(besttype=="cubic"){
        prior.var <- paste("data.in[,q-1] + I(data.in[,q-1]^2) + I(data.in[,q-1]^3)")
        }
          else
          if(besttype=="quartic"){
          prior.var <- paste("data.in[,q-1] + I(data.in[,q-1]^2) + I(data.in[,q-1]^3) + I(data.in[,q-1]^4)")
          }
            else
            if(besttype=="quintic"){
            prior.var <- paste("data.in[,q-1] + I(data.in[,q-1]^2) + I(data.in[,q-1]^3 + I(data.in[,q-1]^4 + I(data.in[,q-1]^5)")
            }
  
  pre.form <- update(pre.form, . ~  .  + prior.var)
  print(pre.form)
  #The next few lines below is a second way to do it that I think is worse. 
  #In short, you can use the update function OR you can paste everything together.
  #I think the update function is more parsimonious.
  #print(data.frame(data.in[,q],data.in[,q-1]))
  #fmla.l <- as.formula(paste("data.in[,1] ~", paste(prior.var),"+ data.in[,q]"))
  #print(fmla.l)
  #print(lm(fmla.l))
  return(pre.form)
        
}

# Hijack the bootstrap method using our own function
lark.boot <- function(data.in, y, q, R = 100, seed = 42, besttype = NULL) {
  if(is.null(besttype)){
    
  names(data.in)[1]<-"y"
  names(data.in)[q]<-"iv"
  
  set.seed
  lin <- boot(data=data.in, statistic = lm.aic, R = R, formula = y ~ iv)
  set.seed(seed)
  quad <- boot(data=data.in, statistic = lm.aic, R = R, formula = y ~ iv + I(iv^2))
  set.seed(seed)
  cub <- boot(data=data.in, statistic = lm.aic, R = R, formula = y ~ iv + I(iv^2) + I(iv^3))
  set.seed(seed)
  quart <- boot(data=data.in, statistic = lm.aic, R = R, formula = y ~ iv + I(iv^2) + I(iv^3) + I(iv^4))
  set.seed(seed)
  quint <- boot(data=data.in, statistic = lm.aic, R = R, formula = y ~ iv + I(iv^2) + I(iv^3) + I(iv^4) + I(iv^5))
  }
    else {
      
  prior.var<-var.fit(data.in,y,q=q,besttype = besttype)
  #OPTION #1: Using the Update function.
  set.seed(seed)
  lin <- boot(data=data.in, statistic = lm.aic, R = R, 
              formula = update(pre.form, . ~  .  + data.in[,q]))
  set.seed(seed)
  quad <- boot(data=data.in, statistic = lm.aic, R = R, 
              formula = update(pre.var, . ~ . + (i+1) + I((i+1)^2)))
  set.seed(seed)
  cub <- boot(data=data.in, statistic = lm.aic, R = R, 
              formula = update(pre.var, . ~ .+ (i+1) + I((i+1) ^2) + I((i+1)^3)))
  set.seed(seed)
  quart <- boot(data=data.in, statistic = lm.aic, R = R, 
              formula = update(pre.var, . ~ .+ (i+1) + I((i+1) ^2) + I((i+1)^3) + I((i+1)^4)))
  set.seed(seed)
  quint <- boot(data=data.in, statistic = lm.aic, R = R, 
                formula = update(pre.var, . ~ .+ (i+1) + I((i+1) ^2) + I((i+1)^3) + I((i+1)^4) + I((i+1)^5)))
  
  
  #OPTION #2: The second way to do it that seems less parsimonious.
  #fmla.l <- as.formula(paste("data.in[,1] ~", paste(prior.var),"+ data.in[,q]"))
  #set.seed(seed)
  #lin <- boot(data=data.in, statistic = lm.aic, R = R, formula = fmla.l)
  #print(fmla.l)
  #print(lin$t)                   
                     
  
  }
  
  aic.table <- data.frame(lin$t, quad$t, cub$t, quart$t, quint$t)
  return(aic.table)
}

# Use the hijacked procedure and return the results to the user for the first variable
# We will need to pass the final model out and make this more general to make it work for vars 2-n.

lark <- function(data, y=1, x=2:ncol(data), R = 1000, seed = 42, jitter = TRUE, besttype=NULL){
  prep <- lark.prep(data.frame(data[y],data[x]))
  data.in <- prep[[1]]
  print(paste("r with", names(data.in)[1]))
  print(prep[[2]])
  print("")
  
  # repeat the whole process for every variable
  for(q in 2:ncol(data.in))  {
 
  aic.table <- lark.boot(data.in, y, q, R = R, seed = seed, besttype = besttype)
  
  
  # report % of time each relationship type was the best (lowest AIC)
  for (i in 1:nrow(aic.table)) {
    max = aic.table[i,1]
    aic.table$bestcol[i] = 1
    for (j in 2:ncol(aic.table)-1) {
      if (aic.table[i,j] < max) { 
        max = aic.table[i,j]
        aic.table$bestcol[i] = j
      }
    }
  }
  
  aic.table$bestcol <- factor(aic.table$bestcol, levels = c(1:5), labels = c(1:5))
  raw.best <- table(aic.table$bestcol)
  stats <- describe(aic.table[,1:ncol(aic.table)-1])
    
  
  low.mean = stats$mean[1]
  besttype = 1
  for (i in 2:nrow(stats)){
    if (stats$mean[i] < low.mean) {
      low.mean = stats$mean[i]
      besttype = i
    }
  }
  
  modeltypes = c("linear", "quadratic","cubic","quartic","quintic")  
  besttype = modeltypes[besttype]
  
  final.table <- matrix(cbind(100*raw.best/R, stats$mean, stats$sd), ncol = 3, byrow = FALSE, dimnames = list(modeltypes, c("% Lowest AIC", "Mean AIC", "AIC sd")))
  
  # tell the user: for this variable, % of time that each function had lowest AIC, mean AIC for each function
  # bootstrapping suggests a XXXXXX relationship
  print(paste("Bootstrapping suggests that the most likely relationship"))
  print(paste("between", names(data.in)[1], "and", 
            names(data.in)[q],"is a", besttype, "one:"))  
  
  print(signif(as.table(final.table)), digits = 2)
  
  #return(aic.table)
  
  #plot.all(data.in, y, i, jitter = jitter)
  }  
}

#To plot the fit curves of the regression models.

plot.all <- function(data.in, y, i, jitter=TRUE) {
  data.in <- lark.prep(data.in)[[1]]
  data.in$y = data.in[,y]
  data.in$x = scale(data.in[,i], scale = F)
  fit.l <- lm(y ~ x,data=data.in)
  fit.q <- lm(y ~ x + I(x^2),data=data.in)
  fit.c <- lm(y ~ x + I(x^2) + I(x^3),data=data.in)
  fit.qr <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4),data=data.in)
  fit.qi <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5),data=data.in)
  
  #set default to jitter
  if(jitter==FALSE){ 
    par(mar=c(5,2,1,5))
    plot(data.in[,"x"],data.in[,"y"],
         ylab=paste(names(data.in)[y]),xlab=paste(names(data.in)[i], "centered"))
  } else {
    par(mar=c(7,5,3,6))
    plot(jitter(data.in[,"x"]),jitter(data.in[,"y"]),
         ylab=paste(names(data.in)[y]),xlab=paste(names(data.in)[i], "centered"))
  }
  
  #creating input data for other models
  xx <- seq(min(data.in[,"x"]),max(data.in[,"x"]),.1)
  
  #plotting linear fit
  yy <- fit.l$coef %*% rbind(1,xx)
  lines(xx,yy,lwd=2,col=2,lty=1)
  
  #plotting quadratic fit
  yy <- fit.q$coef %*% rbind(1,xx,xx^2)
  lines(xx,yy,lwd=2,col=3,lty=2)
  
  #plotting cubic fit
  yy <- fit.c$coef %*% rbind(1,xx,xx^2,xx^3)
  lines(xx,yy,lwd=2,col=4,lty=3)

  #plotting quartic fit
  yy <- fit.qr$coef %*% rbind(1,xx,xx^2,xx^3,xx^4)
  lines(xx,yy,lwd=2,col=5,lty=4)
  
  #plotting quintic fit
  yy <- fit.qi$coef %*% rbind(1,xx,xx^2,xx^3,xx^4,xx^5)
  lines(xx,yy,lwd=2,col=6,lty=5)
    
  #plotting lowess line
  lines(lowess(data.in[,"x"],data.in[,"y"]), col=7,lty=6,lwd=2) # lowess line (x,y)
    
  #adding a legend
  legend("right", 0, c("linear","quadratic","cubic","quartic","quintic","lowess"), 
         col = c(2,3,4,5,6,7),lty=c(1,2,3,4,5,6),lwd=c(2,2,2,2,2,2), cex = .5, xpd = TRUE, inset = c(-.15,0))
}


# Test the final procedure using Mark's output from lark.prep
lark(data.frame(students[41],students[19:23]), 1, 2:3, R = 100, seed = 77)
lark.oneVar(test,1, 2, R = 500, seed = 77)
students[41]
# Testing the function to plot the fit curves of the regression models.
plot.all(data.frame(students[41],students[19:23]), 1, 2, jitter=T)
plot.all(test, 1, 2, jitter = T)

# Test Tanzania tb death rate vs gdp:
load("tb.death.RData")
load("gdp.tz.RData")
lark.boot(data.frame(gdp[11:31,][2], tb.death[2]), 1, 2, R = 500, seed = 42)
lark(data.frame(gdp[11:31,][2], tb.death[2]), R = 500, seed = 42)
plot.all(data.frame(gdp[11:31,][2], tb.death[2]), 1, 2, jitter = TRUE)
lark.oneVar(data.frame(tb.death[2], gdp[11:31,][2]), dv.col = 1, iv.col = 2, R = 500, 		seed = 42, jitter = FALSE)
