if(require(moments)==FALSE){install.packages('moments')}
if(require(ggplot2)==FALSE){install.packages('ggplot2')}
library(ggplot2)
library(moments)

gandw=function(obsmean,obsvar)
{
  #Irina's method as a callable function for tabling
  
  #######Functions (run these first)
  phi=function(x)
  {(1/2)*(1+erf(x/sqrt(2)))}
  
  
  falseposfunc2=function(theta)  # we can no longer use error function for second part
  {
    (1-phi(((c-theta)/sderr)))*dnorm(theta,obsscoremean,truesd)#exp(-(theta^2)/2)
    
  }
  
  falsenegfunc2=function(theta)
  {
    (1-(1-phi(((c-theta)/sderr))))*dnorm(theta,obsscoremean,truesd)
  }
  
  
  
  
  #obsvar=30
  library(pracma)
  obsscoremean=obsmean
  sdobs=sqrt(obsvar)
  truecut=55
  rel=.8
  w=1
  
  
  
  
  sderr=sqrt(1-rel)*sdobs
  truevar=rel*(sdobs)^2
  truesd=sqrt(truevar)
  newmat2=matrix(nrow=100000,ncol=7)
  colnames(newmat2)=c("theta","wce","FP","FN","wcec","FPc","FNc")
  fpvec2=Vectorize(falseposfunc2,'theta')
  fnvec2=Vectorize(falsenegfunc2,'theta')
  
  pshouldfail=pnorm((truecut-obsscoremean)/truesd,0,1)
  pshouldpass=1-pshouldfail
  
  w=1
  c=20
 # end=4*sdobs+truecut
  
  i=1
  while(c<=90)
  {
    
    
    
    
    fp2=integrate(fpvec2,-Inf,truecut)#,rel.tol=1e-13)#,subdivisions)#,rel.tol=.0000000000001)
    fn2=integrate(fnvec2,truecut,Inf)#,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
    fpc=fp2$value/pshouldfail
    fnc=fn2$value/pshouldpass
    
    
    
    wce=(w*fp2$value+fn2$value)
    
    
    
    
    newmat2[i,1]=c
    newmat2[i,2]=wce
    newmat2[i,3]=fp2$value
    newmat2[i,4]=fn2$value
    newmat2[i,5]=(1/2)*(fpc+fnc)
    newmat2[i,6]=fpc
    newmat2[i,7]=fnc
    
    i=i+1
    c=c+.1
  }
  
  outmat=subset(newmat2,newmat2[,2]>0)





#optimal error point
loc=which(outmat[,2]==min(outmat[,2]))

minerrC=outmat[loc,1]

minerrFP=outmat[loc,3]
minerrFN=outmat[loc,4]


#now I need error at the true cut
a=round(as.vector(outmat[,1]),digits=1)  #for some reason we need this... error otherwise
loc=which(a[]==55.0)   #might be a rounding thing

truecutFP=outmat[loc,3]
truecutFN=outmat[loc,4]


outvector=as.matrix(c(truecutFP,truecutFN,minerrC,minerrFP,minerrFN))
rownames(outvector)=c('truecutFP','truecutFN','minerrC','minerrFP','minerrFN')
return(outvector)
}#end g&w function
