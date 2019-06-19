if(require(rstudioapi)==FALSE){ install.packages('rstudioapi')}
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('Source.r')

Finaloutput=matrix(nrow=10000,ncol=21)
colnames(Finaloutput)=c('distribution.degree.of.nonnorm','truescore.degree.of.nonnorm','obs.degree.of.nonnorm','mean.true','mean.obs','var.true','var.obs','actual.fp.at.truecut','actual.fn.at.truecut',
                        'actual.tot.err.truecut','actual.opt.cut','fp.err.at.opt','fn.err.at.opt','toterr.at.opt',
                        'est.fp.at.truecut','est.fn.at.truecut','est.opt.cut','fp.est.at.opt','fn.est.atopt','toterr.est.attrue','toterr.est.atopt')


set.seed(9999)
#kurtosis

#we want to loop from k=3 to k=5.99. I proposed to do 50 total. I have (effectively) 3/50 = a deltax of .06

l=13 # this is just for the pdf() script. Used for figure number in paper
k=3.00
i=1 # i iterates with k, but needs to be integers increasing by 1

while(k<=5.999)
{           #BEGIN KURTOSIS LOOP
  
  #mixture kurtosis
  
  
var2=sqrt(((k*625)/3)-1250+625)+25
var1 = (25-.5*(var2))/.5
  
  
  x=rnorm(5000,50,sqrt(var2)) 
  y=rnorm(5000,50,sqrt(var1)) 
  
  t=c(x,y)
  

#let y be the vector of our observed scores, recall that 2.5 is our std err
z=vector(length=10000)
n=1
while(n<=10000)
{
  z[n]=rnorm(1,t[n],2.5)
  n=n+1
}

trueandobs=cbind(t,z) #true scores on left, obs scores on right


#okay, now we'll search over score ranges, taking steps of .1
findat=NULL#vector(length=4)

#we'll build this row by row, so start with a vector
c=20
while(c<=80)
{
  tableofdata=vector(length=4)
  shouldpass=subset(trueandobs,trueandobs[,1]>=40)
  
  falsenegcount=nrow(subset(shouldpass,shouldpass[,2]<c))
  
  shouldfail=subset(trueandobs,trueandobs[,1]<40)
  
  falseposcount=nrow(subset(shouldfail,shouldfail[,2]>=c))
  
  tableofdata[1]=c
  tableofdata[3]=falseposcount/10000
  tableofdata[4]=falsenegcount/10000
  tableofdata[2]=tableofdata[3]+tableofdata[4] # are you watching? this is column 2, but has to be written after 3 and 4
  
  findat=rbind(findat,tableofdata)
  
  c=c+.1
}

#need to locate point of optimal error
#and also line where truecut=55

loc55=which(round(findat[,1],digits=1)==40)

locmin=min(which(findat[,2]==min(findat[,2])))  # i take the min of the location, cause sometimes we'll get two points with same error/etc.
                                                  #in those rare instances, I'll take the left most (min on horizontal scale) location

#Kurtosis lets always use z for my mixture/convoluation
Finaloutput[i,1]=k  #the non-normality parameter (skew, kurt, or distance between modes)
Finaloutput[i,2]=kurtosis(t)
Finaloutput[i,'obs.degree.of.nonnorm']=kurtosis(z)

Finaloutput[i,'mean.true']=mean(t)
Finaloutput[i,'mean.obs']=mean(z)
Finaloutput[i,'var.true']=var(t)

Finaloutput[i,'var.obs']=var(z)

Finaloutput[i,'actual.fp.at.truecut']=findat[loc55,3]

Finaloutput[i,'actual.fn.at.truecut']=findat[loc55,4]

Finaloutput[i,'actual.tot.err.truecut']=Finaloutput[i,'actual.fn.at.truecut']+Finaloutput[i,'actual.fp.at.truecut']

Finaloutput[i,'actual.opt.cut']=findat[locmin,1]

Finaloutput[i,'fp.err.at.opt']=findat[locmin,3]

Finaloutput[i,'fn.err.at.opt']=findat[locmin,4]

Finaloutput[i,'toterr.at.opt']=Finaloutput[i,'fp.err.at.opt']+Finaloutput[i,'fn.err.at.opt']


#need to call these after the var has been written
Finaloutput[i,15:19]=gandw(Finaloutput[i,'mean.obs'],Finaloutput[i,'var.obs'])

Finaloutput[i,'toterr.est.attrue'] =Finaloutput[i,'est.fp.at.truecut']+Finaloutput[i,'est.fn.at.truecut']
Finaloutput[i,'toterr.est.atopt']= Finaloutput[i,'fp.est.at.opt']+Finaloutput[i,'fn.est.atopt']

if(k==3|k==3.54|k==4.14|k==4.74|k==5.34|k==5.94) 
{
  
  t2=as.data.frame(t)
  title=paste('Figure ', l,': True Score and Observed Score Histograms with Distribution \nTrue Score Kurtosis of',k,'\n\nTrue Score Frequencies with Kurtosis =',round(Finaloutput[i,2],digits=2))
  a=ggplot(t2,aes(t))+geom_histogram(bins=30)
  pdf(paste('40 Truescore frequencies at Kurtosis =',k,'.pdf'))
  print(a+labs(title=title))
  
  z2=as.data.frame(z)
  title=paste('Observed Score Frequencies with Kurtosis',round(Finaloutput[i,'obs.degree.of.nonnorm'],digits=2))
  a=ggplot(z2,aes(z))+geom_histogram(bins=30)
  #pdf(paste(title,'.pdf'))
  print(a+labs(title=title))
  dev.off()
  l=l+1
}

i=i+1

k=round(k+.06,digits=3) 



output=Finaloutput[1:50,]
write.csv(output,file='kurtosis40.csv')


}# END OF Kurtosis lOOP