read.rep <- 
 function(fn)
 {
  # The following reads a report file
  # Then the 'A' object contains a list structure
  # with all the elemements in the report file.
  # In the REPORT_SECTION of the AMDB template use 
  # the following format to output objects:
  #  	report<<"object \n"<<object<<endl;
  #
  # The part in quotations becomes the list name.
  # Created By Steven Martell
  options(warn=-1)  #Suppress the NA message in the coercion to double
  
  
  ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
   ir=match(vnam[i],ifile)
   if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
   dum=NA
   if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
   if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))
   
   if(is.numeric(dum))#Logical test to ensure dealing with numbers
   {
    A[[vnam[i]]]=dum
   }
  }
  options(warn=0)
  
  return(A)
 }


ATF=read.rep("/Users/ingridspies/admbtest/2017_CIE/atf.rep") 

#compare biomass, fsb
par(mfrow=c(1,1))
par(mar=c(5,4,1,1))
plot(seq(1976,2016,1),ATF$fspbio/10000,ylim=c(0,1.2e2),type="l",ylab="Biomass (x10,000 t)",xlab="Year",cex.axis=1.4,cex.lab=1.4,lwd=3,las=1,col="red",lty=2)
lines(seq(1976,2016,1),ATF$pred_bio/10000,col="chartreuse",lwd=3)
legend("topleft",c("Model 15.0","Model 15.0a","Model 15.0b","Model 15.1","Model 15.1a","Model 15.1b"),lwd=3,col=c("black","grey","red","blue","orange","chartreuse"),cex=1.3)

#fishing mortality on fully selected fish [,8]
par(mfrow=c(1,1))
plot(seq(1976,2016,1),round(ATF$Mort_est_fem[,8],3),type="l",ylab="Full selection F",xlab="Year",cex.axis=1.4,cex.lab=1.4,lwd=2)
#lines(seq(1976,2014,1),round(ATF_15_0$Mort_est_fem[,8],3),col="red",lwd=2)
lines(seq(1976,2016,1),round(ATF$Mort_est_fem[,8],3),col="red",lwd=2)
legend("topright",c("2014 Model F", "2016 Model F"),lty=c(1,1),col=c(1,2),lwd=2,cex=1.1)

#recruitment
par(mfrow=c(1,1))
par(mar=c(5,5,2,2))
plot(c(seq(1976,2011,1)),ATF$Numbers_fem[1:36,1],type="l",ylab="Numbers (x1000)",xlab="Year",cex.axis=1.4,cex.lab=1.4,lwd=3,ylim=c(0,3.7e8))
lines(seq(1976,2009,1),ATF$Numbers_fem[1:34,1],lwd=3,col="red")
legend("bottomright",c("Model 15_0","Model 15_1b"),col=c("black","red"),lwd=4,cex=1.4)

#table of recruitments
cbind(seq(1976,2016,1,),ATF$Numbers_fem[,1],ATF$Numbers_fem[,1])

#selectivity 2016
par(mar=c(1,6,2,2))
par(mfrow=c(2,2),mar=c(5,4,2,2),cex.lab=1.2,cex.main=1.2)
plot(c(1:21),ATF$Fishsel_fem/max(ATF$Fishsel_fem),col="red",lwd=3,type="l",xlab="Age",ylab="Proportion selected",main="Fishery selectivity",ylim=c(0,1))
lines(c(1:21),ATF$Fishsel_mal/max(ATF$Fishsel_mal),col="blue",lwd=3,lty=3,xlab="Age",ylab="Proportion selected",ylim=c(0,1))
plot(c(1:21),ATF$Survsel_fem[1,]/max(ATF$Survsel_fem[1,]),col="red",lwd=3,type="l",xlab="Age",ylab="Proportion selected",main="Shelf Survey selectivity",ylim=c(0,1))
lines(c(1:21),ATF$Survsel_mal[1,]/max(ATF$Survsel_mal[1,]),col="blue",lwd=3,lty=3,xlab="Age",ylab="Proportion selected",ylim=c(0,1))
plot(c(1:21),ATF$Survsel_fem[2,]/max(ATF$Survsel_fem[2,]),col="red",lwd=3,type="l",xlab="Age",ylab="Proportion selected",main="Slope Survey selectivity",ylim=c(0,1))
lines(c(1:21),ATF$Survsel_mal[2,]/max(ATF$Survsel_mal[2,]),col="blue",lwd=3,lty=3,xlab="Age",ylab="Proportion selected",ylim=c(0,1))
plot(c(1:21),ATF$Survsel_fem[3,]/max(ATF$Survsel_fem[3,]),col="red",lwd=3,type="l",xlab="Age",ylab="Proportion selected",main="AI Survey selectivity",ylim=c(0,1))
lines(c(1:21),ATF$Survsel_mal[3,]/max(ATF$Survsel_mal[3,]),col="blue",lwd=3,lty=3,xlab="Age",ylab="Proportion selected",ylim=c(0,1))

#what do likelihoods look like?
ATF$catch_like;ATF$surv_like;ATF$sel_like;ATF$age_like_survey1;ATF$rec_like;ATF$length_like_fishery;ATF$length_like_survey
