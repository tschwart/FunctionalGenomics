#Set your working directory
setwd("~/functional genomics/group/trimmed counts/BWA_Counts")


#read in dataset
depth1in=read.table(file="1_In_in_Counts.txt",header=F)
depth1=read.table(file="1_CGATGT_Counts.txt",header=F)
depth2=read.table(file="2_TGACCA_Counts.txt",header=F)
depth3=read.table(file="3_GCCAAT_Counts.txt",header=F)
depth4=read.table(file="4_CGATGT_Counts.txt",header=F)
depth5=read.table(file="5_TGACCA_Counts.txt",header=F)
depth6=read.table(file="6_GCCAAT_Counts.txt",header=F)
depth7=read.table(file="7_CTTGTA_Counts.txt",header=F)
depth8=read.table(file="8_CGATGT_Counts.txt",header=F)
depth9=read.table(file="9_TGACCA_Counts.txt",header=F)
depth10=read.table(file="10_GCCAAT_Counts.txt",header=F)
depth10in=read.table(file="10-In_CTTGTA_Counts.txt",header=F)








#1-IN
#initial plot
#plot(depth$V2,depth$V3,main="chrXL_group3b Coverage",xlab="Position",ylab="Coverage",type="l")
#too much noise, need to smooth the data

#also want to make x-axis nicer by converting positions to Megabases
mb1in=depth1in$V2/10^6

#perform loess smoothing
lo1in=loess(depth1in$V3~depth1in$V2,span=0.01,data.frame(x=depth1in$V2,y=depth1in$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="1_In_in.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb1in,depth1in$V3,main="1_In_CTTGTA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb1in,predict(lo1in),lwd=3,col="red")
dev.off()
#---------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb2=depth2$V2/10^6

#perform loess smoothing
lo2=loess(depth2$V3~depth2$V2,span=0.01,data.frame(x=depth2$V2,y=depth2$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="2_TGACCA.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb2,depth2$V3,main="2_TGACCA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb2,predict(lo2),lwd=3,col="red")
dev.off()
#----------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb1=depth1$V2/10^6

#perform loess smoothing
lo1=loess(depth1$V3~depth1$V2,span=0.01,data.frame(x=depth1$V2,y=depth1$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="1_CGATGT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb1,depth1$V3,main="1_CGATGT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb1,predict(lo1),lwd=3,col="red")
dev.off()
#------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb3=depth3$V2/10^6

#perform loess smoothing
lo3=loess(depth3$V3~depth3$V2,span=0.01,data.frame(x=depth3$V2,y=depth3$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="3_GCCAAT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb3,depth$V3,main="3_GCCAAT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb3,predict(lo3),lwd=3,col="red")
dev.off()
#--------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb4=depth4$V2/10^6

#perform loess smoothing
lo4=loess(depth4$V3~depth4$V2,span=0.01,data.frame(x=depth4$V2,y=depth4$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="4_CGATGT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb4,depth4$V3,main="4_CGATGT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb4,predict(lo4),lwd=3,col="red")
dev.off()
#-------------------------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb5=depth5$V2/10^6

#perform loess smoothing
lo5=loess(depth5$V3~depth5$V2,span=0.01,data.frame(x=depth5$V2,y=depth5$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="5_TGACCA.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb5,depth5$V3,main="5_TGACCA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb5,predict(lo5),lwd=3,col="red")
dev.off()
#-------------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb6=depth6$V2/10^6

#perform loess smoothing
lo6=loess(depth6$V3~depth6$V2,span=0.01,data.frame(x=depth6$V2,y=depth6$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="6_GCCAAT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb6,depth6$V3,main="6_GCCAAT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb6,predict(lo6),lwd=3,col="red")
dev.off()
#---------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb7=depth7$V2/10^6

#perform loess smoothing
lo7=loess(depth7$V3~depth7$V2,span=0.01,data.frame(x=depth7$V2,y=depth7$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="7_CTTGTA.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb7,depth7$V3,main="7_CTTGTA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb7,predict(lo7),lwd=3,col="red")
dev.off()
#------------------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb8=depth8$V2/10^6

#perform loess smoothing
lo8=loess(depth8$V3~depth8$V2,span=0.01,data.frame(x=depth8$V2,y=depth8$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="8_CGATGT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb8,depth8$V3,main="8_CGATGT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb8,predict(lo8),lwd=3,col="red")
dev.off()
#-----------------------------------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb9=depth9$V2/10^6

#perform loess smoothing
lo9=loess(depth9$V3~depth9$V2,span=0.01,data.frame(x=depth9$V2,y=depth9$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="9_TGACCA.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb9,depth9$V3,main="9_TGACCA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb9,predict(lo9),lwd=3,col="red")
dev.off()
#-------------
#also want to make x-axis nicer by converting positions to Megabases
mb10=depth10$V2/10^6

#perform loess smoothing
lo10=loess(depth10$V3~depth10$V2,span=0.01,data.frame(x=depth10$V2,y=depth10$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="10_GCCAAT.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb10,depth10$V3,main="10_GCCAAT Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb10,predict(lo10),lwd=3,col="red")
dev.off()
#--------------------------------------------------------------------------------------------
#also want to make x-axis nicer by converting positions to Megabases
mb10in=depth10in$V2/10^6

#perform loess smoothing
lo10in=loess(depth10in$V3~depth10in$V2,span=0.01,data.frame(x=depth10in$V2,y=depth10in$V3))
#lo=loess(depth$V3~depth$V2,enp.target=tail(depth$V2,1)/1000,data.frame(x=depth$V2,y=depth$V3))
#function will take ~1-2 minutes to run depending on your OS

#make blank plot
#plot(mb,depth$V3,main="1_In_in Coverage",xlab="Position (Mb)",ylab="Coverage",type="n")

#add loess prediction
#lines(mb,predict(lo),lwd=3,col="red")

#can also plot both raw and smoothed data together
png(file="10-In_CTTGTA.png", res = 300, width = 90, height = 100, units = "mm")
plot(mb10in,depth10in$V3,main="10-In_CTTGTA Coverage",xlab="Position (Mb)",ylab="Coverage",type="l",col="grey50")
lines(mb10in,predict(lo10in),lwd=3,col="red")
dev.off()



