library(lme4)
library(nlme)
library(mgcv)
library(BiodiversityR)
library(sciplot)
library(ggplot2)
library(vegan)
library(dplyr)
library(car)
library(RColorBrewer)
library(viridis)
library(Rmisc)
citation("vegan")
citation()
##README####
#this is the version that Adrian corrected that contains NMDS and community analysis
#setwd("C:\Users\spenc\OneDrive\Documents\Grad School\Research\Collected Data\Coding\smithbees")
#Cathy this version contains both my code and Adrian's. Scroll to bottom for
#NMDS and most current ANOVA code is after NMDS
#Initial ANOVAs were done by Adrian for data exploration
#first section is building the matrix from original data


# this reads the original long-format data
dat1<-read.csv(file = "smithbees.csv")
#remove 6 wasps
dat<-filter(dat1, bee.or.wasp!="wasp")

#learn about the data
summary(dat)
dat$trt<-as.factor(dat$severity)
summary(dat$trt)
summary(dat$site_id)
dat$Site<-as.factor(dat$site_id)

# lets keep only those you could categorize as either


#------------------------ building our datasets to analyze for above vs below -------------

# the reshape tool can make counts of values from long-form and turn them into short-form
library(reshape2) 

# I used the 'dcast' function to aggregate (by default it counts whatever is right of '~' in the formula) 
# so here, I counted all the 'IDs' for each 'Site' and 'trt' combination
# this gives us the identical trtxst matrix as in excel, 3 sites x 3 trts, should be 9 rows, perfect
trtxst<-dcast(dat, Site + trt  ~id)  

# sort the data by trt
attach(trtxst)

# here's where I sorted, you could change this to whatever order you want
trtxst <- trtxst[order(trt),]
detach(trtxst)

# calculate abundance (the columns changed a bit)
# ALC adjusted for removed wasp IDs
trtxst$Abund<- rowSums(trtxst[,3:48])

# lets check to see that it adds up to 448 - 6 wasps should bbe 442
# yep!
sum(trtxst$Abund)

# just to illustrate, lets plot by site, year and nest

ggplot(data = trtxst, aes(x=trt, y=Abund)) + 
  geom_boxplot()

# and species richness
trtxst$S<-specnumber(trtxst[,3:48])
summary(trtxst$S)

ggplot(data = trtxst, aes(x=trt, y=S)) + 
  geom_boxplot() +
  theme(legend.position = "right")

# you can add a column for Shannon's diversity index H' that references the species columns
trtxst$H<-diversity(trtxst[,3:48])

summary(trtxst$H)
##Cathy this is rarefaction analysis, not sure if I will use for my paper but
#Adrian (ALC) ran it for me 

# ALC: for rarefaction, you use the lowest common abundance
summary(trtxst$Abund)

# use 2 for rarefaction, it'll subsample 84 bees in each site
summary(trtxst$Abund)

# so the minimun is 13 bees
trtxst$Srare<-rarefy(trtxst[,3:48],13)

# you can see how they're related
par(mfrow=c(1,1))
plot(trtxst$S,trtxst$Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

# and plot the individual based rarefaction curves for each row (site)
#Alaina-got confused at this point, should I change row numbers for HB, LB, UB??
# first let's breakout by year and then make locality the row names (can't have duplicate)

HB <- trtxst[1:3,]
LB <- trtxst[4:6,]
UB <- trtxst[7:9,]

# with these, you could use extrapolated species richness for each site as a response, Chao in this case, extrapolated from the abundance data
# I had to make a transposed matrix (t()) of the chao estimate using th estimateR function
chao<-t(estimateR(trtxst[,3:48],index="chao"))

# You can double check that the observed species rich lines up with the rows from sitexyrxn, yep
# Then just used the second column (Chao richness estimator) as a new column
trtxst$Chao<-chao[,2]
summary(trtxst$Chao)

# Now you can plot all three by year if you want, I used melt to add the columns as factors for plotting
library(reshape2)

fig<- melt(trtxst,id.vars='trt', measure.vars=c('S','Srare','Chao'))

## breaking it out by metric
# ALC: now looks better for rare richness with minimum 13 instead of 2
ggplot(data = fig, aes(x=trt, y=value, fill=variable)) +
  geom_boxplot() +
  theme(legend.position = "right")

# or barcharts
# followed by anovas
par(mfrow=c(1,1))

attach(trtxst)
# ALC: something was screwey here
trtxst$trt1 <- factor(trtxst$trt, levels=c('UB', 'LB', 'HB'))

# ALC: I just use hist to know how to set y-lims
hist(Abund)
####Cathy ignore these ANOVAs, these are Adrian's. I'll include my most
###current code for two way ANOVA w/interaction at the bottom, after the NMDS
# this tests for normality
shapiro.test(Abund)
par(mar=c(5,5,1,1))
bargraph.CI(trt,Abund, xlab="Burn severity", ylab="Flower visitor abundance",ylim=c(0,60),
            cex.axis=1.25, cex.lab=1.75,legend=T)
abline(h=0)
mod<-(lm(Abund~trt))
summary(mod)
a.av <-aov(mod)
summary(a.av)
ab<- ggplot(data = trtxst, aes(x=trt, y=S, fill=trt)) + geom_boxplot(lwd=2, fatten=2) + theme(legend.position = "right") + ylim(0,25)
ab + scale_fill_manual("Fire Severity",values=c("#915E50","#B8AB83","#037171")) + labs( x = "Burn Severity", y = "Flower Visitor Abundance") + theme(text =element_text(size=48, family = "serif"),panel.background=element_blank(), axis.line = element_line(color="black"))
# NS
# note if you run a full model with more variables these would change
text(3,13,"p = 0.09")


########### Adrian's species richness analysis ###
hist(S)
shapiro.test(S)
par(mar=c(5,5,1,1))
bargraph.CI(trt,S, xlab="Burn severity", ylab="Flower visitor richness",
            cex.axis=1.25, cex.lab=1.75,ylim=c(0,20))
abline(h=0)
sr <- ggplot(data = trtxst, aes(x=trt, y=S, fill=trt)) + geom_boxplot(lwd=2, fatten=2) + theme(legend.position = "right") + ylim(0,25)
sr + scale_fill_manual("Fire Severity",values=c("#915E50","#B8AB83","#037171")) + labs( x = "Burn Severity", y = "Flower Visitor Richness") + theme(text =element_text(size=48, family = "serif"),panel.background=element_blank(), axis.line = element_line(color="black"))
windowsFonts()
#Alaina changed to log since non normal
mod<-(lm(log(S)~trt))
summary(mod)
 #Alaina ran TukeyHSD
s.av <-aov(mod)
 summary(s.av)
 TukeyHSD(s.av)
# NS
# ALC: this one is significant p: 0.03263
# text(3,37,"p = 0.13")
text(3,18,"p = 0.033")


############Adrian's Chao or extrapolated richness####
hist(Chao)
shapiro.test(Chao)
# sig diff from normal so log transform
shapiro.test(log(Chao))
par(mar=c(5,5,1,1))
bargraph.CI(trt,Chao, xlab="Burn severity", ylab="Extrapolated richness (Chao2)",
            cex.axis=1.25, cex.lab=1.75, ylim=c(0,40))
abline(h=0)
mod<-(lm(log(Chao)~trt))
summary(mod)
#NS
text(3,56,"p = 0.1246")

# Here's some community analayses
# remember all this data is from a small sample size so while not all significant, interesting patterns exist...

######### Species accumulation curve #########
##Cathy, I'm unsure if I will use this either but Adrian ran it
## these are the opposite of rarefaction, and are hwat estimate etrapolated sp richness Chao

## first, make an object of the curves for the columns of species abundances across all sites
spec<-specaccum(trtxst[,3:48],"random")
spec

# then can estimate the # extrapolated 'species'
specpool(trtxst[,3:48],smallsample=T)

# here, you sampled 46 species across all sites, and it's estimating regional richness as 55.97 =/- 7.3 species
boxplot(spec, col="gray25",cex.axis=1, cex.lab=1.5,lwd=1)
text(15,5,"chao: 56", cex=1.5)

46/56
# so you sampled roughly 82% of the estimated regionalspecies pool, not bad


# lets plot instead of boxplot and change y axis to add in the chao line
par(mar=c(5,5,3,3))
plot(spec, xlim =c(0,18), ylim=c(0,105))

#make a hashed line representing Chao
abline(h=62, lty=2)

# Make a model (here it's an asymptotic model) to fit the data, check the differnt models: 
# I was reading here: https://urldefense.com/v3/__https://rdrr.io/rforge/vegan/man/specaccum.html__;!!NCZxaNi9jForCP_SxBKJCA!SKQkyKvNUZdwREP5bQqeq9DsxI7xAfe97u87pdrlhXGO9ruw0v5wVVTx73o2-M9LWZh8w1XxeDwTM4iDNZty1XTBzqJndZpR$ 
# The permissible alternatives are "arrhenius" (SSarrhenius), "gleason" (SSgleason), "gitay" (SSgitay), 
# "lomolino" (SSlomolino) of vegan package. In addition the following standard R models are available: 
# "asymp" (SSasymp), "gompertz" (SSgompertz), "michaelis-menten") (SSmicmen), "logis" (SSlogis), "weibull" (SSweibull).

mod<-fitspecaccum(spec, "asymp")
plot(mod)

# create a new range of data to predict from the asymptotic model (1-9 every 1), calculate predicted values over that range, and plot it
# ALC: FYI, you could plot the different burns as different curves
newdata=data.frame(x=seq(1,9,1))
pred<-predict(mod, newdata)
plot(pred, xlim=c(0,18), cex=.5,ylim=c(0,105), xlab="sites", ylab="number of species")
plot(spec,ci.type = "polygon", ci.col = rgb(0,0,1,alpha=0.4), col="black",add=T)
abline(h=56,lty=2)
text(1,80,"Chao: 56", cex=1)


########## NMDS (runs randomly so try again if it doesnt converge) #########
# ALC: want converging with 2 dimentions (K) so upped to 3
nmds1<-metaMDS(trtxst[,3:48], k=3,distance="bray")
plot(nmds1,disp="sites")
# report the dimensions (3), the stress (.153...)
nmds1

#make different symbols and colors for factors, important to note the classes (HB,LB, UB, etc.) are read in alphabetical order
# I just did it manually for each line since you just have 9
#Alaina- this didn't work for me
pch1<-c(18,18,18,19,19,19,17,17,17)
col1<-c("#882255","#882255","#882255","#DDCC77","#DDCC77","#DDCC77","#44AA99","#44AA99","#44AA99")

#and a frame of the type names
fr.nm<-c("high","low", "unburned")

# then plot the ordination with different symbols for sites, first the plot with 'null data'
par(mar=c(5,5,1,1))
plot(nmds1,disp="sites", cex.axis=1.0, cex.lab=1.5,lwd=3)
ordiellipse(nmds1, kind="se",conf=0.95,trtxst$trt, border=c("black"),lwd=1,display = "sites", draw=c("polygon"), col=c("gray95","gray60","gray3"), alpha=180,label=F)
# note how much overlap there is?

#let's redo the axis to get all the ellipses
plot(nmds1, disp="sites",type="n",cex.axis=1.0, cex.lab=1.5,lwd=3)

#then add the 95% confidence elipses around grouped sites for fire severity
ordiellipse(nmds1, kind="se",conf=0.95,trtxst$trt, border=c("black"),lwd=1,display = "sites", draw=c("polygon"), col=c("#882255","#DDCC77","#44AA99"), alpha=180,label=T)
legend("topright", legend=c("Unburned", "Low","High"), fill=c("#44AA99","#DDCC77","#882255"))
#add points to the plot
points(nmds1, display = "sites",lwd=1, pch = pch1, cex=2.0, col=col1)

#add the stress and dimensions in from the ordination
text(1,-.4, labels="stress = 0.15",cex=.75)
text(1,-.5, labels="k = 3",cex=.75)

#you can do an ANOVA of centroids of the 95% ellipses (permutational analysis of variance of tpye centroids, not significant (p= 0.274)
adonis(trtxst[,3:48]~trtxst$trt)

text(1,-.6, labels="PERMANOVA, P = 0.195",cex=0.75)

#####Most Current ANOVA Code#####
# ALC: I just use hist to know how to set y-lims
hist(Abund)

# this tests for normality
shapiro.test(Abund)
#Bee abundance is normal 
par(mar=c(5,5,1,1))
bargraph.CI(trt,Abund, xlab="Burn severity", ylab="Flower visitor abundance",ylim=c(0,60),
            cex.axis=1.25, cex.lab=1.75,legend=T)
abline(h=0)
#two way ANOVA with interaction effect abundance
aovA<-aov(Abund~trt*fire)
summary(aovA)
#this matches JMP, not significant
#Only graph severity since no effect of fire identity or interaction 
#publication quality plot of bee abundance by severity
#need to change the colors to colorblind friendly
plot.abun<- ggplot(data = trtxst, aes(x=trt, y=Abund, fill=trt)) + geom_boxplot(lwd=0.5) + theme(legend.position = "top", legend.text=element_text(size=12)) + coord_cartesian(ylim=c(0,50)) + scale_y_continuous(expand= c(0,0), limits=c(0,50)) +
  scale_fill_manual("Burn Severity",values=c("#920000","#ffff6d","#24ff24")) + labs( x = "Burn Severity", y = "Flower Visitor Abundance") + theme(panel.background=element_blank(), axis.line = element_line(color="black")) +
  theme_cowplot(font_size=12, line_size=0.5)+theme(legend.position = "right", legend.text = element_text(size=12)) + theme(text=element_text(size=12, family="Arial")) +
  annotate("text", x='high', y=49, label="a", size=4) + 
  annotate("text", x='low', y=34, label="a", size=4) +
  annotate("text", x='unburned', y=41, label="a", size=4)
plot.abun
#loadpackage for gridding plots into one figure
require(gridExtra)
windowsFonts()
# NS




########### Species Richness ###########
#test species richness for normality
hist(S)
shapiro.test(S)
#species richness is not normal
lS<-log10(S)
shapiro.test(lS)
hist(lS)
#log transformed species richness is normal
par(mar=c(5,5,1,1))
bargraph.CI(trt,S, xlab="Burn severity", ylab="Flower visitor richness",
            cex.axis=1.25, cex.lab=1.75,ylim=c(0,20))
abline(h=0)
#Publication Quality Graph of Species Richness, need to change to colorblind friendly
plot.logsr<- ggplot(data = trtxst, aes(x=trt, y=lS, fill=trt)) + geom_boxplot(lwd=0.5) + theme(legend.position = "top", legend.text=element_text(size=12)) + coord_cartesian(ylim=c(0,1.5)) + scale_y_continuous(expand= c(0,0), limits=c(0,1.5), breaks=seq(0,1.5,0.5)) +
  scale_fill_manual("Burn Severity",values=c("#920000","#ffff6d","#24ff24")) + labs( x = "Burn Severity", y = "log(Flower Visitor Richness)") + theme(panel.background=element_blank(), axis.line = element_line(color="black")) +
  theme_cowplot(font_size=12, line_size=0.5)+theme(legend.position = "right", legend.text = element_text(size=12)) + theme(text=element_text(size=12, family="Arial")) +
  annotate("text", x='high', y=1.40, label="a", size=4) + 
  annotate("text", x='low', y=1.10, label="b", size=4) +
  annotate("text", x='unburned', y=1.17, label="ab", size=4)
plot.logsr
plot.sr<- ggplot(data = trtxst, aes(x=trt, y=S, fill=trt)) + geom_boxplot(lwd=0.5) + theme(legend.position = "top", legend.text=element_text(size=12)) + coord_cartesian(ylim=c(0,25)) + scale_y_continuous(expand= c(0,0), limits=c(0,25))+
  scale_fill_manual("Burn Severity",values=c("#920000","#ffff6d","#24ff24")) + labs( x = "Burn Severity", y = "Flower Visitor Richness") + theme(panel.background=element_blank(), axis.line = element_line(color="black")) +
  theme_cowplot(font_size=12, line_size=0.5)+theme(legend.position = "right", legend.text = element_text(size=12)) + theme(text=element_text(size=12, family="Arial")) +
  annotate("text", x='high', y=23, label="a", size=4) + 
  annotate("text", x='low', y=12, label="b", size=4) +
  annotate("text", x='unburned', y=14, label="ab", size=4)
plot.sr
#Alaina changed to log since non normal
#Two way ANOVA log(S) by fire and severity
aovS<-aov(lS~trt*fire)
summary(aovS)
#significant effects of treatment (severity)
#Alaina ran TukeyHSD on treatment only since only significant variable
TukeyHSD(aovS, 'trt')
#P value is 0.0239, low and high are significantly different


###Species Diversity (H)####
hist(H)
shapiro.test(H)
#species diversity is normal
aovH<-aov(H~trt*fire)
summary(aovH)
#no factors are significant for diversity
#graph diversity
Hplot <- ggplot(data = trtxst, aes(x=trt, y=H, fill=trt)) + geom_boxplot(lwd=0.5) + theme(legend.position = "top", legend.text=element_text(size=12)) + coord_cartesian(ylim=c(0,3.25)) + scale_y_continuous(expand= c(0,0), limits=c(0,3.25), breaks=seq(0,3.25,0.75))+
  scale_fill_manual("Burn Severity",values=c("#920000","#ffff6d","#24ff24")) + labs( x = "Burn Severity", y = "Flower Visitor Diversity") + theme(panel.background=element_blank(), axis.line = element_line(color="black")) +
  theme_cowplot(font_size=12, line_size=0.5)+theme(legend.position = "right", legend.text = element_text(size=12)) + theme(text=element_text(size=12, family="Arial")) +
  annotate("text", x='high', y=3.05, label="a", size=4) + 
  annotate("text", x='low', y=2.4, label="a", size=4) +
  annotate("text", x='unburned', y=2.3, label="a", size=4)
Hplot
###Create Figure with all Plots####
library(ggpubr)
Figure1<-ggarrange(plot.abun, plot.logsr,Hplot,legend="none")
Figure1
###Create Figure1####

###Works Now, makes pdf plot##
ggsave("Figure1.pdf", width=10, height =7, unit= 'in', dpi=300, plot=Figure1, device=cairo_pdf)
print(Figure1)
dev.off()

 