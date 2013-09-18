setwd("~/Documents/git/multivariate")
library("ggplot2")
library("reshape")
source("test.funcs.R")
source("sim.funcs.R")

#extract subset of vector with particular name
v.subset = function(x,n){
  return(x[names(x)==n])
}

#examples for ASHG poster/paper

rho = 0.7
beta1 = 1.6*c(0,-0.1,0.1/rho) #1.6 just chosen to give nice separation
beta2 = 1.6*c(0.2,0.2,0.1)
f=0.5
set.seed(150294)


#done to produce biplots for paper (big effect size to give nice separation)
res=data.frame()
mm=5
for(j in 1:length(rho)){
  for(i in 1:length(beta1)){
s = simset(500,c(mm*beta1[i],mm*beta2[i]),rho[j],f,nrep=1)
d = data.frame(beta1[i],beta2[i],rho[j],cbind(t(s$s[[1]]$Y),s$s[[1]]$X,i))
res = rbind(res,d)
}
}
names(res)=c("beta1","beta2","rho","Height","Weight","X","simtype")

res$Genotype = factor(res$X)
res$scenario = c("a","b","c")[res$simtype]

pdf("biplot_paperexamples.pdf",height=3,width=7)
#png("biplot_ashgexamples.png",height=240,width=960)
p=ggplot(data=res,aes(Height,Weight,colour=Genotype,group=1))+geom_point(shape=16)+ facet_grid(.~scenario)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()


#examples, with beta1 = (-0.1,0,.2/rho[2],.2/rho[3]) beta2 = 0.2
mm=5
rho = c(0,0.4,0.7,0.9)
beta1 = c(-0.1,0,0.2*rho[2],0.2*rho[3])
beta2 = c(0.2,0.2,0.2,0.2)
f=0.2
set.seed(150294)

res1=data.frame()
for(j in 1:length(rho)){
  for(i in 1:length(beta1)){
s = simset(1000,c(mm*beta1[i],mm*beta2[i]),rho[j],f,nrep=1)
d = data.frame(round(beta1[i],2),round(beta2[i],2),rho[j],cbind(t(s$s[[1]]$Y),s$s[[1]]$X))
res1 = rbind(res1,d)
}
}
names(res1)=c("beta1","beta2","rho","Y1","Y2","Genotype")
res1$Genotype = as.factor(res1$Genotype)
#res1$Genotype = factor(res1$X)

#res$mycolor= factor(res$X,labels=c("red","green","blue"))

png("biplot_examples_subset_rho04_3beta.png",height=240,width=960)
p=ggplot(data=subset(res1,rho==0.4 & beta1<.1),aes(Y1,Y2,colour=Genotype))+
  geom_point(alpha=1,size=1)+
  facet_grid(rho~beta1)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()


png("biplot_examples_subset_rho07.png",height=240,width=960)
p=ggplot(data=subset(res1,rho==0.7),aes(Y1,Y2,colour=Genotype,group=1))+geom_point(aes(alpha=1),size=1,shape=1) + facet_grid(rho~beta1)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()


png("biplot_examples.png",height=960,width=960)
p=ggplot(data=res1,aes(Y1,Y2,colour=Genotype,group=1))+geom_point(aes(alpha=1),size=1,shape=1) + facet_grid(rho~beta1)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()

res = data.frame()
for(j in 1:length(rho)){
  for(i in 1:length(beta1)){
    s = simset(1000,c(beta1[i],beta2[i]),rho[j],f,nrep=50)
    d = data.frame(lbfall=v.subset(unlist(s$t.01),"lbfall"),lbfuni=v.subset(unlist(s$t.01),"lbfuni"),lbfav = v.subset(unlist(s$t.01),"lbfav"),lbfuni.1=v.subset(unlist(s$t.01),"lbfuni.1"), lbfuni.2=v.subset(unlist(s$t.01),"lbfuni.2"), lbf2given1=v.subset(unlist(s$t.01),"lbf2given1"), lbf1given2=v.subset(unlist(s$t.01),"lbf1given2"),rho=rho[j],beta1=round(beta1[i],2),beta2=round(beta2[i],2))
    res = rbind(res,d)
  }
}

res$ref = res$lbfav
res$ref = ifelse(res$beta1==-0.1, res$lbfall,res$ref)
res$ref = ifelse(res$beta1==0, res$lbf2given1,res$ref)
res$ref = ifelse(res$beta1>0.01,res$lbfuni.2,res$ref)
#names(res)[1] = "lbfall"
#names(res)[2] = "lbfuni"
#names(res)[3] = "lbfav"

res2  = melt(res,id.vars=c("rho","beta1","beta2","ref","lbf2given1","lbf1given2","lbfuni.2","lbfuni.1"),variable_name="Test")
res2 = res2[sample(1:nrow(res2)),]

#note res1 was simulated with effect sizes mm=5 times bigger than res2
png("BF_vs_ref_rho04_panela.png",height=240,width=720)
p=ggplot(data=subset(res1,rho==0.4 & (beta1 < 0.1) ),aes(Y1,Y2,colour=Genotype))+
 #   opts(title="Y1 effect size") +
  geom_point(size=1)+
  facet_grid(rho~beta1)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()

png("BF_vs_ref_rho07_panela.png",height=240,width=720)
p=ggplot(data=subset(res1,rho==0.7 & (beta1!=0.08) ),aes(Y1,Y2,colour=Genotype))+
#    opts(title="Y1 effect size") +
  geom_point(size=1)+
  facet_grid(rho~beta1)
#+scale_colour_gradient(low="red",high="blue")
print(p+scale_y_continuous(limits=c(-4,6))+scale_x_continuous(limits=c(-4,6))+coord_equal(ratio=1))
dev.off()


png("BF_vs_ref_rho04_panelb.png",height=240,width=720)
p=ggplot(data = subset(res2,rho==0.4 & (beta1 < 0.1)), aes(ref,value,color=Test,group=1))+
  scale_color_manual(values=c("red","blue","black")) +
  xlab("log(BFref)")+
  ylab("log(BF)")+
  geom_abline(colour = "black")+
  geom_point(alpha = 1,size=2)+#,colour="black")
  facet_grid(rho~beta1) 
print(p + coord_equal(ratio=1))
dev.off()

png("BF_vs_ref_rho07_panelb.png",height=240,width=720)
p=ggplot(data = subset(res2,rho==0.7 & (beta1 != 0.08)), aes(ref,value,color=Test,group=1))+
  scale_color_manual(values=c("red","blue","black")) +
  xlab("log(BFref)")+
  ylab("log(BF)")+
  geom_abline(colour = "black")+
  geom_point(alpha = 1,size=2)+#,colour="black")
  facet_grid(rho~beta1) 
print(p + coord_equal(ratio=1))
dev.off()




pdf("uni_vs_all.pdf")
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = res, aes(lbfuni, lbfall))+
  xlab("BFuni")+
  ylab("BFall")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()

png("uni_vs_all_subset_rho04.png",height=240,width=960)
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = subset(res,rho==0.4), aes(lbfuni, lbfall))+
  xlab("BFuni")+
  ylab("BFall")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()

png("uni1_vs_all_subset_rho04.png",height=240,width=960)
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = subset(res,rho==0.4), aes(lbfuni.1, lbfall))+
  xlab("BFuni.1")+
  ylab("BFall")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()


png("bf2given1_vs_uni_subset_rho04.png",height=240,width=960)
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = subset(res,rho==0.4), aes(lbfuni, lbf2given1))+
  xlab("BFuni")+
  ylab("BF2given1")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()

png("bf2given1_vs_all_subset_rho04.png",height=240,width=960)
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = subset(res,rho==0.4), aes(lbfall, lbf2given1))+
  xlab("BFall")+
  ylab("BF2given1")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()


png("uni2_vs_all_subset_rho04.png",height=240,width=960)
#qplot(lbfuni,lbfall,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFall")
p=ggplot(data = subset(res,rho==0.4), aes(lbfuni.2, lbfall))+
  xlab("BFuni.2")+
  ylab("BFall")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()


pdf("uni_vs_av.pdf")
#qplot(lbfuni,lbfav,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFav")
p=ggplot(data = res, aes(lbfuni, lbfav))+
  xlab("BFuni")+
  ylab("BFav")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()

pdf("all_vs_av.pdf")
#qplot(lbfuni,lbfav,data=res,facets = rho ~ beta1,xlab="BFuni",ylab="BFav")
p=ggplot(data = res, aes(lbfall, lbfav))+
  xlab("BFall")+
  ylab("BFav")+
  geom_abline(colour = "red")+
  geom_point(alpha = 0.6,colour="black")+
  facet_grid(rho~beta1) 
print(p)
dev.off()


