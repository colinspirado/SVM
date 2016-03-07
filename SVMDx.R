#Support vector machine classification analysis to identify EMCI at risk for conversion from baseline diagnosis within 3 years of baseline
#Training data set to include LMCI
#Model performance to be evaluated on EMCI
#Inputs:
	#Ventricle and whole brain volume as a (scaled) proportion of intracranial volume
	#Hippocampal, Entorhinal, Fusiform and MidTemporal volume as a proportion of whole brain volume

library(ADNIMERGE)
library(ggplot2)
library(lme4)
library(e1071)

#Identify converters and non-converting completers; reversion events ignored
dxconv=subset(adnimerge,VISCODE%in%'bl'&DX.bl%in%c('EMCI','LMCI')&
	RID%in%unique(subset(adnimerge,VISCODE%in%c('m03','m06','m12','m18','m24','m30','m36')&
		DX%in%c('NL to Dementia','MCI to Dementia','Dementia'))$RID),
	c(RID,DX.bl))
comp=subset(adnimerge,VISCODE%in%'m36'&DX.bl%in%c('EMCI','LMCI')&!RID%in%dxconv$RID,c(RID,DX.bl))
table(dxconv$DX.bl)
table(comp$DX.bl)

#Create training and testing data sets
train=subset(adnimerge,VISCODE%in%'bl'&DX.bl%in%c('LMCI')&RID%in%c(dxconv$RID,comp$RID),
	c(RID,DX.bl,Ventricles.bl,Hippocampus.bl,WholeBrain.bl,Entorhinal.bl,Fusiform.bl,MidTemp.bl,ICV.bl))
	train$cv=ordered(as.numeric(train$RID%in%comp$RID)+1)   #1=T, 2=F
	train[,3:9]=lapply(train[,3:9],as.numeric)
	train=na.omit(train)
test=subset(adnimerge,VISCODE%in%'bl'&DX.bl%in%'EMCI'&RID%in%c(dxconv$RID,comp$RID),
	c(RID,DX.bl,Ventricles.bl,Hippocampus.bl,WholeBrain.bl,Entorhinal.bl,Fusiform.bl,MidTemp.bl,ICV.bl))
	test$cv=ordered(as.numeric(test$RID%in%comp$RID)+1)   #1=T, 2=F
	test[,3:9]=lapply(test[,3:9],as.numeric)
	test=na.omit(test)

#Process volumes as proportions
k=1
for(j in c(4,6:8,3,5)){for(t in c('test','train'))
	{
	foo=eval(as.name(t))
	foo[,j]=foo[,j]/foo[,c(5,5,5,5,9,9)[k]]
	assign(eval(t),foo)
	}; k=k+1}; rm(k,foo)

#Scatterplot matrix
pairs(formula=~Ventricles.bl+Hippocampus.bl+WholeBrain.bl+Entorhinal.bl+Fusiform.bl+MidTemp.bl,data=train,cex=0.7,lwd=0.5,
	labels=c('V','H','WB','E','F','MT'),col=c('red','blue')[as.numeric(train$cv)])   #red = cv

#Fit svms and catpure results in train
	#fitl=linear kernel tuned for cost; 71% in train; more misses than hits in test
	#fitp=polynomial kernel tuned for cost and degree; 73% in train; worst in test with notable bias toward predicted conversion relative to linear
	#fitr=radial kernel tuned for cost and gamma; 72% in train; slight bias toward predicted conversion relative to linear in test
	#fits=sigmoid kernel tuned for cost and gamma; 71% in train; equivalent to linear in test
fitl=svm(formula=cv~.,data=subset(train,,-c(RID,DX.bl,ICV.bl)),kernel='linear',cost=0.25,cross=10)
fitp=svm(formula=cv~.,data=subset(train,,-c(RID,DX.bl,ICV.bl)),kernel='polynomial',degree=2,cost=0.03125,gamma=1,coef0=1,cross=10)
fitr=svm(formula=cv~.,data=subset(train,,-c(RID,DX.bl,ICV.bl)),kernel='radial',cost=1000,gamma=0.001,cross=10)
fits=svm(formula=cv~.,data=subset(train,,-c(RID,DX.bl,ICV.bl)),kernel='sigmoid',cost=10,gamma=0.01,coef0=0,cross=10)
fit=fitr
train$fx=fit$decision.values[,1]   #values < 0 indicate predicted conversion
train$pred=fit$fitted
train$sv=row.names(train)%in%row.names(fit$SV)
with(train,table(cv,pred))

#Plot f(x) against inputs for train
op=par(mfrow=c(2,3))
for(j in 3:8)
	{
	plot(train$fx~train[,j],xlab=names(train)[j],ylab='f(x)',
		col=c(NA,'blue','black','red')[as.numeric(train$cv)+as.numeric(train$pred)],
		pch=c(1,3)[train$sv+1])
	abline(h=0)
	}
par(op)

#Plot scatterplot matrix color coded by predictions
pairs(formula=~Ventricles.bl+Hippocampus.bl+WholeBrain.bl+Entorhinal.bl+Fusiform.bl+MidTemp.bl,data=train,cex=0.7,lwd=0.5,
	labels=c('V','H','WB','E','F','MT'),col=c('red','blue')[as.numeric(train$pred)],pch=c(1,3)[with(train,as.numeric(cv)==as.numeric(pred))+1])

#Apply model to test
pred=predict(fit,subset(test,,-c(RID,DX.bl,ICV.bl)),decision.values=T)
test$fx=attr(pred,'decision.values')[,1]
test$pred=pred

#Summarize and plot results for test
with(test,table(cv,pred))
op=par(mfrow=c(2,3))
for(j in 3:8)
	{
	plot(test$fx~test[,j],xlab=names(test)[j],ylab='f(x)',
		col=c(NA,'blue','black','red')[as.numeric(test$cv)+as.numeric(test$pred)],
		pch=1)
	abline(h=0)
	}
par(op)
pairs(formula=~Ventricles.bl+Hippocampus.bl+WholeBrain.bl+Entorhinal.bl+Fusiform.bl+MidTemp.bl,data=test,cex=0.7,lwd=0.5,
	labels=c('V','H','WB','E','F','MT'),col=c('red','blue')[as.numeric(test$pred)],pch=c(1,3)[with(test,(as.numeric(cv)==as.numeric(pred))+1)])
