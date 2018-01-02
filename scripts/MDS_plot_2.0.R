argu=commandArgs(trailingOnly=TRUE)
k=argu[1]
file=argu[2]
vfarg=argu[3:length(argu)]
outpart=sub(".csv", "", file)
date=Sys.Date()

pdf(paste("MDS_plot_", outpart, "_", date ,".pdf", sep=""), height=4, width=10)
par(mfrow=c(1,3))

for(i in 0:k)
	{
	# Read in data for euclidian and manhattan
	if(i==0)        {       filnamn="all.csv"       }
        else    {       filnamn=paste(i, ".csv", sep="")        }
	rawdata1<-read.csv(filnamn, sep=',', header=F)
	rawdata<-rawdata1[,-20]
	rm(rawdata1)
	
	# Read in data for pearson
	if(i==0)	{	filnamn="all.csv"	}
	else	{	filnamn=paste(i, ".csv", sep="")	}
        rawdata_P<-read.csv(filnamn, sep=',')
        data_P<-rawdata_P[,-1]
        sim<-cor(data_P, method="pearson")
        avst_P<-(1-sim)
        rm(data_P, rawdata_P, sim)

	# Set row and column names
	samples1<-rawdata[1,]
	samples2<-samples1[-1]
	samples<-as.vector(as.matrix(samples2))
	rm(samples1)
	rm(samples2)
	genes1<-rawdata[,1]
	genes<-genes1[-1]
	rm(genes1)

	# Transpose data
	inter<-rawdata[-1,-1]
	data<-t(inter)
	rm(rawdata, inter)
	colnames(data)<-genes
	rownames(data)<-samples

	# Do analysis

	for(j in c("euclidean", "manhattan", "pearson"))
		{
		if(j!="pearson")	{	avst<-dist(data, method=j)	}
		else	{	avst<-avst_P	}
		analys<-cmdscale(avst, eig=TRUE, k=2)
		x<-analys$points[,1]
		y<-analys$points[,2]

		if(i==0)	{	plotname=paste("All genes (", j, ")",  sep="")	}
		else	
			{
			num=((1000*i)+1)-1000	
			plotname=paste("Genes ", num, " to ", 1000*i, " (", j, ")", sep="")
			}
		plot(x,y, xlab="Coord_1", ylab="Coord_2", main=plotname, type="n")
		text(x,y, labels=row.names(data), cex=0.6, col=vfarg, font=2)
		}
	}
dev.off()
