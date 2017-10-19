# Usage: "Rscript box_scatter_MA_plots_2.0.R expression_matrix.csv outname png_or_pdf rows cols width height resolution samplename1 samplename2..."

argu=commandArgs(trailingOnly=TRUE)

infile=argu[1]
outname=argu[2]
filetype=argu[3]
rows=as.numeric(argu[4])
cols=as.numeric(argu[5])
vidd=as.numeric(argu[6])
hojd=as.numeric(argu[7])
reso=as.numeric(argu[8])

sample_names=argu[9:(length(argu))]
date=Sys.Date()
data<-read.csv(infile, sep=',', header=F)
colnames=c("gene", sample_names)
num_samp=length(sample_names)
num_stop=1+num_samp

library(grid)

################ Make summary plot

print("Making summary plot...")
print("")

if(filetype=="pdf") {
pdf(paste("Summary_plot_",outname,"_",date,".",filetype, sep=""), width=12, height=8)
} else if(filetype=="png") {
png(paste("Summary_plot_",outname,"_",date,".",filetype, sep=""), width=12, height=8, units="in", res=reso)
}
boxplot(log2(data[,2:ncol(data)]+1), xaxt='n', ylab="log2(expr+1)")
axis(1, at=1:num_samp, labels=sample_names, las=2, cex.axis=0.5)
dev.off()






################ Make scatterplots
num_on_page=(rows*cols)

##### If PDF
if(filetype=="pdf") {

pdf(paste("Scatter_plots_",outname,"_",date,".pdf", sep=""), width=vidd, height=hojd)
par(mfrow=c(rows,cols))

# For every sample in the set...
for(i in (2:num_stop))
        {
        sample=i-1
        print(paste("Scatterplots sample ", sample, sep=""))

        # Plot it against every other sample in the set (including itself)
        for(j in (2:num_stop))
                {
                plot(x=log2(data[,i]+1), y=log2(data[,j]+1), xlim=c(0,20), ylim=c(0,20), xlab=paste("log2(expr+1)_", colnames[i],sep=""), ylab=paste("log2(expr+1)_",colnames[j],sep=""))
                abline(a=0, b=1, col="blue")
                }

	if(num_on_page-num_samp>0){
	for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
	}
	}
dev.off()
}




##### If PNG
if(filetype=="png") {

# For every sample in the set...
for(i in (2:num_stop))
        {
        sample=i-1
        print(paste("Scatterplots sample ", sample, sep=""))

	# Set up a png file 
	png(paste("Scatter_plots_",outname,"_",date,"_",sample,".png", sep=""), width=vidd, height=hojd, units="in", res=reso)	
	par(mfrow=c(rows,cols))

        # Plot it against every other sample in the set (including itself)
        for(j in (2:num_stop))
                {
                plot(x=log2(data[,i]+1), y=log2(data[,j]+1), xlim=c(0,20), ylim=c(0,20), xlab=paste("log2(expr+1)_", colnames[i],sep=""), ylab=paste("log2(expr+1)_",colnames[j],sep=""))
                abline(a=0, b=1, col="blue")
                }
	if(num_on_page-num_samp>0){
        	for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }	
		}
	dev.off()
	}	

}

print("")



################ Make MA-plots

##### If PDF
if(filetype=="pdf") {
pdf(paste("MA_plots_",outname,"_",date,".pdf", sep=""), width=vidd, height=hojd)
par(mfrow=c(rows,cols))

# For every sample in the set...
for(i in (2:num_stop))
        {
        sample=i-1
        print(paste("MA_plots sample ", sample, sep=""))

        # Plot it against every other sample in the set (including itself)
	for(j in (2:num_stop))
		{
                m = log2((data[,i]+1)/(data[,j]+1))
#                a = log2((data[,i]+data[,j])/2)
		a = ((log2(data[,i]+1))+(1+log2(data[,j]+1)))/2

                plot(x=a, y=m, xlim=c(0,20), ylim=c(-10,10), xlab=paste("Avr.expr ", colnames[i], "+", colnames[j], sep=""), ylab= paste("Fold change ", colnames[i], " vs ",colnames[j], sep="")  )
                abline(h=0,col='red')
                abline(h=0.5, col='blue')
                abline(h=-0.5, col='blue')
		}
	if(num_on_page-num_samp>0){
        	for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
		}
        }
dev.off()
}


##### If PNG
if(filetype=="png") {
# For every sample in the set...
for(i in (2:num_stop))
        {
        sample=i-1
        print(paste("MA-plots sample ", sample, sep=""))

        # Set up a png file
        png(paste("MA_plots_",outname,"_",date,"_",sample,".png", sep=""), width=vidd, height=hojd, units="in", res=reso)
        par(mfrow=c(rows,cols))

        # Plot it against every other sample in the set (including itself)
        for(j in (2:num_stop))
                {
                m = log2((data[,i]+1)/(data[,j]+1))
#                a = log2((data[,i]+data[,j])/2)
		a = ((log2(data[,i]+1))+(log2(data[,j]+1)))/2
                plot(x=a, y=m, xlim=c(0,20), ylim=c(-10,10), xlab=paste("Avr.expr ", colnames[i], "+", colnames[j], sep=""), ylab= paste("Fold change ", colnames[i], " vs ",colnames[j], sep="")  )
                abline(h=0,col='red')
                abline(h=0.5, col='blue')
                abline(h=-0.5, col='blue')
                }
	if(num_on_page-num_samp>0){
        	for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
		}
        dev.off()
        }
}
