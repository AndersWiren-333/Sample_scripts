# Description: Makes three kinds of plots (boxplot, pairwise scatterplots, pairwise MA-plots) to visualise global gene expression levels in a data set.
# Input arguments:
#	infile.csv (a csv file containing the data set)
#	outname (desired basename for resulting plot files)
#	filetype (which filetype the output plots should be. "pdf" or "png")
#	rows (number of rows of individual plots on an output page)
#	cols (number of columns of individual plots on an output page)
#	vidd (width of output page)
#	hojd (height of output page)
#	reso (resolution of output page)
#	sample_name_1, sample_name_2 etc. (names of samples in the dataset, will be used for labels and axes)
# Usage: "Rscript box_scatter_MA_plots_2.0.R expression_matrix.csv outname png_or_pdf rows cols width height resolution samplename1 samplename2..."
# end description

########################################## Universal R script header ##########################################

# Load libraries (= packages = modules) that this script depends on
library(grid)

# Read command line arguments
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

# Get date
date=Sys.Date()

# end header

########################################## Define local functions ##########################################

# end functions

########################################## Processing ##########################################

data<-read.csv(infile, sep=',', header=F)
numgenes <- nrow(data)
colnames=c("gene", sample_names)
num_samp=length(sample_names)
num_stop=1+num_samp

library(grid)

# Set parameters for output plots
cex_main <- 7
cex_axis <- 5
cex_n <- 5
cex_lab <- 5
summar <- c(25,10,10,5)		# Summary plot margins
scatmar <- c(10,10,10,10)	
mamar <- c(10,10,10,10)

################ Make summary plot

print("Making summary plot...")
print("")

if(filetype=="pdf")
	{
	pdf(paste("Summary_plot_",outname,"_",date,".",filetype, sep=""), width=65, height=40)
	} else if(filetype=="png")
	{
	png(paste("Summary_plot_",outname,"_",date,".",filetype, sep=""), width=65, height=40, units="in", res=reso)
	}

par(mar=summar)
boxplot(log2(data[,2:ncol(data)]+1), xaxt='n', xlab="Samples", ylab="log2(expr+1)", main=paste("Log2 expression levels of ", numgenes, " genes", sep=""),
cex.main=cex_main, cex=cex_n, cex.lab=cex_lab, cex.axis=cex_axis)
axis(1, at=1:num_samp, labels=sample_names, las=2, cex=cex_n, cex.lab=cex_lab, cex.axis=cex_axis)

dev.off()

################ Make scatterplots
num_on_page=(rows*cols)

##### If PDF
if(filetype=="pdf")
	{
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
			par(mar=scatmar)
			plot(x=log2(data[,i]+1), y=log2(data[,j]+1), xlim=c(0,20), ylim=c(0,20), xlab=paste("log2(expr+1)_", colnames[i],sep=""), ylab=paste("log2(expr+1)_",colnames[j],sep=""), cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			abline(a=0, b=1, col="blue")
			#axis(1, cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			#axis(2, cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			}
			
		if(num_on_page-num_samp>0)
			{
			for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
			}
		}
	dev.off()
	}

##### If PNG
if(filetype=="png")
	{

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
			par(mar=scatmar)
			plot(x=log2(data[,i]+1), y=log2(data[,j]+1), xlim=c(0,20), ylim=c(0,20), xlab=paste("log2(expr+1)_", colnames[i],sep=""), ylab=paste("log2(expr+1)_",colnames[j],sep=""), cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			abline(a=0, b=1, col="blue")
			#axis(1, cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			#axis(2, cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			}

		if(num_on_page-num_samp>0)
			{
			for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }	
			}
			
		dev.off()
		}	
	}

print("")


################ Make MA-plots

##### If PDF
if(filetype=="pdf")
	{
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
#	     	a = log2((data[,i]+data[,j])/2)
			a = ((log2(data[,i]+1))+(1+log2(data[,j]+1)))/2

			par(mar=mamar)
			plot(x=a, y=m, xlim=c(0,20), ylim=c(-10,10), xlab=paste("Avr.expr ", colnames[i], "+", colnames[j], sep=""), ylab= paste("Fold change ", colnames[i], " vs ",colnames[j], sep=""),
			cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			abline(h=0,col='red')
			abline(h=0.5, col='blue')
			abline(h=-0.5, col='blue')
			}
			
		if(num_on_page-num_samp>0)
			{
			for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
			}
        }
	dev.off()
	}


##### If PNG
if(filetype=="png")
	{
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
#          	a = log2((data[,i]+data[,j])/2)
			a = ((log2(data[,i]+1))+(log2(data[,j]+1)))/2
			
			par(mar=mamar)
			plot(x=a, y=m, xlim=c(0,20), ylim=c(-10,10), xlab=paste("Avr.expr ", colnames[i], "+", colnames[j], sep=""), ylab= paste("Fold change ", colnames[i], " vs ",colnames[j], sep=""),
			cex=cex_n, cex.axis=cex_axis, cex.lab=cex_lab)
			abline(h=0,col='red')
			abline(h=0.5, col='blue')
			abline(h=-0.5, col='blue')
			}
			
		if(num_on_page-num_samp>0)
			{
			for(k in (1:(num_on_page-num_samp)))     {       plot.new()      }
			}
			
        dev.off()
        }
	}

# end processing
