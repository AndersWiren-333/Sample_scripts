lista=c("testA.csv", "testB.csv")
list=c("qplot_01_B24_1.fastq.csv",
"qplot_02_B24_2.fastq.csv",
"qplot_03_B24_3.fastq.csv",
"qplot_04_B96_1.fastq.csv",
"qplot_05_B96_2.fastq.csv",
"qplot_06_B96_3.fastq.csv",
"qplot_07_B96_Low.fastq.csv",
"qplot_08_BLay_1.fastq.csv",
"qplot_09_BLay_2.fastq.csv",
"qplot_10_BLay_3.fastq.csv",
"qplot_11_O24_1.fastq.csv",
"qplot_12_O24_2.fastq.csv",
"qplot_13_O24_3.fastq.csv",
"qplot_14_O96_1.fastq.csv",
"qplot_15_O96_2.fastq.csv",
"qplot_16_O96_3.fastq.csv",
"qplot_17_OLay_1.fastq.csv",
"qplot_18_OLay_2.fastq.csv",
"qplot_19_OLay_3.fastq.csv",
"qplot_20_FB24_1.fastq.csv",
"qplot_21_FB24_2.fastq.csv",
"qplot_22_FB24_3.fastq.csv",
"qplot_23_FB96_1.fastq.csv",
"qplot_24_FB96_2.fastq.csv",
"qplot_25_FB96_3.fastq.csv",
"qplot_26_FBLay_1.fastq.csv",
"qplot_27_FBLay_2.fastq.csv",
"qplot_28_FBLay_3.fastq.csv")

for(i in list)
	{
	infile=i
	print(infile)
	rawdata<-read.csv(infile, sep=',', header=F)
	newdata=apply(rawdata, 2, summary)
	data=newdata[-4,]

	pdf(paste("plot_", i, ".pdf", sep=""), width=25, height=20)
	boxplot(data[,1:50], xaxt='n', main=paste("Quality scores per base, ", i,  sep=""), cex.main=2, cex.axis=1.5)
	axis(1, at=1:50, labels=c(1:50), las=1, cex.axis=1.5)
	dev.off()
	}
