
plotHistograms <- function(list1, list2, title, lxlab, filename,lenInSame, lenNotInSame){
    
    outFile = paste("~/data/DB/PomEnzNet/RMS_buffering/plots/pairEnz/",filename,".H2.all", sep="")
    png(outFile)
    
    hist(list1 , main=title, freq=FALSE, xlab=lxlab, breaks=50, col=rgb(1,0,0,0.5))
    
    hist(list2, freq=FALSE, breaks=100, col=rgb(0,0,1,0.5), add=T)
    
    legend("topleft", col=c("red", "blue"), legend=c(paste("InSameRMS: ",lenInSame), paste("NotInSameRMS: ",lenNotInSame)), lty=1, cex=0.8)
    
    res = t.test(list1, list2 )
    legend("topright", legend=paste("p-value: ", signif(res$p.value, digits = 3 ) ), cex=0.8)
    
    dev.off()
    
    return()
}



df <- read.table("~/data/DB/PomEnzNet/data/allEnzPairDataRNA.DB", header=FALSE)


lenInSame = length(df[ df$V15=='True', ]$V3)
lenNotInSame = length(df[ df$V15=='False', ]$V3)

# gmean N
plotHistograms(list1=df[ df$V15=='True', ]$V3, list2=df[ df$V15=='False', ]$V3, title="gmean sum RNA enzyme pairs - RMS H2 N", 
               lxlab = "log2(gmean(enz1+enz2))", filename="gmean_RMS_H2_N.png", lenInSame, lenNotInSame )

#sd N
plotHistograms(list1=df[ df$V15=='True', ]$V4, list2=df[ df$V15=='False', ]$V4, title="SD sum RNA enzyme pairs - RMS H2 N", 
               lxlab = "log2(SD(enz1+enz2))", filename="sd_RMS_H2_N.png", lenInSame, lenNotInSame )

# correlation N
plotHistograms(list1=df[ df$V15=='True', ]$V5, list2=df[ df$V15=='False', ]$V5, title="Correlation RNA enzyme pairs - RMS H2 N", 
               lxlab = "Pearson correlation score", filename="corr_RMS_H2_N.png", lenInSame, lenNotInSame )


# gmean H
plotHistograms(list1=df[ df$V15=='True', ]$V6, list2=df[ df$V15=='False', ]$V6, title="gmean sum RNA enzyme pairs - RMS H2 H", 
               lxlab = "gmean(log2(enz1+enz2))", filename="gmean_RMS_H2_H.png", lenInSame, lenNotInSame )

#sd H
plotHistograms(list1=df[ df$V15=='True', ]$V7, list2=df[ df$V15=='False', ]$V7, title="SD sum RNA enzyme pairs - RMS H2 H", 
               lxlab = "SD(log2(enz1+enz2))", filename="sd_RMS_H2_H.png", lenInSame, lenNotInSame )

# correlation H
plotHistograms(list1=df[ df$V15=='True', ]$V8, list2=df[ df$V15=='False', ]$V8, title="Correlation RNA enzyme pairs - RMS H2 H", 
               lxlab = "Pearson correlation score", filename="corr_RMS_H2_H.png", lenInSame, lenNotInSame )


#delta SD
plotHistograms(list1=df[ df$V15=='True', ]$V9, list2=df[ df$V15=='False', ]$V9, title="Delta SD - RMS H2", 
               lxlab = "SD(log2(enz1+enz2))_N - SD(log2(enz1+enz2))_H", filename="deltaSD.png", lenInSame, lenNotInSame )


# sd_reduction score N
plotHistograms(list1=df[ df$V15=='True', ]$V10, list2=df[ df$V15=='False', ]$V10, title="reduction score SD - RMS H2 N",
               lxlab = "SD(log2(enz1+enz2)) - min(SD(log2(enz1)), SD(log2(enz2)) )", filename="sdReductionScore_N.png", lenInSame, lenNotInSame )

# sd_reduction score H
plotHistograms(list1=df[ df$V15=='True', ]$V11, list2=df[ df$V15=='False', ]$V11, title="reduction score SD - RMS H2 H",
               lxlab = "SD(log2(enz1+enz2)) - min(SD(log2(enz1)), SD(log2(enz2)) )", filename="sdReductionScore_H.png", lenInSame, lenNotInSame )


#sd FC
plotHistograms(list1=df[ df$V15=='True', ]$V12, list2=df[ df$V15=='False', ]$V12, title="SD fold change enzyme pairs - RMS H2", 
               lxlab = "SD(FC(log2(enz1+enz2)))", filename="sdFC_RMS_H2.png", lenInSame, lenNotInSame )

#sdreduction score FC
plotHistograms(list1=df[ df$V15=='True', ]$V13, list2=df[ df$V15=='False', ]$V13, title="SD reduction score in fold change enzyme pairs - RMS H2", 
               lxlab = "SD(FC(log2(enz1+enz2))) - min[SD(FC(log2(enz1))), SD(FC(log2(enz2)))]", filename="sdReduction_FC_RMS_H2.png", lenInSame, lenNotInSame )

#mean FC
plotHistograms(list1=df[ df$V15=='True', ]$V14, list2=df[ df$V15=='False', ]$V14, title="Mean fold change enzyme pairs - RMS H2", 
               lxlab = "mean(FC(log2(enz1+enz2)))", filename="meanFC_RMS_H2.png", lenInSame, lenNotInSame )



