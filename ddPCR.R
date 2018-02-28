setwd("/Users/kuanghe/Desktop/R/ddPCR")

allfiles <- list.files()

csvfiles <- allfiles[grep(pattern="csv", x=allfiles)]

newrun <- read.csv(file = csvfiles, header=F, stringsAsFactors=F)

#assign headers
colnames(newrun) = newrun[1,]
#delete header row
newrun <- newrun[-1,]
#convert characters to numeric values for copy numbers
newrun$CopiesPer20uLWell <- as.numeric(newrun$CopiesPer20uLWell)



temps <- c(63, 62.4, 61.1, 59.2, 57, 55.1, 53.7, 53)
alltemps <-c()
for (i in 1: 4)
{
  alltemps <- c(alltemps,temps)  
}

#add a column of annealing temperatures
newrun$AnnealTemp <- alltemps

#subset the 4 columns in a 96-well plate
mLINE62bp_mouse <- newrun[1:8,]
mLINE62bp_blank <- newrun[9:16,]
mLINE58bp_mouse <- newrun[17:24,]
mLINE58bp_blank <- newrun[25:32,]

#plot copy number vs annealing temperature into a single pdf file
pdf(paste0("kh2-111", ".pdf"))

plot(y=mLINE62bp_mouse$CopiesPer20uLWell, x=mLINE62bp_mouse$AnnealTemp, ylim=c(0,5000), xlab = "Annealing Temperature", ylab = "Copy Number", pch=16, cex=2)
lines(y=mLINE62bp_mouse$CopiesPer20uLWell, x=mLINE62bp_mouse$AnnealTemp, type = "l")
title(main = "Mouse LINE-1 62 bp Target \n 20 ng HaeIII-digested mouse gDNA")
plot(y=mLINE62bp_blank$CopiesPer20uLWell, x=mLINE62bp_blank$AnnealTemp, xlab = "Annealing Temperature", ylab = "Copy Number", pch=16, cex=2)
lines(y=mLINE62bp_blank$CopiesPer20uLWell, x=mLINE62bp_blank$AnnealTemp, type = "l")
title(main = "Mouse LINE-1 62 bp Target \n No Template Control")
plot(y=mLINE58bp_mouse$CopiesPer20uLWell, x=mLINE58bp_mouse$AnnealTemp, ylim=c(0,5000), xlab = "Annealing Temperature", ylab = "Copy Number", pch=16, cex=2)
lines(y=mLINE58bp_mouse$CopiesPer20uLWell, x=mLINE58bp_mouse$AnnealTemp, type = "l")
title(main = "Mouse LINE-1 58 bp Target \n 20 ng HaeIII-digested mouse gDNA")
plot(y=mLINE58bp_blank$CopiesPer20uLWell, x=mLINE58bp_blank$AnnealTemp, xlab = "Annealing Temperature", ylab = "Copy Number", pch=16, cex=2)
lines(y=mLINE58bp_blank$CopiesPer20uLWell, x=mLINE58bp_blank$AnnealTemp, type = "l")
title(main = "Mouse LINE-1 58 bp Target \n No Template Control")

#calculate signal to noise ratios for each annealing temperature. S/N=copy numbers of target in 20 ng mouse gDNA divided by those of same target in IDTET
S_to_N_62bp <- c()
S_to_N_58bp <- c()
S_to_N_62bp <- c(S_to_N_62bp, mLINE62bp_mouse$CopiesPer20uLWell/mLINE62bp_blank$CopiesPer20uLWell)
S_to_N_58bp <- c(S_to_N_58bp, mLINE58bp_mouse$CopiesPer20uLWell/mLINE58bp_blank$CopiesPer20uLWell) 

#plot signal to noise ratios
plot(y=S_to_N_62bp, x=temps, ylim = c(0, 200), xlab = "Annealing Temperature", ylab = "Signal to Noise Ratio", pch=16, cex=2)
lines(y=S_to_N_62bp, x=temps, type = "l")
title(main = "Mouse LINE-1 62 bp Target")
plot(y=S_to_N_58bp, x=temps, ylim = c(0, 200), xlab = "Annealing Temperature", ylab = "Signal to Noise Ratio", pch=16, cex=2)
lines(y=S_to_N_58bp, x=temps, type = "l")
title(main = "Mouse LINE-1 58 bp Target")

dev.off()