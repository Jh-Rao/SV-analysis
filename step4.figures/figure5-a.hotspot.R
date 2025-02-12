library(karyoploteR)
library(primatR)

args<-commandArgs(TRUE)

dfile <- args[1]
ifile <- args[2]
oup <- args[3]

dnb <- read.table(dfile, head = F, col.names = c("chr","start","end","strand"))
ill <- read.table(ifile, head = F, col.names = c("chr","start","end","strand"))

Gd <- makeGRangesFromDataFrame(dnb)
Gi <- makeGRangesFromDataFrame(ill)

seqlengths(Gd) <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)
seqlengths(Gi) <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560)

Hd <- hotspotter(Gd, bw=200000, num.trial = 1000)
Hi <- hotspotter(Gi, bw=200000, num.trial = 1000)

# output hotspot
write.table(data.frame(Hd), file = paste(oup, ".DNBSEQ.hotspot.xls",sep=""), sep = "\t")
write.table(data.frame(Hi), file = paste(oup, ".Illumina.hotspot.xls",sep=""), sep = "\t")

# draw
pdf(paste(oup, ".pdf", sep = ""), height = 4, width = 3)
chrlist <- paste0("chr", c(1:22, "X"))
kp <- plotKaryotype(plot.type = 2, genome = "hg19", chromosomes=chrlist)
kpPlotRegions(kp, Hd, data.panel = 1, col = rgb(248/255,118/255,109/255))
kpPlotRegions(kp, Hi, data.panel = 2, col = rgb(0,191/255,196/255))
dev.off()

