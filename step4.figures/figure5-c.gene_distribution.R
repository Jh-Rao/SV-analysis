library(ggplot2)

args<-commandArgs(TRUE)

dfl <- args[1]
ifl <- args[2]
oup <- args[3]


dd <- read.table(dfl, head=T)
di <- read.table(ifl, head=T)

n=length(unique(dd$Gene))
m=2*n

frame = data.frame(platform = rep(0, m),svtype = rep(0, m), value = rep(0, m),rate=rep(0,m))

for (i in 1:length(unique(dd$Gene))){
  frame$platform[i] = "DNBSEQ"
  frame$svtype[i] = names(table(dd$Gene))[i]
  frame$value[i]  = as.numeric(table(dd$Gene))[i]
  frame$rate[i] = as.numeric(table(dd$Gene))[i] / sum(as.numeric(table(dd$Gene))) * 100

  j=i+n
  frame$platform[j] = "Illumina"
  frame$svtype[j] = names(table(di$Gene))[i]
  frame$value[j]  = as.numeric(table(di$Gene))[i]
  frame$rate[j] = as.numeric(table(di$Gene))[i] / sum(as.numeric(table(di$Gene))) * 100

}
frame$svtype = factor(frame$svtype, levels = c("non-coding", "5-UTR", "3-UTR", "promoter-TSS", "TTS", "exon", "intron", "Intergenic"))

p <- ggplot(frame, aes(svtype, weight = rate, fill = platform)) +
     geom_bar(width = .7, position = 'dodge') +
     labs(y = '% of SVs', x = '') +
     theme_classic() +
     theme(axis.text.x = element_text(size = 7, angle = 60, vjust = 0.5))
ggsave(oup, height = 2, width = 5)

