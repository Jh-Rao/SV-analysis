library(ggplot2)

args<-commandArgs(TRUE)

dfl <- args[1]
ifl <- args[2]
oup <- args[3]


dd <- read.table(dfl, head=T)
di <- read.table(ifl, head=T)

n = 10
m = n * 2

frame = data.frame(platform = rep(0, m),svtype = rep(0, m), value = rep(0, m))

for (i in 1:n){
  frame$platform[i] = "DNBSEQ"
  frame$svtype[i] = names(table(dd$Repeat))[i]
  frame$value[i]  = as.numeric(table(dd$Repeat))[i] / sum(as.numeric(table(dd$Repeat))) * 100

  j=i+n
  frame$platform[j] = "Illumina"
  frame$svtype[j] = names(table(di$Repeat))[i]
  frame$value[j]  = as.numeric(table(di$Repeat))[i] / sum(as.numeric(table(di$Repeat))) * 100
}

frame$svtype <- factor(frame$svtype, levels =  c("TRF", "STR", "Alu", "L1", "LTR", "SVA", "HERV", "Complex", "Other", "NoRepeat"))

p <- ggplot(frame, aes(svtype, weight = value, fill = platform)) +
     geom_bar(width = .7, position = 'dodge') +
     labs(y = '% of SVs', x = '') +
     theme_classic() +
     theme(axis.text.x = element_text(size = 7, angle = 60, vjust = 0.5)) +
     ylim(0, 50)
ggsave(oup, height = 2, width = 5)

