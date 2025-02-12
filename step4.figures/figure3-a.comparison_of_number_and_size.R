library(ggplot2)
options(digits=4)

args<-commandArgs(TRUE)

input <- args[1]
outpt <- args[2]

data <- read.table(input, head = T)
# Software        Algorithm       Key     Type    DNBSEQ          Illumina
# BASIL-ANISE     RP-AS           Number  INS     6020.6250       4246.0000

data$Algorithm <- factor(data$Algorithm,
                      levels = c("RP", "SR", "RD", "AS", "RP-SR", "RP-RD", "RP-AS", "RP-SR-AS", "RP-SR-RD", "Other"))
data$Key <- factor(data$Key, levels = c("Number", "Length"))
data$Type <- factor(data$Type, levels = c("DEL", "DUP", "INS", "INV", "TRA"))

#data$Type = gsub('DEL', 'DEL (n = 33)', data$Type)
#data$Type = gsub('DUP', 'DUP (n = 22)', data$Type)
#data$Type = gsub('INS', 'INS (n = 24)', data$Type)
#data$Type = gsub('INV', 'INV (n = 17)', data$Type)
#data$Type = gsub('TRA', 'TRA (n = 8)', data$Type)

# calculate p and cor
id = 0
hp2 = hn = ha = hs = hk = hc = hp = rep(0, length(unique(data$Type)) * length(unique(data$Key)))
for (svt in unique(data$Type)){
  for (k in unique(data$Key)){
    dm  = subset( subset(data, Type == svt), Key == k)$DNBSEQ
    di  = subset( subset(data, Type == svt), Key == k)$Illumina
    hnk = length(dm)

    if(length(dm) * length(di) > 0){
      cor = signif(cor.test(dm, di, method = "spearman")$estimate, 4)
      p2  = signif(cor.test(dm, di, method = "spearman")$p.value, 4)
      p   = signif(wilcox.test(dm, di)$p.value, 4)
    }else{
      cor = "-"
      p = "-"
    }

      id = id + 1
      hs[id] = svt
      hk[id] = k
      hc[id] = cor
      hp[id] = p
      ha[id] = "RP"
      hn[id] = hnk
      hp2[id] = p2
  }
}
statdata = data.frame(Key = hk, Type = hs, Cor = hc, p = hp, Algorithm = ha, Num = hn, Corp = hp2)
statdata$label = paste("n = ", statdata$Num, "\nrho = ", statdata$Cor, "\n", sep = "")
statdata$label[id] = ""
statdata$Key <- factor(statdata$Key, levels = c("Number", "Length"))
statdata$Type <- factor(statdata$Type, levels = c("DEL", "DUP", "INV", "INS", "TRA"))
write.table(statdata, paste(outpt, ".cor.xls", sep = ""))
#for (i in 1:nrow(data)){
#  for (j in 1:nrow(statdata)){
#    if (statdata[j, 1] == data$Key[i] && statdata[j, 2] == data$Type[i]){
#      data$label = paste("cor = ", statdata[j, 3], " P = ", statdata[j, 4], sep = "")
#    }
#  }
#}

bk = c(0, 2, 5, 8)
lb = c("0", "100", "100k", "100M")

p <- ggplot(data, aes(x = DNBSEQ, y = Illumina, colour = Algorithm)) +
     geom_point(size = 2) +
     stat_smooth(formula=y~x,aes(group=1),method="lm",se=FALSE, linetype = "solid", size = 0.5, color = "grey") +
     geom_text(x = 1, y = 8.5, aes(label = label), data = statdata, color = "black", vjust = 1, hjust = 0) +
     facet_grid(Key ~ Type) +
     scale_x_log10() +
     scale_y_log10() +
     theme(strip.text = element_text(size = 14),
	   legend.position = "bottom",
	   axis.title.x = element_text(size = 14),
	   axis.title.y = element_text(size = 20)) +
     guides(fill = guide_legend(ncol = 5))
ggsave(outpt, height = 6, width = 10)

library(tidyr)
nd <- data %>% gather(Platform, Value, -Software, -Algorithm, -Key, -Type)
p <- ggplot(nd, aes(x = Platform, y = Value)) +
     geom_boxplot(aes(fill = Platform, color = Platform), width = 0.1) +
     facet_grid(Key ~ Type, scales = "free")+
     xlab("") + ylab("") +
     scale_y_log10()
ggsave(paste(outpt, ".boxplot.pdf", sep = ""), height = 6, width = 10)

