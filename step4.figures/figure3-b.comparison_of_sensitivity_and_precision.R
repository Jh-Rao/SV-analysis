library(ggplot2)

options(digits=4)
options(scipen = 1)

args<-commandArgs(TRUE)

input <- args[1]
outpt <- args[2]

data <- read.table(input, head = T, col.names = c("Software", "Alogrithm", "Key", "Type", "DNBSEQ", "Illumina"))
# Software        Alogrithm  Key     Type      DNBSEQ     Illumina
# BASIL-ANISE     RP-AS   Call    INS       6020.62 4246.00

#data$Key = gsub('Recall', 'Sensitivity', data$Key)
#data$Key = gsub('Precis', 'Precision', data$Key)


data$Alogrithm <- factor(data$Alogrithm,
		      levels = c("RP", "SR", "RD", "AS", "RP-SR", "RP-RD", "RP-AS", "RP-SR-AS", "RP-SR-RD", "Other"))
data$Type   <- factor(data$Type,
		      levels = c("DEL", "DUP", "INS", "INV", "TRA",
				 "DEL_SS", "DEL_S", "DEL_M", "DEL_L", "DUP_S", "DUP_M", "DUP_L", "INS_A", "INV_A"))
data$Key    <- factor(data$Key,
		      levels = c("Call", "Sensitivity", "Precision", "F1"))

# call
nd <- subset(data, Key == "Call")
p <- ggplot(nd, aes(x = DNBSEQ, y = Illumina, colour = Alogrithm )) +
     geom_point(size = 1.8) +
     geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey") +
     stat_smooth(formula=y~x,aes(group=1),method="lm",se=FALSE, linetype = "solid", color = "blue") +
     facet_grid(. ~ Type) +
     scale_x_log10() +
     scale_y_log10() 
ggsave(paste(outpt, ".call.pdf", sep = ""), height = 3.5, width = 8)

# no call
nd <- subset(data, Key != "Call")

#  calculate p
id = 0
hf = ha = hc = hsv = hke = hp = hdf = hmm = him = rep(0, length(unique(nd$Type)) * length(unique(nd$Key)))
for (svt in unique(nd$Type)){
  for (ke in unique(nd$Key)){
    dm = subset(subset(nd, Type == svt), Key == ke)$DNBSEQ
    di = subset(subset(nd, Type == svt), Key == ke)$Illumina
    if (length(dm) > 0 && length(di) > 0){
      wilcoxtest = wilcox.test(dm, di)
      p = wilcoxtest$p.value
      df = log(mean(dm)/mean(di)) / log(2)
      k = kappa(dm, di)
      cor = cor.test(dm,di)$estimate
    }else{
      p = 1
      df = 0
    }

    id = id + 1
    hsv[id] = svt
    hke[id] = ke
    hc[id]  = floor(cor * 10000 + 0.5) / 10000
    hp[id]  = floor(p * 10000 + 0.5) / 10000
    hdf[id] = df
    hmm[id] = mean(dm)
    him[id] = mean(di)
    ha[id] = "RP"
    if(p<0.01){
      hf[id] = " **"
    }else if(p < 0.05){
      hf[id] = " *"
    }else{
      hf[id] = ""
    }
  }
}
statdata = data.frame(Type = hsv, Key = hke, DNBSEQM = hmm, ILLM = him, difflog = hdf, p = hp, Cor = hc, Algorithm = ha, flag = hf)
statdata$label = paste("rho = ", statdata$Cor, sep = "")
write.table(statdata, paste(outpt, ".nocall.wilcoxtestp.xls", sep = ""))

statdata$Type   <- factor(statdata$Type,
                     levels = c("DEL", "DUP", "INS", "INV", "TRA"))

#                     levels = c("DEL_SS", "DEL_S", "DEL_M", "DEL_L", "DUP_S", "DUP_M", "DUP_L", "INS_A", "INV_A"))
statdata$Key    <- factor(statdata$Key,levels = c("Call", "Sensitivity", "Precision", "F1"))


p <- ggplot(nd, aes(x = DNBSEQ, y = Illumina, colour = Alogrithm )) +
     geom_point(size = 1.8) +
     geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey") +
     geom_text(x = 1, y = 100, aes(label = label), data = statdata, color = "black", vjust = 1, hjust = 0, size = 3) +
     facet_grid(Key ~ Type) +
     xlim(0, 100) +
     ylim(0, 100) +
     theme(strip.text = element_text(size = 14),
           legend.position = "bottom",
           axis.title.x = element_text(size = 14),
           axis.title.y = element_text(size = 20)) +
     guides(fill = guide_legend(ncol = 5))
ggsave(paste(outpt, ".nocall.pdf", sep = ""), height = 6, width = 10)

