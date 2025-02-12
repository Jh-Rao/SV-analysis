library(ggplot2)
library(RColorBrewer)

args<-commandArgs(TRUE)

inp <- args[1]
oup <- args[2]

d <- read.table(inp, head = T)
col = rev(c(brewer.pal(9, "Set1"), "black"))

d <- subset(d, svtype != "All")
d$reptype <- factor(d$reptype, levels = rev(c("TRF", "STR", "Alu", "L1", "SVA", "HERV", "LTR", "Complex", "Other", "NoRepeat")))

p <- ggplot(d, aes(x = gcbin, y = repcount, fill = reptype)) +
     geom_bar(width = 0.75, position = "stack", stat = 'identity') +
     facet_wrap(source ~ svtype, scales = "free_y", ncol = 4) +
     xlim(-1, 21) +
     scale_fill_manual(values = col) +
     labs(x = "GC (%)", y = "Count of SVs") +
     scale_x_continuous(breaks = c(0,4,8,12,16,20),
                        labels = c("0","20", "40", "60", "80","100"),)
ggsave(oup, height = 5, width = 12)

