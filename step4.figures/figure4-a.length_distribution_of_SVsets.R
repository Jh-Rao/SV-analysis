library(ggplot2)

args<-commandArgs(TRUE)

inp <- args[1]
oup <- args[2]

d <- read.table(inp, head =T)

p <- ggplot(d, aes(x = svsize, color = source)) +
     geom_density() +
     facet_grid(. ~ svtype) +
     scale_x_log10() +
     geom_vline(xintercept = c(300, 6000), col = "grey", size = 0.5, linetype = 2)

ggsave(oup,height = 4, width= 12)

d <- subset(d, svtype == "All")
p <- ggplot(d, aes(x = svsize, color = source)) +
     geom_density() +
     scale_x_log10() +
     geom_vline(xintercept = c(300, 6000), col = "grey", size = 0.5, linetype = 2) +
     labs(x = '', y = 'Density of SV (%)') +
     theme_classic()
ggsave(paste(oup, ".all.pdf", sep = ""), height = 2, width = 4)

