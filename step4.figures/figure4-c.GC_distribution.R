library(ggplot2)

args<-commandArgs(TRUE)

inp <- args[1]
oup <- args[2]

d <- read.table(inp, head =F, col.names = c("id", "length", "gc", "source" ,"svtype"))

p <- ggplot(d, aes(x = gc, color = source)) +
     geom_density() +
     facet_grid(. ~ svtype) +
     theme(strip.text = element_text(size = 12)) +
     labs(x = "GC (%)", y = "Density of SVs") +
     scale_x_continuous(breaks = seq(0, 100, 20), labels = seq(0, 100, 20)) 

ggsave(oup,height = 3, width= 12)


