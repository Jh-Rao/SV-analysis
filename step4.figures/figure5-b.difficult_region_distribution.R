library(ggplot2)

args<-commandArgs(TRUE)

dfl <- args[1]
ifl <- args[2]
oup <- args[3]


dd <- read.table(dfl, head=T)
di <- read.table(ifl, head=T)

n = length(unique(dd$DifficultRegion))
m = n * 2 * 5

frame = data.frame(platform = rep(0, m),difftype = rep(0, m), svtype = rep(0, m), value = rep(0, m))

for (i in 1:n){
  frame$platform[i] = "DNBSEQ"
  frame$svtype[i] = "All"
  frame$difftype[i] = names(table(dd$DifficultRegion))[i]
  frame$value[i]  = as.numeric(table(dd$DifficultRegion))[i] / sum(as.numeric(table(dd$DifficultRegion))) * 100

  j=i+n
  frame$platform[j] = "Illumina"
  frame$svtype[j] = "All"
  frame$difftype[j] = names(table(di$DifficultRegion))[i]
  frame$value[j]  = as.numeric(table(di$DifficultRegion))[i] / sum(as.numeric(table(di$DifficultRegion))) * 100
}
for (k in 1:4){
  sdd = subset(dd, SVtype == unique(dd$SVtype)[k])
  sdi = subset(di, SVtype == unique(dd$SVtype)[k])

  for (i in 1:n){
    m = i + k * 2 * n
  frame$platform[m] = "DNBSEQ"
  frame$svtype[m] = unique(dd$SVtype)[k]
  frame$difftype[m] = names(table(sdd$DifficultRegion))[i]
  frame$value[m]  = as.numeric(table(sdd$DifficultRegion))[i] / sum(as.numeric(table(sdd$DifficultRegion))) * 100

    l = m + n
  frame$platform[l] = "Illumina"
  frame$svtype[l] = unique(dd$SVtype)[k]
  frame$difftype[l] = names(table(sdi$DifficultRegion))[i]
  frame$value[l]  = as.numeric(table(sdi$DifficultRegion))[i] / sum(as.numeric(table(sdi$DifficultRegion))) * 100
  }
}

frame = subset(frame, difftype == "Diff")

p <- ggplot(frame, aes(svtype, weight = value, fill = platform)) +
     geom_bar(width = .7, position = 'dodge') +
     labs(y = '% of SVs', x = '') +
     theme_classic() 
ggsave(oup, height = 2, width = 5)

