#Done: 1,2,3,4,5,6, 10, 8, 11
#Need to do: 7
#Non-programming tasks: 9,12 

#1 Harness the power of the Internet

#2 Read table
d = read.table("data/HNSC.txt", header = T, row.names = 1, stringsAsFactors = F)
dim(d)
d = as.matrix(d)
all_tumors = grep("Tumor",colnames(d),ignore.case = T)
all_controls = grep("Control",colnames(d),ignore.case = T)

#3 Heatmap of all data
library(pheatmap)
d = d[rowSums(d) > 0,]
col.features = data.frame(Type = factor(rep(c("Tumor","Control"), each=length(all_tumors))))
rownames(col.features) = colnames(d)
pheatmap(log2(d + 1),
         cluster_cols = F,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("darkblue", "white", "red2")) (256),
         border_color = NA,
         show_rownames = F,
         show_colnames = F,
         xlab="Sample",
         main="Heatmap of HNSC",
         gaps_col=length(all_tumors),
         annotation_col=col.features,
         annotation_colors=list(Type = c(Tumor="red",Control="green")))

#4 Calculate mean, sd, and ttest p value
means = vector()
sds = vector()
pvals = vector()
for(row in rownames(d)) {
  means = c(means, mean(d[row,]))
  sds = c(sds, sd(d[row,]))
  pvals = c(pvals, t.test((d[row,all_tumors]), (d[row,all_controls]))$p.value)
}
stats = cbind(means, sds, pvals)
rownames(stats) = rownames(d)
colnames(stats) = c("Mean", "Standard Deviation", "P-Value")

#5 Find 30 most significant miRNA names
ordered.pvals = order(pvals)
sorted.pvals = pvals[ordered.pvals]
rows.top30 = rownames(d)[ordered.pvals[1:30]]
d.top30 = d[rows.top30,]
d.top30 = d.top30[rowSums(d.top30) > 0,]
d.sig = d[which(pvals > 0.05),]

#6 Make scatterplot
plot(d["hsa-miR-200a-3p",], d["hsa-miR-429",],
     col=c(rep("red", length(all_tumors)), rep("blue", length(all_tumors))),
     main="Sample Expression for hsa-miR-200a-3p and hsa-miR-429",
     xlab="hsa-miR-200a-3p",
     ylab="hsa-miR-429")
abline(lm(d["hsa-miR-429",]~d["hsa-miR-200a-3p",]))
cor(d["hsa-miR-200a-3p",], d["hsa-miR-429",])

#7 Heatmap for top 30
col.features = data.frame(Type = factor(rep(c("Tumor","Control"), each=length(all_tumors))))
rownames(col.features) = colnames(d.top30)
pheatmap(log2(d.top30 + 1),
         cluster_cols = F,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("darkblue", "white", "red2")) (256),
         border_color = NA,
         show_rownames = F,
         show_colnames = F,
         xlab="Sample",
         main="Heatmap of HNSC",
         gaps_col=length(all_tumors),
         annotation_col=col.features,
         annotation_colors=list(Type = c(Tumor="red",Control="green")))

#8 correlation
#a correlation between top 30 and 29 other, all samples
top30rownames = rownames(d)[ordered.pvals[1:30]]
corAll = cor(d[top30rownames,],d[31:60,])
#b correlation among top 30, all samples, control, tumor
corTC = cor(d[top30rownames,all_tumors],d[top30rownames,all_controls])
corAT = cor(d[top30rownames,],d[top30rownames,all_tumors])
corAC = cor(d[top30rownames,],d[top30rownames,all_controls])

#9 patterns between mirnas



#10 save mirna names significantly up/down regulated
d.up = data.frame()
d.down = data.frame()
for(i in 1:dim(d.top30)[1]) {
  d.tum = mean(d.top30[i,all_tumors], na.rm=T)
  d.con = mean(d.top30[i,all_controls], na.rm=T)
  d.ratio = d.tum / d.con
  if(!is.nan(d.ratio)) {
    if(d.ratio > 1) {
      d.up = rbind(d.up, d.top30[i,])
      rownames(d.up)[dim(d.up)[1]] = rownames(d.top30)[i]
    } else {
      d.down = rbind(d.down, d.top30[i,])
      rownames(d.down)[dim(d.down)[1]] = rownames(d.top30)[i]
    }
  }
}
colnames(d.up) = colnames(d.top30)
colnames(d.down) = colnames(d.top30)

#11
target = read.table("data/h_mirtb_strongevid_target.txt", stringsAsFactors = F)
target = as.matrix(target)

# Get target genes for upregulated
target.up = vector()
for(rna in rownames(d.up)) {
  for(row in which(target[,1] == rna, arr.ind=T)) {
    target.up = c(target.up, target[row,2])
  }
}

# Get target genes for downregulated
target.down = vector()
for(rna in rownames(d.down)) {
  for(row in which(target[,1] == rna, arr.ind=T)) {
    target.down = c(target.down, target[row,2])
  }
}

#12
write(target.up, "target_up.txt")
write(target.down, "target_down.txt")
# https://string-db.org/cgi/input.pl