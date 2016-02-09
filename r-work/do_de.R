source("func.R")
source("KL_func.R")

# load data from previous R scrips (do_annotation.R)
load(".RData")

### CONDENSE miRNAs ###
setkey(counts.norm, read)
setkey(read.annotm, read)
counts.annot <- read.annotm[counts.norm,nomatch=0]
counts.annot[,mir.region := toupper(mir.region)]
counts.annot2 <- counts.annot[, lapply(.SD, sum), .SDcols=c( bee.samplenames.rm), by=list(read=mir.region) ]

cmp.list <- list(c("EW", "LW"), c("EQ", "LQ"), c("EW", "EQ"), c("LW", "LQ"));
comparisons <- sapply(cmp.list, paste,collapse=".")

#### DE CALCULATION #####
de <- direction.combo(counts.annot2[,c("read", bee.samplenames.rm), with=F], 
    cmp.list,
    bee.samplenames.rm, "(\\w\\w).\\d", "\\w\\w.(\\d)", offsets=data.table(sample=bee.samplenames.rm, offset=20))

setkey(de$fc, read)
setkey(read.region.unique, mir.region)
read.region.unique[,read := NULL]
read.region.unique[,mir.region := toupper(mir.region)]
defc.annot <- read.region.unique[de$fc,nomatch=0][,list(read, stype, source, name, comparison, d, ofc, ref, obs, A)]

# Cross plots for combinations of comparisons
g <- de_xplot(defc.annot, "EW.LW", "EQ.LQ") + labs(x="EW to LW", y="EQ to LQ", colour="Annotation")
ggsave2("de_xplot_EW.LW-EQ.LQ.png", g, w=5, h=5)

g <- de_xplot(read.annotm[de$fc], "EW.EQ", "LW.LQ") + labs(x="EW to EQ", y="LW to LQ", colour="Annotation")
ggsave2("de_xplot_EW.EQ-LW.LQ.png", g, w=5, h=5)

g <- de_xplot(read.annotm[de$fc], "EW.LW", "LW.LQ")+ labs(x="EW to LW", y="LW to LQ", colour="Annotation")
ggsave2("de_xplot_EW.LW-LW.LQ.png", g, w=5, h=5)

g <- de_xplot(read.annotm[de$fc], "EQ.LQ", "EW.EQ")+ labs(x="EQ to LQ", y="EW to EQ", colour="Annotation")
ggsave2("de_xplot_EQ.LQ-EW.EQ.png", g, w=5, h=5)

defc.mirs <- defc.annot[source=="miRNA"]
defc.mirs.d <- dcast.data.table(defc.mirs, read + name ~ comparison, value.var="d")
defc.mirs.fc <- dcast.data.table(defc.mirs, read + name ~ comparison, value.var="ofc")
defc.mirs.fc2 <- defc.mirs.fc[defc.mirs.d[,apply(.SD, 1, function(x) sum(x != "S") > 1), .SDcols=comparisons]]
defc.mirs.fc.sig <- defc.mirs.fc2[defc.mirs.fc2[,apply(.SD, 1, function(x) sum(abs(x) > 1) > 0), .SDcols=comparisons]]
docsv(defc.mirs.fc.sig, "miRNAs_significant.csv")

# Preparation of signficant DE miRNAs
# Preparation of results tables
# - all predicted miRNAs with sifnificant fold changes on at least one comparison
# - counts of all predicted miRNAs

# All miRNAs where total > 100
setkey(read.region.unique, mir.region)
read.region.unique[,read := NULL]
read.region.unique[,mir.region := toupper(mir.region)]
setkey(counts.annot2, read)
counts.annot3 <- counts.annot2[read.region.unique,nomatch=0]
counts.mirs <- counts.annot3[source=="miRNA"]
counts.mirs[,total := apply(.SD, 1, sum), .SDcols=bee.samplenames.rm]

# add top mature sequence for each region
mirna.regions[, mirname:=paste(mirname, gsub("(\\d).+", "\\1p", arm), sep="-")]
setkey(mirna.regions, region.seq)
setkey(counts.mirs,read)
counts.mirs2 <- mirna.regions[,
                 list(mirname, region.seq, top.mature, precursor.seq, pre.start=start, pre.end=end, seqid, strand,source)
                 ][counts.mirs[,c("read", bee.samplenames.rm, "total"), with=F]]

source.names <- c("MC", "MD", "MC/MD", "MM", "MC/MM", "MD/MM", "MC/MD/MM")
counts.mirs2[,source.name := source.names[source]]
mirs.out <- counts.mirs2[total > 100][,c("miRNA"="mirname", "mature.region"="region.seq", "precursor"="precursor.seq", 
                                        "seqid", "strand", "pre.start", "pre.end", source="source.name", bee.samplenames.rm),with=F]
setnames(mirs.out, c("name", "mature sequence", "precursor sequence", "seqid", "strand", "precursor start", "precursor end", "source tools", bee.samplenames.rm))
mirs.dups <- mirs.out[, .N, by=region.seq]
mirs.dups[N>6]
docsv(mirs.out, "all_miRNAs.csv")

# DE bar plots of Northern blotted miRNAs
noblo.mirs <- c("miR-6001-5p", "miR-6001-3p", "miR-11-3p", 
                "miR-12-5p", "miR-9a-5p", "miR-71-5p", 
                "miR-184-3p", "miR-275-3p", "miR-283-5p", "miR-315-5p")
noblo.mirs.formal <- paste(letters[1:10], paste("bte-",c("miR-6001-5p", "miR-6001-3p", "miR-11", 
                "miR-12", "miR-9a", "miR-71", 
                "miR-184", "miR-275", "miR-283", "miR-315"), sep=""), sep=") ")

mirs.tobp <- mirs.out[name %in% noblo.mirs]
mirs.tobp[,name := factor(name, levels=noblo.mirs)]
mirs.tobp <- mirs.tobp[order(name)]
mirs.tobp[,formal := noblo.mirs.formal]
for(m in noblo.mirs){
    g <- plot_de2(mirs.out[name==m], bee.samplenames.rm, bee.groupnames, bee.groupnames)
    pdf(paste0(m,".pdf"),width=4, height=4)
    print(g + labs(title=paste(m, mirs.out[name==m, read], sep="\n", x="Phenotype")) 
          + theme(plot.title=element_text(size=10))
          )
    dev.off()
}

# all in one handy pdf
g <- plot_de2(mirs.tobp[name%in%noblo.mirs], 
              bee.samplenames.rm, bee.groupnames, bee.groupnames, readcol="formal")

pdf(paste0("noblo_mirs.pdf"),width=6, height=12)
print(g + facet_wrap(~ formal, scale="free", ncol=2)
      + labs(x="Phenotype") 
      + theme(plot.title=element_text(size=10), strip.background=element_blank(), strip.text=element_text(hjust=0, size=rel(1.2)))
      )
dev.off()
