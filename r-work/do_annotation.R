source("func.R")
# External input sources
files <- commandArgs(T)
mature.mirnas.pat.file <- defArgs(files, "../miRNAs/mature_miRNAs.pat")
files <- files[-1]
mirna.regions.file <- defArgs(files, "../miRNAs/mature_miRNA_regions.csv")
files <- files[-1]
trna.bed.file <- defArgs(files, "../mapped/larval_trnas_mapped.bed")
files <- files[-1]
rfam.annotation.file <- defArgs(files, "../annotation/rfam/larval-rfam_e1.csv")

# load data from previous scripts (do_qc.R)
load(".RData")

# Combine annotation mappings with counts
read.dt <- genome.counts[,list(read)]

# mirs from mir regions
mir.reads <- data.table(read.table(mature.mirnas.pat.file,col.names=c("mirid", "read", "start", "end", "strand", "mismatches"), stringsAsFactors=F))
mir.reads <- mir.reads[strand == "+"]
mirna.regions <- data.table(read.csv(mirna.regions.file, stringsAsFactors=F))
mir.reads[, mirname := gsub("(.+)-\\dp", "\\1", mirid)]
mir.reads[, arm := gsub(".+-(\\d)p", "\\1_prime", mirid)]
setkey(mir.reads, mirname, arm)
setkey(mirna.regions, mirname, arm)
mir.reads <- mirna.regions[mir.reads][,list(mirid, read, source, mir.region = region.seq)]

mir.reads.unique <- mir.reads[order(-source)][!duplicated(read)]
setkey(mir.reads.unique, read)
read.annot <- mir.reads.unique[read.dt][,list(read,mirid, mir.region)]
setnames(read.annot, "mirid", "miRNA")
read.annot[!is.na(miRNA), type := "miRNA"]

trnas <- data.table(read.table(trna.bed.file, col.names=c("trna", "start","end","read","count", "strand"), stringsAsFactors=F))
trnas[,count := NULL]
trnas.unique <- trnas[!duplicated(read)]
setkey(trnas.unique, read)
read.annot2 <- trnas.unique[,list(read,trna)][read.annot]
setnames(read.annot2, "trna", "tRNAScan")
read.annot2[!is.na(tRNAScan), type := "tRNA"]

rfam <- as.data.table(read.csv(rfam.annotation.file, stringsAsFactors=F))
rfam.unique <- rfam[order(mismatches)][!duplicated(seq)]

rfam.names <- rfam.unique[,.N,by=hit][order(-N)]
rfam.names[grepl("mir|MIR|let-7|bantam",hit), type:="miRNA"]
rfam.names[grepl("rRNA",hit), type:="rRNA"]
rfam.names[grepl("SNO",hit), type:="snoRNA"]
rfam.names[grepl("7SK|^U\\d$",hit), type:="snRNA"]
rfam.names[grepl("tRNA",hit), type:="tRNA"]
rfam.names[grepl("SRP|6S",hit), type:="SRP"]
rfam.names[is.na(type), type := "other"]

setkey(rfam.names, hit)
setkey(rfam.unique,hit)
rfam.unique <- rfam.unique[rfam.names[,list(hit, type)]]

setkey(rfam.unique, seq)
read.annot3 <- rfam.unique[,list(seq,hit,type)][read.annot2]
read.annot3[is.na(i.type), i.type:=type]
read.annot3[,type:=NULL]
setnames(read.annot3, "i.type", "type")
setnames(read.annot3, c("seq", "hit"), c("read", "Rfam"))

read.annot3[is.na(tRNAScan) & is.na(miRNA) & is.na(Rfam), unannotated:="unannotated"]
read.annot3[!is.na(unannotated), type := unannotated]
 
# Condense miRNAs in mir regions
read.annot3[is.na(mir.region), mir.region := read]
annotation.sources <- c("tRNAScan", "miRNA", "Rfam", "unannotated")
read.annotm <- melt(read.annot3, measure.vars=annotation.sources, variable.name="source", value.name="name")
read.annotm <- read.annotm[!is.na(name)]
read.annotm[,source := factor(source, levels=annotation.sources)]
read.annotm <- read.annotm[order(source)][!duplicated(read)]

read.region.unique <- read.annotm[order(source)][!duplicated(mir.region)]

setkey(read.annotm,read)

annotation.types <- unique(read.annotm[,type])
annotation.types.short <- c(  miRNA="miRNA", tRNA="t/rRNA", rRNA="t/rRNA", other="other ncRNA", "snRNA" = "other ncRNA", "snRNA"="other ncRNA", 
    SRP="other ncRNA", "snoRNA"="other ncRNA", unannotated="unannotated")
read.annotm[,stype:=annotation.types.short[type]]
read.annotm[,stype:=factor(stype, levels=unique(annotation.types.short))]

# Mapping statistics and bar plot
annotated.mapping.stats <- genome.counts[read.annotm][,lapply(.SD, sum),.SDcols=bee.samplenames,by=type]
total <- bees.counts[,lapply(.SD, sum), .SDcols=bee.samplenames]
total.unmapped <- total-annotated.mapping.stats[,lapply(.SD, sum),.SDcols=bee.samplenames]
annotated.mapping.stats <- rbind(annotated.mapping.stats, cbind(type="unmapped", total.unmapped))
thesis.mapping.stats <- rbind(annotated.mapping.stats[type %in% c("unannotated", "unmapped", "miRNA")], 
                              cbind(type="other", annotated.mapping.stats[type %in% c("other ncRNA", "snRNA", "tRNA", "rRNA", "SRP", "snoRNA"), 
                                    lapply(.SD, sum), .SDcols=bee.samplenames]))
thesis.mapping.stats <- melt(thesis.mapping.stats, measure.vars=bee.samplenames, variable.name="Sample", value.name="count")
thesis.mapping.stats[,Sample := factor(Sample, levels=rev(bee.samplenames))]

g <- ggplot(thesis.mapping.stats, aes(Sample, count, fill=type)) + 
geom_bar(position="stack", stat="identity") +
    coord_flip() + labs(y="Redundant count", fill="Annotation") + theme(axis.text=element_text(size=rel(1.5)), title=element_text(size=rel(1.6)), 
                                                                    legend.text=element_text(size=rel(1.5)))
ggsave2("annotation_mapped_bars.pdf", g, 6, 6) 

annotated.mapping.stats.scd <- genome.counts[read.annotm][,lapply(.SD, sum),.SDcols=bee.samplenames,by=list(type, size=nchar(read))]
annot.mapstats.scd.m <- melt(annotated.mapping.stats.scd, measure.vars=bee.samplenames, variable.name="Sample", value.name="count")
annot.mapstats.scd.m[,treatment:=gsub("(\\w\\w).\\d", "\\1", Sample)]
annot.mapstats.scd.m[,replicate:=gsub("\\w\\w.(\\d)", "\\1", Sample)]
g <- ggplot(annot.mapstats.scd.m, 
    aes(factor(size), count, group=Sample, colour=treatment)) + 
geom_line() + labs(y="Redundant count") + facet_wrap(~ type, ncol=2, scales="free_y")
ggsave2("annotation_scd.png", g, w=10, h=7)

# save for future scripts
save.image()
