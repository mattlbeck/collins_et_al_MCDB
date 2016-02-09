source("func.R")
source("norm_func.R")
library(plyr)

# External data to load. Edit the default file paths if loading from an R environment
files <- commandArgs(T)
mapped.file <- defArgs(files, "../mapped/mapped.bed")
files <- files[-1]
combined.counts.file <- defArgs(Files, "../preprocess/combined.csv")
files <- files[-1]


# declare various versions of sample names
bees.libnames <- paste("S", 1:16, sep="")
bee.groupnames <- c("EW","LW","EQ","LQ")
bee.groupnames.full <- c("Early Worker Destined", "Late Worker Destined", "Early Queen Destined", "Late Queen Destined")
bee.samplenames <- paste(rep(bee.groupnames, each=4), rep(1:4, times=4), sep=".")
bee.labels <- c("Early Worker", "Late Worker", "Early Queen", "Late Queen")

t.pat <- c("(\\w+).\\d")
r.pat <- c("\\w+.(\\d)")

# load mapping file
bees.bed <- as.data.table(readBed(mapped.file))
bees.bed[,count:=NULL]

# load all counts
bees.counts <- as.data.table(read.csv(combined.counts.file), header=T, stringsAsFactors=F)

# original column names if the run.sh script was used to assign 1:16 as sample names
setnames(bees.counts, paste("X", 1:16, sep=""), bee.samplenames)

setkey(bees.bed, read)
setkey(bees.counts, read)
genome.counts <- bees.counts[read %in% unique(bees.bed[,read]),]

mapping.table.r <- mapping.summary(bees.counts, genome.counts)
mapping.table.nr <- mapping.summary(bees.counts, genome.counts, r=F)

# MA replicates
counts.m <- melt_counts(genome.counts, bee.samplenames)

counts.m[, treatment := str_extract(sample, "\\w\\w")]
counts.m[, replicate := paste0("R",gsub("\\w\\w.(\\d)", "\\1", sample))]
rep.names <- c("R1", "R2", "R3", "R4")
counts.r <- dcast.data.table(counts.m, read + treatment ~ replicate, value.var="count")
counts.r <- cleancounts(counts.r, count.cols=rep.names)
counts.r[, treatment := factor(treatment, levels=bee.groupnames)]

repc.names <- c();
for(i in 1:length(rep.names)){
    for(j in  1:length(rep.names)){
        if(j!=i)
       repc.names <- c(repc.names, paste(rep.names[c(i,j)], collapse=",")) 
    }
}

counts.repc <- spBind("R1,R2"=counts.r[,list(read,treatment,ref=R1,obs=R2)], 
                      "R1,R3"=counts.r[,list(read,treatment,ref=R1,obs=R3)], 
                      "R1,R4"=counts.r[,list(read,treatment,ref=R1,obs=R4)], 
                      "R2,R3"=counts.r[,list(read,treatment,ref=R2,obs=R3)], 
                      "R2,R4"=counts.r[,list(read,treatment,ref=R2,obs=R4)], 
                      "R3,R4"=counts.r[,list(read,treatment,ref=R3,obs=R4)],
                      libname="Comparison")

counts.repma <- data.table(counts.repc, getMA(counts.repc[,list(ref,obs)], base=2, offset=1))

g <- ggplot(counts.repma, aes(A,M)) + 
    geom_point(size=0.4, alpha=0.4) + labs(x=expression("A ("*log[2]*" average abundance)"), y=expression("M ("*log[2]*" fold change)")) +
    facet_grid(treatment~Comparison) +
    geom_hline(yintercept=0) + 
    theme(axis.text=element_text(size=rel(1.25)), strip.text=element_text(size=rel(1.25)), axis.title=element_text(size=rel(1.25)))
ggsave2("paper_maqc.png", g, 12, 8)

genome.jaccard <- dt_jaccard_table(genome.counts,threshold=500)
jbeenames <- gsub("\\.", "", bee.samplenames)
genome.jaccard[,Var1 := factor(gsub("\\.","",Var1), levels=jbeenames)]
genome.jaccard[,Var2 := factor(gsub("\\.","",Var2), levels=jbeenames)]
g <- dt_jaccard_heat_map(genome.jaccard,  xrot=T, grouping.lines=c(4,8,12), diagx=T) + coord_fixed()
ggsave2("paper_jaccard.pdf", g, 10,10)


# Quality check filter: removal of poor quality replicates and squeezing of the size class range (if also low quality)
# removing two replicates
# 4: erroneous size class with no miRNA peak
# 14: very high miRNA peak and erroneous fold change at igher size class end
rep.remove <- c(8,14) 
bee.samplenames.rm <- bee.samplenames[-rep.remove]

# size class squeeze. Keep 16-27 but ignore longer reads because they are causing slight fold change issues
sc.keep <- 16:27

counts.rm <- genome.counts[, c("read",bee.samplenames.rm),with=F][nchar(read) < 28]

# normalise using the Bootstrap method
counts.norm.all <- all_norm(counts.rm, c("btsp", "qnorm"))
counts.norm.btsp <- counts.norm.all[method=="btsp"]
counts.norm <- cleancounts(counts.norm.btsp, bee.samplenames.rm)
counts.norm.m <- melt_counts(counts.norm, bee.samplenames.rm)
#counts.norm[, total := apply(.SD, 1, sum), .SDcols=bee.samplenames.rm]

sizecounts.r <- counts.norm.m[,list(count=sum(count)),by=list(size=nchar(read),sample)] 
sizecounts.r[,c("replicate", "treatment") := list(gsub(r.pat, "\\1", sample), gsub(t.pat, "\\1", sample))]
sizecounts.mean <- sizecounts.r[,list(mean=mean(count)),by=list(size,treatment)]
sizecounts.sd <- sizecounts.r[,list(sd=sd(count)),by=list(size,treatment)]
sizecounts.bardt <- data.table(sizecounts.mean, sd=sizecounts.sd[,sd])
sizecounts.bardt[, title:=factor(treatment, levels=bee.groupnames, labels=paste0(letters[1:4], ") ",bee.groupnames))]
sizebp <- function(data){
    ggplot(data, aes(factor(size), mean)) + geom_bar(stat="identity", colour="black", fill="grey") + 
        facet_wrap(~ title, ncol=2, scales="free") + 
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5) +
        scale_y_continuous(labels=comma, expand=c(0,0)) + labs(x="Read length (nt)", y="Normalized count") + 
        theme(panel.border=element_blank(), panel.grid=element_blank(), strip.background=element_blank(), 
        axis.line=element_line(), panel.margin=unit(0.25, "in"),
        axis.text=element_text(size=rel(1.25)), axis.title=element_text(size=rel(1.25)), strip.text=element_text(size=rel(1.25), hjust=0))
}
bplist <- list()
for(i in unique(sizecounts.bardt[,title])){
    bplist[[i]] <- sizebp(sizecounts.bardt[title==i])
}
pdf("paper_sizeclass_bar.pdf", width=9, height=6.5)
multiplot(plotlist=bplist, layout=matrix(c(1,3,2,4), nrow=2))
dev.off()
# ggsave2("paper_sizeclass_bar.pdf", g, 10,10)

# write out counts
docsv(counts.rm, "all_counts.csv")

# Save data for future R scripts
save.image()
