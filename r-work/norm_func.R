library(ggplot2)
library(edgeR)
library(data.table)
#source("/scratch/mbeckers/mouse/ro60/r-work/func.R")

# TMM functions
source("~/scripts/R/tmm.R")

# QNORM functions
source("~/scripts/R/qnorm.R")

# bootstrap functions
source("~/scripts/R/bootstrap.R")

# deseq functions for general matrixes (functions stripped of DEseq lib ecosystem)
source("~/scripts/R/deseq.R")

f$nzcounts <- function(counts, count.cols){ 
    # Remove sequences in a count matrix that are zero in 1 or more samples
    # args:
    #  counts: data.frame of counts with columns "read", sample counts
    nzseqs <- apply(as.data.table(counts)[,count.cols,with=F],1,function(x) sum(x == 0) == 0)
    nzcounts <- counts[nzseqs,]
    return(nzcounts)
}

f$nzrow <- function(counts, count.cols){
    apply(as.data.table(counts)[,count.cols,with=F],1,function(x) sum(x == 0) == 0)
}

f$cleancounts <- function(counts, count.cols=colnames(counts)[-1]){
    # Remove count rows that contain all zeros. Useful after bootstrapping or splitting
    # an expression matrix
    nonzeros <- apply(as.data.table(counts)[,count.cols,with=F],1, function(x) sum(x != 0) > 0)
    cleaned <- counts[nonzeros,]
    return(cleaned)
}


f$all_norm <- function(raw.counts, norm.methods=c("tc", "btsp", "tmm", "qnorm2", "uq", "med", "deseq"), readcol="read", othercols=c(), ...){
    # Normalise using given normalisation methods
    # args:
    #   matrix of raw counts and read column

    # need to start with a data.table of doubles
    reads <- raw.counts[,c(readcol, othercols),with=F]
    raw.counts <- raw.counts[,!c(readcol,othercols),with=F]
    raw.counts <- as.data.table(lapply(raw.counts, as.numeric))
    norm.counts <- cbind(reads, raw.counts, method="raw")
    
    if("rpm" %in% norm.methods)
        stop("Please change tpm to tc from now on")
    if("tc" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, ...), method="tc"))
    if("tmm" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,tmm_norm_counts(raw.counts, tmm_norm_factors(raw.counts), ...), method="tmm")) 
    if("med" %in% norm.methods)
        #norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) median(x[x > 0]), upscale=1), method="med"))
        norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) persum(x, 0.5)), method="med"))
    if("logmed" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) median(log10(x[x > 0]))), method="logmed"))
    if("uq" %in% norm.methods)
        #norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) quantile(x[x > 0], 0.75), upscale=1), method="uq"))
        norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) persum(x, 0.75)), method="uq"))
    if("loguq" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,scale_norm_counts(raw.counts, FUN=function(x) quantile(log10(x[x > 0]), 0.75)), method="loguq"))
    if("qnorm" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,dt_q_norm_counts(raw.counts, keep.zeros=F, normalise.equal.values=F), method="qnorm"))
    if("qnorm2" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads,dt_q_norm_counts(raw.counts), method="qnorm2"))
    if("btsp" %in% norm.methods){
        btsp.counts <- cbind(norm_bootstrap(cbind(reads[,readcol,with=F],raw.counts), readcol=readcol), method="btsp")
        btsp.counts <- cleancounts(btsp.counts, count.cols=colnames(raw.counts))
        if(length(othercols) > 0)
           btsp.counts <- cbind(btsp.counts,  reads[,othercols,with=F])
        norm.counts <- rbind(norm.counts, btsp.counts)
    }
    if("deseq" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads, scale_factor_counts(raw.counts, estimateSizeFactors(raw.counts)), method="deseq"))

    # deseq's reference set derived from geometric means over samples results in too many zero's from log(1). +1 pseudocount may help resolve this:w
    if("deseq2" %in% norm.methods)
        norm.counts <- rbind(norm.counts, cbind(reads, scale_factor_counts(raw.counts, estimateSizeFactors(raw.counts+1)), method="deseq2"))

    # order method factor by input vector
    norm.counts[, method := factor(method, levels=c("raw", norm.methods))]
    return(norm.counts) 
}
f$readnorm <- function(counts, norm.func, ..., read.cols="read")
{
    # Wrapper to replace read column after using a normalisation function
    data.table(counts[,read.cols,with=F], norm.func(counts[,!read.cols, with=F], ...))
}

f$norm_counts <- function(counts, read.colname="read"){
    norms <- apply(as.data.table(counts)[,!read.colname,with=F],2,function(x)10^6*as.numeric(x)/sum(as.numeric(x)));
    return(data.table(read=counts[,get(read.colname)], norms))
}

f$get_global_scale <- function(counts, global.scale="mean", modifiers = rep(1, ncol(counts)))
{
    if(is.numeric(global.scale))
        gs <- global.scale
    if(global.scale == "mean"){
        gs <- mean(counts[, apply(.SD, 2, sum)] * modifiers)
        print(paste0("Using a global upscaling factor of mean of libraries: ",gs))
    }
    return(gs)
}

f$scale_norm_counts <- function(counts, FUN=sum, global.scale="mean"){
    gs <- get_global_scale(counts, global.scale)

    # scale counts using a global constant scale gs and a library scale FUN(x)
    as.data.table(apply(counts, 2, function(x) x * ( gs / FUN(x) ) ))
}

persum <- function(x, p, nz=T)
{
    if(nz)
        x <- x[x>0]
    sumfrom <- quantile(x, p)
    sum(x[x >= sumfrom])
}

f$scale_factor_counts <- function(counts, factors){
    counts <- copy(counts)
    for(i in seq_along(counts))
        counts[, i := counts[[i]]/factors[i], with=F]
    return(counts)
}

## Normalisation Measuring methods
calc_covariance <- function(norm.counts)
{
    apply(norm.counts,1,sd)/rowMeans(norm.counts)
}


calc_all_cov <- function(norm.counts, sample.combos=list(), readcols="read", nzc=F)
{
# Parameters are like so
    # (dt, list(control=1:2, treat=3:4))
    # nzc if only nonzero rows should be used
    #dt <- norm.counts[,readcols,with=F]
    #dt <- data.table(lib=character(0), covariance=double(0))
    dt <- c()
    for(i in names(sample.combos))
    {
        combo <- sample.combos[[i]]
        ncounts <- norm.counts
        if(nzc){
            #ncounts <- readnorm(norm.counts, nzcounts, count.cols=combo, read.cols=readcols)
            ncounts <- ncounts[nzrow(ncounts[,combo,with=F])]
        }
        reads <- ncounts[,readcols,with=F]
        #dt <- cbind(dt, calc_covariance(norm.counts[, combo, with=F]))
        print(head(ncounts[,combo,with=F]))
        dt <- rbind(dt, data.table(reads, lib=i, covariance=calc_covariance(ncounts[, combo, with=F])))
    }
#    dt <- data.table(dt);
#    setnames(dt, c(readcols,names(sample.combos)))
#    return(melt(dt, measure.vars=names(sample.combos), variable.name="lib", value.name="covariance"))
return(dt)
}

f$sample_counts <- function(reads, counts, size, replacement=F){
    seq.list <- sample(rep(reads, counts), size, replace=replacement)
    summary.df <- as.data.table(table(seq.list)) 
    setnames(summary.df, c("read", "count"))
    return(summary.df)
}

# XPLOTS - sig DE scatter between norms
gg_norm_xplot <- function(dd, nx, ny){
    # dd: dcast table of M values for each norm (column) and each read (row)
    dd[, sig := ifelse(abs(dd[[nx]])>1, ifelse(abs(dd[[ny]])>1, "both", nx), ifelse(abs(dd[[ny]])>1, ny, "none"))]
    dd[, sig := factor(sig, levels=c(nx, ny, "both", "none"))]

    ggplot(dd, aes_string(nx, ny)) +
        geom_hline(yintercept=c(-1,1), linetype="dashed") +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=c(-1,1), linetype="dashed") +
        geom_vline(xintercetp=0) +
        geom_point(size=1, aes(colour=sig)) +
        scale_colour_manual(values=c("red", "blue", "black", "grey50"), drop=FALSE) 
}

f$multisampling <- function(readcounts, factors, bins, nowidth=F, replacement=F, pseudocount=F){
    # Args:
    #   counts: sequence count matrix
    #   count.col: column on matrix to use as counts
    #   factors: list of sample sizes as a function of population size 
    #   bins: sequence of bin cuts to use 
    #   pseudocount: if TRUE, includes zero-rate sampling and changes factor calculation to total + 1 / sample + 1
    samples <- c()
    reads <- readcounts[[1]]
    counts <- readcounts[[2]]
    original <- data.table(read=reads, ocount=counts)
    total.counts <- sum(counts)
    for (f in factors){
        # find the sample size to use
        ssize <- round(total.counts*f)

        # sample from the population using size
        scounts <- sample_counts(reads, counts, ssize, replacement)
        print(total.counts)
        setkey(scounts, read)
        setkey(original, read)
        
        # add column for the population counts
        # Ignores seqs that were never sampled...
        smerge <- list()
        if (pseudocount){
            smerge <- merge(scounts, counts[counts[,count.col] > 1,count.col,drop=FALSE], by.x="read", by.y=0, all.y=TRUE)
            smerge[is.na(smerge)] <- 0
            smerge$factors <- (smerge[,3] + 1)/(smerge[,2] + 1)
            smerge$scale <- (smerge[,2] + 1)/(smerge[,3] + 1)

        }
        else{
            #smerge <- merge(scounts, counts[,count.col,with=F], by.x="read", by.y=0)
            smerge <- original[scounts] 
            # find the factor required to scale each count back to its original count
            smerge[, factors := count/ocount]

       }

        if(nowidth){
            smerge <- smerge[smerge[,3] %in% bins,] 
            smerge$bin <- factor(smerge[,3])
        }
        else{
            smerge$bin <- cut(smerge[[3]], c(0, bins, +Inf), c(bins, paste(tail(bins, 1)+1, "+", sep="")))
        }
        samples <- rbind(samples, cbind(smerge, resample=f))
    }
    return(samples)
}

# detach functions and then reattach to update
while("f" %in% search())
    detach("f")
attach(f)
