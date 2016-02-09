
# functions relating to finding offsets based
# KL divergence measures
find_strand_bias <- function(bedcounts)
{
    # bedcounts: bed table with a count column and window id
    # must have columns (seqid, window, count, strand)

    ## window abundance totals
    lib.by.window <- bedcounts[, list(total = sum(count)), keyby=list(seqid, window)]

    ## Find strand bias of each window
    lib.by.window.strand <- bedcounts[,list(strand.total = sum(count)), keyby=list(seqid, window, strand)]
    lib.by.window.strand <- lib.by.window[lib.by.window.strand]
    lib.by.window.strand[, strand.proportion := abs(0.5 - strand.total/total)]
    lib.by.window <- lib.by.window.strand[, list(strand.bias=ifelse(length(strand.proportion) > 1, 
                                                           sum(strand.proportion), # case where there is another strand presense on this window
                                                           strand.proportion + 0.5) # case where there is no other strand presense - sb = 1
                                                ), keyby=list(seqid, window, total)]
    return(lib.by.window)
}

calc_kl_sb <- function(sb.out, N)
{
    # calculate KL divergence based on strand bias variations
    # in a sliding sequence window
    # sb.out: result from find_strand_bias

    # Find N bins of strand bias
    sb.out[, sb.bin := cut(strand.bias, N)]

    # number of windows in each bin
    lib.by.sbBin <- sb.out[,list(num.windows=length(window)), keyby=sb.bin]

    # number of windows for each total
    lib.by.total <- sb.out[,list(windows.per.total=length(window)), keyby=total]

    # number of windows per abundance level, per bin
    lib.by.total.bin <- sb.out[,list(windows.per.bin.total=length(window)), keyby=list(total, sb.bin)]

    # caluclate kl divergence measure
    by.total.bin <- lib.by.total.bin[lib.by.total]
    by.total.bin[, pi := windows.per.bin.total/windows.per.total] 
    by.total.bin[, kl.by.bin := pi*log(pi/(1/N))]
    by.total <- by.total.bin[, list(kl = sum(kl.by.bin)), keyby=total]

    return(by.total)
}


findOffsets <- function(counts, bed, read.cols="read", window.length, N, max.window.abundance=500, 
    value=c("mean", "max", "median", "all"), lo=T, lo.span=0.3)
{
    offsets <- data.table(sample=character(0), offset=numeric(0))
    if(value[[1]] == "kl")
        kl <- data.table(total=numeric(0), kl=numeric(0), lo=numeric(0), lib=numeric(0))
    
    for(i in colnames(counts[, !read.cols, with=F]))
    {
        # extract sample
        lib <- counts[,c("read",i), with=F]

        # join with bed columns
        setkey(bed, read)
        lib.bed <- bed[lib, nomatch=0]

        # rename count col
        setnames(lib.bed, i, "count")
        lib.bed <- lib.bed[count > 0]

        if(value[[1]] == "kl"){
            thiskl <- find_offset_for_sample(lib.bed, window.length, N, max.window.abundance, kl=T)
            kl <- rbind(kl, cbind(thiskl,lib=i))
        }
        else{
            # find offset
            of <- find_offset_for_sample(lib.bed, window.length, N, max.window.abundance, lo=T, lo.span=lo.span)

            # find abundance at global minimum of kl divergence
            offsets <- rbind(offsets, list(sample=i, offset=of))
        }
        
   }

   if(value[[1]] == "kl")
     return(kl)

   if(value[[1]] == "all")
     return(offsets)
   else
     # any other value interpreted as a summary function on offsets
     return(do.call(value[[1]], list(offsets)))
}

find_offset_for_sample <- function(bedcounts, window.length, N, max.window.abundance=500, kl=F, lo=T, lo.span=0.2)
{
    # finds an offset for one sample
    # bedcounts: bed with counts
    ## Find a sliding window of reads, length window.length
    bedcounts[, window := findWindow(bedcounts, window.length)]

    # find strand bias per window
    lib.by.window <- find_strand_bias(bedcounts)[total <= 500]
    lib.by.window[,total := floor(total)]

    # kl divergence on strand bias
    kl.by.total <- calc_kl_sb(lib.by.window, N) 

    # Find smoothed value
    kl.by.total[,lo := predict(loess(kl ~ total, .SD, span=lo.span))]
        
    # return offset at minimum
    if(kl)
        return( kl.by.total[,list(total,kl,lo)])

    varmin <- "kl"
    if(lo)
        varmin <- "lo"
    return( kl.by.total[which.min(kl.by.total[[varmin]]), total] )
}

findWindow <- function(mapping, winlength)
{
    halfway <- (mapping$start + mapping$end)/2
    win <- floor(halfway / winlength)
    return(win)
}

test_offsets <- function(counts.norm, offsets=c(1, seq(10,1000,10)))
{
    # test the effect of different offsets on library comparison
    # counts.norm: data.table of read column and two normalised libraries
    ma <- data.table(read=character(0), M=numeric(0), A=numeric(0), count.mean=numeric(0), offset=numeric(0))
    sample.names <- colnames(counts.norm[,!"read",with=F])

    for(o in offsets)
    {
        ma <- rbind(ma, cbind(read=counts.norm[,read], getMA(counts.norm[,!"read",with=F], base=2, m.offset=o, a.offset=1), 
                              count.mean=counts.norm[,apply(.SD,1,mean),.SDcols=sample.names], offset=o))
    }
    return(ma)
}

test_sig_offset_curve <- function(counts.norm, offsets=1:100)
{
    ma <- data.table(sig=numeric(0), offset=numeric(0))
    for(o in offsets)
    {
        m <- getM(counts.norm[,!"read",with=F], b=2, o=o)
        ma <- rbind(ma, data.table(sig=sum(abs(m)>1), offset=o))
    }
    return(ma)
}

test_sbKL_params <- function(bedcounts, max.window.abundance=500, window.sizes=seq(40,500,10), Ns=10, value=c("min", "KL") ) 
{
    # Experiment on the effect of window size and bin size on the minimum abundance
    # from KL divergence

    # value: min return minimum abundance for each experiment
    #         KL return KL distribution for each experiment
    paramvar <- data.table(total=numeric(0), kl=numeric(0), window.length=numeric(0), N=numeric(0), lo=numeric(0))
    for(window.size in window.sizes)
    {
        msg <- paste0("Window size: ", window.size)
        for(N in Ns)
        {
            print(paste0(msg, " Bin: ", N))
            paramvar <- rbind(paramvar, cbind(find_offset_for_sample(bedcounts, window.size, N, max.window.abundance,  kl=T),
                                              window.length=window.size, N=N))
        }
    }
    if(value[1] == "min")
        return (paramvar[, list(offset=total[which.min(kl)]), by=list(window.length, N)])
    else
        return (paramvar)
    return(paramvar)
}

test_sb_windowsize <- function(bedcounts, window.sizes=seq(40,500,10))
{
    dt <- data.table(strand.bias=numeric(0), seqid=character(0), window=numeric(0), total=numeric(0), window.size=numeric(0))
    for(winsize in window.sizes)
    {
        print(paste0("Window size: ",winsize))
        bedcounts[, window := findWindow(bedcounts, winsize)]
        windows <- find_strand_bias(bedcounts)
        dt <- rbind(dt, cbind(windows, window.size=winsize))
    }
    return(dt)
}

extract_offsets_from_kltest <- function(kltest){
    # discern offsets from a kltest where kl values are written out
    kltest.offsets <- kltest[,list(Smoothed=total[which.min(lo)], Unsmoothed=total[which.min(kl)]), by=list(window.length)]
    kltest.offsets.m <- melt(kltest.offsets, measure.vars=c("Smoothed", "Unsmoothed"))
    kltest.offsets.m[,variable := factor(variable, levels=c("Unsmoothed", "Smoothed"))]
}

gg_kltest <- function(kltest, klmin) {
        ggplot(kltest, aes(total, kl )) + 
    geom_line(aes(group=interaction(factor(N), factor(window.length)), colour=factor(window.length), linetype=factor(N))) +
    geom_vline(data=klmin, aes(xintercept=minKl, colour=factor(window.length)), linetype="dashed")
}

gg_kl <- function(kl, filename){
    g <- ggplot(kl, aes(x=total)) +
         geom_line(aes(y=kl), colour="grey") + geom_line(aes(y=lo), colour="blue") +
         geom_vline(aes(xintercept=total[which.min(lo)]), linetype="dashed") +
         labs(x="Expression Level", y="KL Divergence") + theme_nar
    ggsave(filename, g, w=3.3, h=2)
}

gg_kltest_lo <- function(kltest, klmin=NULL, continuous.colour.scale=F) {
        klmin <- kltest[,list(minKl=total[which.min(kl)], minLo=total[which.min(lo)]), by=list(window.length, N)]

    if(!continuous.colour.scale){
        kltest$window.length <- factor(kltest$window.length, levels=kltest[order(window.length), unique(window.length)])
        klmin$window.length <- factor(klmin$window.length, levels=klmin[order(window.length), unique(window.length)])
    }
    ggplot(kltest, aes(total, lo, group=window.length,)) + 
    geom_point(aes(y=kl, colour=window.length), size=0.5)  +
    geom_line( aes(colour=window.length) ) +
    geom_vline(data=klmin, aes(xintercept=minLo, colour=window.length), linetype="dashed") + 
    labs(colour="window length", x="Abundance level", y="KL divergence")# +
    #scale_colour_gradientn(colours=rainbow(5))
}

gg_kltest_offsets <- function(kltest.offsets){
    ggplot(kltest.offsets, 
        aes(window.length, value, colour=variable, group=variable)) + geom_line(size=0.3) +
        scale_colour_manual(values=c("grey", "blue"), guide=F) +
        labs(colour="Curve", x="Window length", y="Derived offset") + scale_x_continuous(labels=comma) 
}

run_kltest <- function(bed, counts, sample.names, normalisations, win.lengths=seq(50,6000,50), plotted.windows=c(100,300,500,1000,3000,6000), 
                       max.abundance=200, kldir="kltest" , returnOffset=TRUE){
    dir.create(kldir)
    setwd(kldir)
    tests <- list()
    kl <- list()
    for(thissample in sample.names){
        normtest <- list()
        normkl <- list()
        for(norm in normalisations){
            print(c(thissample, norm));
            bedcounts <- bed[counts[method %in% norm][,c("read", thissample), with=F]]
            setnames(bedcounts, thissample, "count")
            bedcounts <- bedcounts[count>0]
            kltest <- test_sbKL_params(bedcounts, max.abundance, win.lengths, 100, value=F)
            offsets <- extract_offsets_from_kltest(kltest)

            normtest[[norm]] <- offsets
            normkl[[norm]] <- kltest

            # plot offsets by window length
            g <- gg_kltest_offsets(offsets) + theme_nar + scale_x_continuous(breaks=seq(1000,10000,1000))
            ggsave(paste0("kltest_",norm,"_",thissample,"_offsets.pdf"), g, w=3.3, h=2)

            # plot selected curves
            g <- gg_kltest_lo(kltest[window.length %in% plotted.windows], continuous.colour.scale=F) + theme_nar +
                    guides(colour=guide_legend(ncol=3))
            ggsave(paste0("klLoCurve_",norm,"_",thissample,".png"), g, w=nar.col.width, h=nar.col.width*0.8, units="mm", dpi=600)
        }
        tests[[thissample]] <- spBind(spList=normtest, libname="norm")
        kl[[thissample]] <- spBind(spList=normkl, libname="norm")
    }

    if(returnOffset)
        return(spBind(spList=tests, libname="sample"))
    else
        return(spBind(spList=kl, libname="sample"))

    setwd("../")
}

# detach functions and then reattach to update
while("f" %in% search())
    detach("f")
attach(f)
