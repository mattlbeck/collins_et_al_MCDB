nar.col.width <- 84
nar.col.double.width <- 178
require(ggplot2)
require(reshape2)
require(scales)
require(data.table)
require(stringr)
theme_set(theme_bw())
theme_nar <- theme(legend.position="bottom", text=element_text(size=8), panel.background=element_rect(colour="black"), panel.grid=element_blank())
# Create enviroment for these functions
f <- new.env()

f$theme_min <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(), panel.grid=element_blank(), axis.ticks.x=element_blank())

defArgs <- function(args, defaultVar)
{
    # fail-safe method that returns the first argument of the vector args,
    # unless it is empty, where it returnts the default
    if(length(args) == 0){
        return(defaultVar)
    }
    return(arg[1])
}

# Quick, convenient saving of a high res png ggplot
f$ggsave2 <- function(f,g,w=7,h=7){ggsave(f,g,width=w,height=h,dpi=600)}

f$readBed <- function(file){
    # Reads a BED formatted file and ouputs a table with the correct column
    # classes and extra information:
    # - read size
    # 
    # Args:
    #   file: filename for the BED file
    # 
    # Returns:
    #   data.frame with columns "seqid, start, end, read, count, strand, size"
    table <-  read.table(file, sep="\t", header=F, colClasses=c("character", 
        "numeric", "numeric", "character", "numeric", "character"), stringsAsFactors=F)
    colnames(table) <- c("seqid", "start", "end", "read", "count", "strand")
    table$size <- nchar(table$read)
    return(table)
}



f$find_overlaps <- function(x,y,colname="ID"){
    overlaps <- y[with(y, seqname == x["seqid"] & (end > x["end"]) & (start < x["end"])),] 
    if(nrow(overlaps) < 1){
        overlaps <- data.frame("ID"=NA)
        colnames(overlaps)[1] <- colname
    }
    overlaps$read <- x[c("read", "start", "end", "seqid")]
    return (overlaps[,c("read", "start", "end", "seqid", colname)])
}


f$read_counts <- function(file, seqs=1){
    # Read a table of counts, putting the reads as the row names.
    #  Args:
    #  file=the filename
    #  seqs=header name/number for the seq column
    #
    read.csv(file, header=TRUE, row.names=seqs);
}

f$read_patman <- function(patfile, seq.columns=c("hit", "query"), counts=FALSE){
    # Read a patman file
    # Args:
    #   patfile=filename
    #   seq.columns=names for the hit and query columns
    #   counts=TRUE if there are counts in the query column. Will add these 
    #          to a seperate column
    tab <- read.table(patfile, sep="\t", header=FALSE, 
        col.names=c(seq.columns, "start", "end", "strand", "mismatches"),
        colClasses=c("character","character","numeric","numeric","character","numeric") )
    if(counts){
        tab$count <- as.numeric(gsub("\\w+\\((\\d+)\\)", "\\1", tab[,2]))
        tab[,2] <- gsub("(\\w+)\\(\\d+\\)", "\\1", tab[,2]);
    }
    tab$size <- nchar(tab[,2]);
    return(tab[tab$size > 0,]);
}

f$getMappingTable <- function(mapped, unmapped){
    # Given two count tables, derive mapping statistics in one dataframe
    m <- getTotals(mapped)
    u <- getTotals(unmapped)
    df <- cbind(mapped=m, unmapped=u, percent=round(100*m/(m+u), 1))
    return(df)
}

f$getTotals <- function(counts){
    data.frame(redundant=apply(counts, 2, sum), unique=apply(counts, 2, function(x) length(x[x > 0])))
}

f$get_sizecounts <- function(readcounts, r=TRUE, read.col=0){
    # Derive a redundant OR non-redundant distribution from
    # any number of count columns.
    # Args:
    #   readcounts=data frame of counts and reads to be used
    #   r=redundant counts if TRUE, non-redundant if FALSE
    #   read.col=column containing reads, 0 if row names
    r.fun <- sum;
    if(!r){
        r.fun <- function(x) length(x[x > 0]); 
    }
    reads <- c()
    counts <- c()
    if(read.col == 0){
        reads <- rownames(readcounts)
        counts <- readcounts
    }
    else{
        reads <- readcounts[,read.col]
        counts <- readcounts[,-read.col]
    }
    counts <- aggregate(counts, by=list(nchar(as.character(reads))),FUN=r.fun)
    colnames(counts)[1] <- "size";
    return(counts);
}

f$dt_get_sizecounts <- function(readcounts, r=TRUE, read.col=0){
    # Derive a redundant OR non-redundant distribution from
    # any number of count columns.
    # Args:
    #   readcounts=data frame of counts and reads to be used
    #   r=redundant counts if TRUE, non-redundant if FALSE
    #   read.col=column containing reads, 0 if row names
    r.fun <- sum;
    if(!r){
        r.fun <- function(x) length(x[x > 0]); 
    }
    sdt <- readcounts[,lapply(.SD, r.fun), .SDcols=(1:ncol(readcounts))[-1], by=nchar(as.character(read))]
    setnames(sdt, 1, "size");
    return(sdt[order(size)]);
}

f$get_sdt <- function(counts, readcol="read")
{
    # returns a melted data.table of r, nr, and complexity
    # for all samples in the count table
    m <- melt_counts(counts, names(counts[,!readcol, with=F]))
    nrfun <- function(x) length(x[x > 0]) 
    sdt <- m[,list(r=sum(count), nr=nrfun(count)), by=list(size=nchar(as.character(read)), sample)]
    sdt[,comp := nr/r]
    return(sdt);
}

f$get_aggregate_counts <- function(readcounts, r=TRUE, read.col=0){
    r.fun <- sum;
    if(!r){
        r.fun <- function(x) length(x[x > 0]); 
    }
    reads <- c()
    counts <- c()
    if(read.col == 0){
        reads <- rownames(readcounts)
        counts <- readcounts
    }
    else{
        reads <- readcounts[,read.col]
        counts <- readcounts[,-read.col]
    }
    counts <- aggregate(counts, by=list(as.character(reads)),FUN=r.fun)
    colnames(counts)[1] <- "read";
    return(counts);
}

f$pool_adjacent_lengths <- function(lengths.tab,namepattern="(.*)"){
    lengths.pooled <- data.frame(length = lengths.tab$length);
    pooled.cols <- pool_adjacent_columns(lengths.tab[,-1],namepattern);
    return (as.data.frame(cbind(lengths.pooled,pooled.cols)));
}

f$pool_adjacent_reads <- function(reads,namepattern="(.*)"){
    pooled.cols <- pool_adjacent_columns(reads,namepattern);
    return(as.data.frame(pooled.cols));
}

f$pool_adjacent_columns <- function(columns, namepattern="(.*)"){
    pooled <- c();
    
    for (c in seq(1,ncol(columns),2)){
        thesepooled <- apply(columns[,c(c,c+1),with=F],1,sum);
        thisname <- colnames(columns)[[c]];
        if(length(namepattern) > 1){
            thisname <- namepattern[[c]]
        }
        else{
            thisname <- sub(namepattern, "\\1", thisname);
        }
        pooled <- cbind(pooled, thesepooled);
        colnames(pooled)[c(ncol(pooled))] <- thisname;
    }
    return(as.data.table(pooled));
}

# MAplot
# First two are for using in the pairs function
f$getMA <- function(x, y=NULL, base=NA, offset=NA, m.offset=offset, a.offset=offset, prop=F, prop.split=NA){
    # Get MA values comparing two count samples
    # args:
    #  x, y vectors of counts
    #  b log base (base now REQUIRED)
    #  offset: offset to add before log to prevent Inf and NaN.
    #           offset is now REQUIRED to force me to know what I am doing
    #  prop: should values be done using proportions? (x/sum(x))?
    if(is.na(offset) & (is.na(m.offset) | is.na(a.offset)))
        stop("offset and m.offset or a.offset not set. Please provide and offset!")
    if(is.na(base))
        stop("Please provide a base")

    M <- getM(x,y,base,m.offset,prop)
    A <- getA(x,y,base,a.offset,prop)
    return(data.table(M=M, A=A))
}

f$getM <- function(x, y=NULL, b, offset=0, prop=F)
{
    if(is.null(y)){
        y <- x[[2]]
        x <- x[[1]]
    }
    x <- x + offset
    y <- y + offset
    if(prop){
        x <- x/sum(x)
        y <- y/sum(y)
    }
    # NOTE: This function is assuming fold change 
    # wants to be used as x -> y i.e. a higher
    # count in y is an upregulation from x -> y
    # and should be positive
    log(y, b) - log(x, b)
}

f$getA <- function(x, y=NULL, b, offset=0, prop=F)
{
    if(offset > 1)
        warning("Offset is > 1 when calculatin A values. Are you SURE this is what you want?")

    if(is.null(y)){
        y <- x[[2]]
        x <- x[[1]]
    }
    x <- x + offset
    y <- y + offset
    if(prop){
        x <- x/sum(x)
        y <- y/sum(y)
    }
    (log(y, b) + log(x, b)) / 2

}

getCiMin <- function(counts, nSD=1, forceSD=F){
    if(length(counts)>2 | forceSD){
        sdc <- sd(counts)
        mc <- mean(counts)
        return(mc-(nSD*sdc))
    }
    else{
        return (min(counts))
    }
}

getCiMax <- function(counts, nSD=1, forceSD=F){
    if(length(counts)>2 | forceSD){
        sdc <- sd(counts)
        mc <- mean(counts)
        return(mc+(nSD*sdc))
    }
    else{
        return (max(counts))
    }
}

getChebychevMax <- function(counts, alpha=0.005, forceSD=F){
    if(length(counts)>2 | forceSD){
        sdc <- sd(counts)
        mc <- mean(counts)
        return(mc+(nSD*sdc))
    }
    else{
        return (max(counts))
    }
}

getChebychev <- function(counts, lower=F, alpha=0.01, forceSD=F)
{
    if(length(counts)>2 | forceSD){
        sdc <- sd(counts)
        mc <- mean(counts)
        alphalim <- sqrt(1/(alpha/2)-1)
        cl <-alphalim * (sdc/sqrt(length(counts)))

        return(ifelse(lower, mc-cl, mc+cl))
    }
    else{
        return(ifelse(lower, min(counts), max(counts)))
    }
}

getChebychevMin <- function(counts, alpha=0.005, forceSD=F){
    if(length(counts)>2 | forceSD){
        sdc <- sd(counts)
        mc <- mean(counts)
        return(mc+(nSD*sdc))
    }
    else{
        return (max(counts))
    }
}

from_zscore <- function(z=1.644852, x){
    return(z*sd(x)+mean(x))
}

f$direction.combo <- function(counts, combos, sample.names, treatment.pat, replicate.pat, offsets, read.cols=c("read"), ...)
{
    default.method <- "none"
    offsets <- offsets
    # offsets: data table of sample, offset, and optional method
    if(!"method" %in% colnames(offsets))
    {
        offsets[,method := default.method]
    }
    offsets[,treatment := gsub(treatment.pat, "\\1", sample)]
    #offsets <- vapply(offsets, function(x) ceiling(mean(x)), numeric(1))
    #out <- data.table(read=character(0), comparison=character(0), d=character(0), ofc=numeric(0), A=numeric(0), ref=numeric(0), obs=numeric(0))
    out <- c()
    # setkey(out,read)
    counts.m <- melt_counts2(counts, sample.names, treatment.pat, replicate.pat)
    setkeyv(counts.m, read.cols)
       bycols <- c("read", "method")

    ci <- counts.m[, c(read.cols, "treatment", "replicate", "count"), with=F]
    if(!"method" %in% colnames(counts)){
        ci[,method := default.method]
    }
    
    readcol.dt <- counts[,read.cols,with=F] 
    setkey(readcol.dt,read)

    read.cols <- unique(c(read.cols, "method"))

    cim <- ci[, list(min=getCiMin(count, ...), max=getCiMax(count, ...)), by=c(read.cols,"treatment")]
    cim[,mid := (min+max)/2]

    combo.names <- c()
    for(i in 1:length(combos))
    {
       combo <- combos[[i]]
       thisO <- offsets[combo %in% treatment][,list(offset=ceiling(mean(offset))), by=method] 
       offset <- thisO[,offset]
       names(offset) <- thisO[,method]

       thiscim <- cim[treatment %in% combo]
       thiscim[, treatment:=factor(treatment, levels=combo)]
       setkeyv(thiscim, bycols)

       # Calculate directions based on CIs
       print("cid")
       cid <- thiscim[treatment==combo[[1]]][thiscim[treatment==combo[[2]]]]
       cid[,d:="S"]
       cid[min > i.max, d:="D"]
       cid[max < i.min, d:="U"]
       cid <- cid[,list(read, method, d)]

       #cid <- thiscim[, list(d=ifelse(min[1]>max[2], "D", ifelse(max[1]<min[2], "U", "S"))), by=bycols]
       setkeyv(cid, bycols)

       plist <- c("U", "U", "D", "D", "S", "S")
       names(plist) <- paste0(rep(combo, times=3), c("max", "min", "min", "max", "mid", "mid"))
       thiscim2 <- melt(thiscim[cid], measure.vars=c("min", "max", "mid"), variable.name="ci", value.name="exp")
       thiscim2[, dir := plist[paste0(treatment, ci)]]
       cip <- dcast.data.table(thiscim2[d==dir], read+method~treatment, value.var="exp")
       setnames(cip, combo, c("ref", "obs"))

       print("cip")
       #cid2 <- cid[d!="S"]
       #cip <- thiscim[cid][,list(ref=ifelse(d[1]=="U", max[1], ifelse(d[1]=="D", min[1], mid[1])), 
       #                          obs=ifelse(d[1]=="U", min[2], ifelse(d[1]=="D", max[2], mid[2])) ), by=bycols]

      # cif <- cim[, list(ofc=getM(mean(c(min[1], max[1])), mean(c(min[2], max[2])), b=2, offset=offsets[i])), by=list(read)]

#      numUniqueReads <- nrow(cip) %% length(offsets[i])
#      
#      if(numUniqueReads !=  0)
#        stop("Number of offsets does not match number of extra variable values")
#      offsetColumn <- rep(offsets[[i]], length.out=nrow(cip))#each=numUniqueReads)

      cip[,offset:=offset[method]]

      print("MA")
      cip[, ofc:=getM(ref,obs, b=2, offset=offset)]
      cip[, A:=getA(ref,obs, b=2, offset=1)]
      # setkey(cif, "read")
      print("out")
      cij <- cid[cip]
       
      #setnames(ci, "V1", paste(i, collapse=".")) 
      combo.names[[i]] <- paste(combo, collapse=".")
      out <- rbind(out, cbind(comparison=combo.names[[i]], cij)) 
    }
    out[,comparison := factor(comparison, levels=combo.names)]
    setkeyv(out, read)
    out <- out[readcol.dt]

    patterns <- out[,list(p=paste0(d, collapse="")), by=bycols][, list(read, p)]

    if(length(out[,unique(method)]) == 1)
    {
        out[,method := NULL]
        cim[,method := NULL]
    }

    setkey(patterns,read)
    patterns <- readcol.dt[patterns]

    #setkey(patterns, read, method)
    #setkey(cim, read, method)
    #samples <- cim[patterns]

    out[,d := factor(d, levels=c("S", "U", "D"))]
    
    
    return(list(fc=out, count=cim, patterns=patterns))
}

direction_cutoff <- function(lofc, cutoff=1){
    d <- list(fc=data.table(lofc$fc), count=data.table(lofc$count), patterns=data.table(lofc$patterns)) # copy data table 
    S.reads <- d$fc[,d=="S"]
    d$fc[,d := factor(ifelse(ofc > cutoff, "U", ifelse(ofc < 0-cutoff, "D", "S")), levels=c("S","U","D"))]
    # set all original S back to S
    d$fc[S.reads, d:="S"]
    bycols <- setdiff(colnames(d$fc), c("comparison", "d", "ref", "obs", "offset", "ofc", "A"))
    d$patterns <- d$fc[,list(p=paste0(d, collapse="")), by=bycols]
    return(d)
}

f$fcpatterns  <- function(counts, sample.names, treatment.pat, replicate.pat, offsets)
{
    counts.m <- melt_counts2(counts, sample.names, treatment.pat, replicate.pat)
}


f$combine_replicates <- function(reads.tab,replicates,columnames){
    sumlibs <- list();
    for( i in 1:length(replicates)){
        sumlib <- apply(reads.tab[, replicates[[i]] ],1,sum);
        sumlibs <- cbind(sumlibs,sumlib);
    }
    colnames(sumlibs) <- columnames;
    return(sumlibs);
}

f$stat_replicates <- function(counts, read.col=1){
# Turn replicate data into mean and standard deviation
df <- data.frame(read=counts$read, mean=apply(counts[,-read.col], 1, mean),
    sd=apply(counts[,-read.col], 1, sd ))
return(df)
}

f$get_stat_replicates <- function(counts, read.col=1){
    names(counts)[read.col] <- "read"
    statdf <- f$stat_replicates(counts, read.col)
    return(statdf[,-1])
}

f$write_counts <- function(counts, file){
    write.csv(counts, file=file, quote=FALSE, row.names=TRUE);
}

f$offset_foldchange <- function(read.counts, control.col, offset=20){
    # compute diff expression using ofc statistic
    control.counts <- unlist(read.counts[,control.col]); 
    ofc <- apply(read.counts, 2, ,
                 c=as.numeric(control.counts), o=offset)

    ofc <- as.data.frame(ofc);
    rownames(ofc) <- rownames(read.counts);
    return(ofc);
}

f$offset_foldchange_combo <- function(read.counts, combinations, colnames, offset=20, log=FALSE, base=2){
    out.frame <- data.frame(matrix(nrow=nrow(read.counts), ncol=length(combinations)) ,row.names=rownames(read.counts));
    for(i in 1:length(combinations)){
        com <- combinations[[i]];
        xc <- as.numeric(read.counts[,com[[1]] ]);
        yc <- as.numeric(read.counts[,com[[2]] ]);
        if(log){
            out.frame[,i] <- log_ofc(xc,yc,base=base,o=offset);
        }
        else{
            out.frame[,i] <- ofc(xc,yc,o=offset);
        }
    }
    colnames(out.frame) <- colnames;
    return(out.frame);
}

f$dt_offset_foldchange_combo <- function(read.counts, combinations = list(), col.names=c(), read.cols="read", offset=20, log=FALSE, base=2){
    # generate all possible combinations
    out.frame <- data.table(matrix(nrow=nrow(read.counts), ncol=length(combinations)), read.counts[,read.cols,with=F]);
#    out.frame[,read.cols:=read.counts[,read.cols,with=F]]
    colv <- col.names
    for(i in 1:length(combinations)){
        com <- combinations[[i]];
        if (length(col.names) == 0){
            colv <- c(colv, paste(colnames(read.counts)[com], collapse="."))
        }
        xc <- as.numeric(read.counts[[com[[1]]]])
        yc <- as.numeric(read.counts[[com[[2]]]])
        if(log){
            out.frame[,i] <- log_ofc(xc,yc,base=base,o=offset);
        }
        else{
            out.frame[,i] <- ofc(xc,yc,o=offset);
        }
    }
    setnames(out.frame, c(colv,read.cols));
    return(out.frame[,c(read.cols,colv),with=F]);
}

f$repvote_pairwise_lofc <- function(countsx, countsy, sig=1, get.consensus=T, ...){

    # columns are different combinations of replicates to find fold change
    # rows are different sequences
    lofc <- matrix(nrow=nrow(countsx), ncol=(ncol(countsx) * ncol(countsy)))

    mcol <- 1
    for(i in 1:ncol(countsx)){
        for(j in 1:ncol(countsy)){
            lofc.list <- log_ofc(countsx[[i]], countsy[[j]], ...)
            lofc[,mcol] <- lofc.list
            mcol <- mcol+1
        }
    }

    if(get.consensus){
        # substitute lofc values for votes.
        # > 1 = 1
        # < -1 = -1
        # anything else = 0
        votes <- ifelse(lofc > 1, 1, ifelse(lofc < -1, -1, 0))
    
        votetable <- apply(votes, 1, function(x) table(factor(x, levels=c(0,1,-1))) )

        consensus <- apply(votetable, 2, function(x){ 
                c1 <- names(x)[x > sum(x)/2]
                if(length(c1) > 1 || length(c1) == 0) 
                    c1 <- 0 
                return(c1)
            })
        return (consensus)
    }
    else{
        return (lofc)
    }
}

f$repvote_lofc <- function(counts, sample.names=NA, combinations=list(),col.names=NA, reps="([\\w+\\.]+)\\.(\\d)", read.cols="read", ...)
{
    outdt <- data.table(counts[,read.cols,with=F])
    # reps: requires a pattern specifying samplename then repicate number
    if(is.na(sample.names[[1]])) stop("Pattern not yet implemented")
    nsamples <- length(sample.names)
    
    for(ci in 1:length(combinations)){
                combo <- combinations[[ci]]
                si <- combo[[1]]
                sj <- combo[[2]]
                icols <- grep(si, colnames(counts), value=T)
                jcols <- grep(sj, colnames(counts), value=T)

                pairwise <- repvote_pairwise_lofc(counts[,c(icols), with=F], counts[,c(jcols), with=F], ...)
                if(is.na(col.names))
                    dtcolname <- paste(si,sj,sep="_")
                else
                    dtcolname <- col.names[[ci]]
                outdt[,c(dtcolname):=pairwise]
    }
    return(outdt)
}

f$ofc <- function(x,c,o){ 
    return( (pmax(as.numeric(x),c)+o) / (pmin(as.numeric(x),c)+o) )
}

f$log_ofc <- function(x,y, o=20, base){
    # log ofc forces a decision for the base
    return( log( (y+o) / (as.numeric(x)+o) ,base) )
}

f$log_offset_foldchange <- function(read.counts, control.col, offset=20){
    control.counts <- unlist(read.counts[,control.col]); 
    ofc <- apply(read.counts, 2, function(x,c,o){ 
            log( (c+o) / (as.numeric(x)+o) ,2)
    },
    c=as.numeric(control.counts), o=offset)
    ofc <- as.data.frame(ofc);
    rownames(ofc) <- rownames(read.counts);
    return(ofc)
}

f$average_repeats <- function(counts, pattern){
    f$summarise_repeats(counts, pattern)
}
f$summarise_repeats <- function(counts, repeats=list(1:3,4:6,7:9,10:12), pattern, fun=mean){
    # Take the mean of replicates, specified as a list and a pattern for each name
    # Args:
    #   counts=counts table
    #   repeats=a list of replicate groups
    #   pattern=the naming pattern to extract for each replicate group
    #   fun=the function to combine replicates with. The mean by default
    if (length(pattern)==1){
        names <- unique(gsub(pattern, "\\1", colnames(counts))); 
    }
    else{
        names <- pattern;
    }
    meancounts <- data.frame(matrix(ncol=length(names),nrow=nrow(counts)));
    colnames(meancounts) <- names;
    rownames(meancounts) <- rownames(counts);

    for(r in 1:length(repeats)){
        i <- repeats[[r]];
        this <- apply(counts[,i],1,fun)
        meancounts[,r] <- this;
    }
    return(meancounts);
}
f$dt_summarise_repeats <- function(counts, groups=NULL, repeats=NULL, pattern=NULL, fun=mean, readcol="read"){
    # Args:
    #   groups: names that can be extracted from the sample names to derive replicates
    #   replicates: sample names grouped by replicates in a list. Used if no groups
    #   pattern: complex regular expression used to extract the grouped sample names

    if (!is.null(pattern)){
        names <- unique(gsub(pattern, "\\1", colnames(counts[,!readcol,with=F]))); 
    }
    else{
        names <- groups
    }
    reads <- counts[,readcol,with=F]

    if(!is.null(groups))
    {
       repeats <- lapply(groups, function(x) grep(x, names(counts[,!readcol,with=F]))) 
    }
    else{
        if(is.null(pattern)){
            stop("Both groups and pattern are undefined. Please define one of them")
        }
    }

    meancounts <- data.table(matrix(ncol=length(names),nrow=nrow(counts)));
    setnames(meancounts, names);
    thesecounts <- as.matrix(counts[,!readcol,with=F])
    for(i in 1:length(repeats)){
        r <- c(repeats[[i]]);
        if(as.character(quote(fun)) == "mean")
            this <- rowMeans(thesecounts[,r,drop=F]) 
        else
            this <- apply(thesecounts[,r,drop=F],1,fun)

        meancounts[,i := this, with=F];
    }
    return(data.table(reads, meancounts))
}

# functions for quality analysis using ggplot2
#  READ IN DATA
f$read_data <- function(data){
   return( read.table(data, header=TRUE, row.names=1) ); 
}

f$read_lengths <- function(file){
   return( read.table(file, header=TRUE) );
}

f$plot_sample_counts <- function(samples, reads, gridx=3, gridy=3, names=c()){
    par(mfrow=c(gridy,gridx));
    for(i in seq(1,ncol(reads),samples)){
        for (x in i:(i+(samples-2))){
            for (y in (x+1):(i+(samples-1))){
                #print(c(i,x,y));
                xlabel <- "";
                ylabel <- "";
                if(length(names) > 0){
                    xlabel <- names[[x]];
                    ylabel <- names[[y]];
                }
                else{
                    xlabel <- paste("log s",x,sep="");
                    ylabel <- paste("log s",y,sep="");
                }
                plot(reads[,x], reads[,y], log ="xy", pch=".", xlab=xlabel, ylab=ylabel); 
            }
        }
    }
}

f$plot_sample_matrix <- function(samples, reads, plotlog=TRUE,groupnames=c()){
    groups <- seq(1,ncol(reads),samples); 
    for(x in 1:length(groups)){
        i <- groups[[x]];
        print(c(i,i+(samples-1)));
        readsunique <- unique(reads[,i:(i+(samples-1))]);
        print("Plotting...");
        if(plotlog){ 
            readsunique <- log10(readsunique);
        }
        groupname <- c();
        if(length(groupnames) == 0){
            groupname <- paste("Group", x, sep=" ")
        }
        else{
            groupname <- groupnames[[x]];
        }
        plot(readsunique, pch=".", xaxt="n", yaxt="n", main = groupname);
    }
}

f$plot_length_complexity <- function(rlengths, nrlengths, replicates=4, groupnames=c(), maintitle="Complexity Index Across Size Classes",lengthcol="length"){
    print("Hello");
    clengths <- nrlengths[,-1]/rlengths[,-1];
    clengths$length <- rlengths[,lengthcol];
    mclengths <- melt_with_groups(clengths, replicates, groupnames)

    if(replicates > 1){
        g <- ggplot(mclengths, aes(factor(length), value, group=replicate, colour=replicate)) + geom_line() + scale_x_discrete("Read Length") + scale_y_continuous("Complexity Ratio") 
            + facet_wrap(~ groups, ncol=1, scales="free_x") + scale_colour_discrete("Replicate") + labs(title=maintitle) 
    }
    else{
        g <- ggplot(mclengths, aes(factor(length), value, group=groups, colour=groups)) + geom_line() + scale_x_discrete("Read Length") + scale_y_continuous("Complexity Ratio") 
    }
        print(g);
}        

# Plot stacked length distributions across samples
f$stacked_lengths <- function(lengths,graphtitle){
    mlengths <- melt(lengths,id.vars = "length");
    print( ggplot(subset(mlengths, length %in% 20:24), aes(variable, value, fill=factor(length))) + geom_bar(stat="identity",position="fill") + 
        scale_x_discrete("Samples") + scale_y_continuous("Proportion of Reads") + scale_fill_discrete("Length") + labs(title=graphtitle) );
}


# "Quick glance" page of all length distributions
f$length_dists <- function(rlengths, nrlengths, intervals=2, graphtitle=""){
    cols <- c("length", 1:(ncol(rlengths)-1));
    colnames(rlengths) <- cols;
    colnames(nrlengths) <- cols;
    combined.lengths <- rbind(melt(rlengths,id.vars="length"), melt(nrlengths,id.vars="length")) 
    combined.lengths$count <- rep(c("redundant","unique"),each=(nrow(combined.lengths)/2));
    breaks <- seq(rlengths[1,1], rlengths[nrow(rlengths),1],intervals);
    print( ggplot(combined.lengths, aes(factor(length),value, fill=count))+geom_bar(stat="identity",position="dodge")+facet_wrap(~variable) 
    + scale_x_discrete("Read Lengths", breaks=breaks) + scale_y_continuous("Frequency") + labs(title=graphtitle)   );
}


# Line plots of Length distributions compared over samples
f$length_plot_over_samples <- function(lengths, sizecol = "size"){
    lengths.percent <- length_percentages(lengths);
    colnames(lengths.percent) <- c(sizecol, 1:(ncol(lengths)-1));
    mlp <- melt(lengths.percent,id.vars=sizecol);
    #colnames(mlp)[1] <- "size"
    print( ggplot(mlp,aes(factor(variable),value,group=factor(size),colour=factor(size))) + geom_line() + geom_text(data=mlp[mlp$variable==max(as.numeric(mlp$variable)),],aes(label=size)) + 
        scale_x_discrete("Samples") + scale_y_continuous("Proportion of Read Counts") + labs(title="Changes in Size Class Composition Over All Samples") + scale_colour_discrete("Size Class") );
}

# Size class distribution using line plots with all samples
f$length_dist_line_plot <- function(lengths, repeats=3, groupnames=c(), maintitle="Size Class Distribution Over All Samples"){
    mlp <- melt_with_groups(lengths, repeats, groupnames);
    gg <- gg_length_dist_line_plot(mlp, maintitle);
}

f$gg_length_dist_line_plot <- function(mlp, maintitle="")
{
    return( ggplot(mlp,aes(factor(size),value, group=factor(variable), colour=factor(groups))) + geom_line() +
        scale_x_discrete("Read length (nt)") + scale_y_continuous("Read count", labels=scales::comma) + labs(title=maintitle)) ;
}

f$length_dist_line_plot_faceted <- function(lengths, repeats=3, groupnames=c(), maintitle="Size Class Distributions"){
    mlp <- melt_with_groups(lengths, repeats, groupnames);
    gg <- gg_length_dist_line_plot_faceted(mlp, maintitle) 
}

f$gg_length_dist_line_plot_faceted <- function(mlp, maintitle=""){
    return( ggplot(mlp,aes(factor(size),value,group=factor(variable),colour=factor(replicate))) + 
        geom_line() + facet_wrap(~groups,ncol=1,scales="free_x") +
        scale_x_discrete("Read length (nt)") + scale_y_continuous("Read Counts") + 
        scale_colour_discrete("Replicate") + labs(title=maintitle) );
}

f$melt_counts <- function(counts, countcols){
    # shorthand to melt count columns. variable and value columns
    # are automatically labelled "sample" and "count"
    melt(counts, measure.vars=countcols, variable.name="sample", value.name="count")
}
melt_counts2 <- function(counts, sample.names, treatment.pat, replicate.pat)
{
    #melt_counts with extra processing to extract treatment and replicate names
    counts.m <- melt_counts(counts, sample.names)
    tnames <- gsub(treatment.pat, "\\1", sample.names)
    names(tnames) <- sample.names
    rnames <- paste0("R", gsub(replicate.pat, "\\1", sample.names))
    names(rnames) <- sample.names
    counts.m[, treatment := factor(tnames[sample], levels=unique(tnames))] 
    counts.m[, replicate := factor(rnames[sample], levels=unique(rnames))]
    #counts.m[,sample := NULL]
    return(counts.m)
}

f$reshape_to_reps <- function(counts, sample.names, treatment.pat, replicate.pat, othercols="read")
{
    # reshape a table of counts to replicate columns with variable treatment column based
    # on pattern of current column names
    formula.str <- paste0(paste(c(othercols,"treatment"), collapse=" + "), " ~ replicate")

    counts.m <- melt_counts2(counts, sample.names, treatment.pat, replicate.pat)

    rep.names <- unique(counts.m[,replicate])

    counts.r <- dcast.data.table(counts.m, as.formula(formula.str), value.var="count")
    counts.r <- cleancounts(counts.r, count.cols=rep.names)
#    counts.r[, treatment := factor(treatment, levels=sample.names)]

    return(counts.r)
}

f$melt_with_groups <- function(lengths, repeats=3, groupnames=c()){
    groups <- (ncol(lengths)-1)/repeats;
    if(length(groupnames)==0){
        groupnames = 1:groups;
    }
    groupnames <- factor(groupnames,levels=groupnames) # making sure their order is the same as the order of input
    mlp <- melt(lengths,id.vars="size");
    mlp$groups <- rep( groupnames[1:groups], each=length(unique(mlp$size))*repeats);
    numgroups <- ncol(lengths)/repeats;
    numlengths <- nrow(lengths);
    #mlp$replicate <- paste("R",rep( rep(1:repeats, each=numlengths), times=groups),sep=""); 
    mlp$replicate <- paste("",rep( rep(1:repeats, each=numlengths), times=groups),sep=""); 
    return(mlp);
}

# Functions to retrieve a percentage table from a length distribution table
f$length_percentages <- function(lengthdist){
    return( as.data.frame(cbind(size=lengthdist$size,sapply(lengthdist[,-1],list_prop))) );
}

f$list_prop <- function(alist){
    total <- sum(alist);
    return( sapply(alist, function(x,y){ (x/y)*100 }, y=total) );
}

f$jaccard <- function(counts, sample.pattern="(\\w+_\\w)(\\d)", num.points=100){
    # sample.pattern indicates where the sample and replicate labelling is in the names
    # assumes the sample labelling comes first
    
    ranked.reads <- lapply(counts[,!"read",with=F], rank_reads, reads=counts[,read])
    names(ranked.reads) <- names(counts[,!"read",with=F])

    dt <- list()
    for(i in 2:length(ranked.reads)){
        for(j in 1:i){
            if(i!=j){
                namei <- names(ranked.reads)[i]
                namej <- names(ranked.reads)[j]

                m <- str_match(namei, sample.pattern) 
                samplei <- m[,2]
                repi <- m[,3]

                m <- str_match(namej, sample.pattern) 
                samplej <- m[,2]
                repj <- m[,3]

                ctype <- ifelse(samplei==samplej, samplei, paste(samplei, samplej, sep=".")) 
                
                thisdt <- cbind(pairwise_jaccard(ranked.reads[[i]], ranked.reads[[j]], num.points), comparison=paste(namei,namej,sep="."), 
                    comparison.type=ctype) 
                dt <- rbind(dt, thisdt)
            }
        }
    }
    dt <- as.data.table(dt)
    dt[,comparison.group := ifelse(str_detect(comparison.type, "\\."), "treatment comparison", "replicate comparison")]
    return(dt)
}

f$gg_jaccard_series <- function(jaccard){
    ggplot(jaccard, aes(n, jac)) + geom_line(aes(group=comparison, colour=comparison.type)) + 
        facet_wrap(~ comparison.group, ncol=1) + 
        labs(x="Number of top sequences", y="Jaccard index", colour="Comparison type")
}

f$rank_reads <- function(counts, reads){
    reads[counts!=0][order(-counts[counts!=0])] 
}

f$pairwise_jaccard <- function(rreadsx, rreadsy, num.points=100){
    xo <- rreadsx
    yo <- rreadsy
    #xo <- reads[countx!=0][order(-countx[countx!=0])]
    #yo <- reads[county!=0][order(-county[county!=0])]

    maxN <- max(length(xo),length(yo))
    intervals <- floor(maxN/num.points)

    jac <- data.table()
    #print(maxN)
    #print(intervals)
    for(n in seq(intervals, maxN, intervals)){
        thisx <- xo[1:n]
        thisy <- yo[1:n]
        jac <- rbind(jac, list(n=n, jac=length(intersect(thisx, thisy))/length(union(thisx,thisy))))    
    }
    return(jac)
}

# Heat Map of Jaccard Index
f$jaccard_table <- function(reads, threshold=500){
    dframe <- data.frame(values=names(reads));
    #get list of top 2000 reads per sample
    top.tab <- list();
    for (r in 1:length(reads)){
        top.tab[[ names(reads[r]) ]] <- rownames(head(reads[order(-reads[r]),],n=threshold));
    }

    for (r in 1:length(reads)){
        S1 <- names(reads[r]);
        values <- c();
        for (s in 1:length(reads)){
            S2 <- names(reads[s]);
            u <- union(top.tab[[S1]], top.tab[[S2]]); 
            n <- intersect(top.tab[[S1]], top.tab[[S2]]);
            jaccard <- length(n)/length(u); 
            values <- c(values,jaccard);
        }
        dframe <- cbind(dframe,as.numeric(values));
    }
    colnames(dframe) <- c("values",names(reads));
    return(dframe)
}

f$dt_jaccard_table <- function(reads, threshold=500){
     readnames <- reads[,read]
     samplenames <- colnames(reads)[-1]
     ranked.reads <- data.table(apply(reads[,-1,with=F], 2, function(x) head( readnames[order(-x)], n=threshold )))

    jac <- function(x) {
        r1 <- ranked.reads[[x[[1]]]]
        r2 <- ranked.reads[[x[[2]]]]
        u <- union(r1,r2)
        n <- intersect(r1,r2)
        jaccard <- length(n)/length(u); 
        return(jaccard)
    }
     #melted <- t(  apply(combn(samplenames,2), 2, jac) )
     jacct <- data.table(expand.grid(samplenames,samplenames))
     
     jacct[,value := apply(.SD, 1,function(x) jac(x))]
     return(jacct)
}

f$dt_jaccard_heat_map <- function(dt.jtable, title=NULL, text=T, xrot=F, limits=c(0,1), grouping.lines=c(), diagx=F, diagx.size=10){
    g <- ggplot(dt.jtable[Var1!=Var2], aes(Var1, Var2)) + geom_tile(aes(fill=value), colour="white") + 
        scale_fill_gradient(low="white", high="black", limits=limits, guide=F) + 
        labs(x="Libraries", y="Libraries", fill="Jaccard\nindex", title=title) + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
        theme(panel.border=element_blank(), panel.grid=element_blank()) 

    if(text){
        g <- g + geom_text(aes(label=sprintf("%.0f", value*100)), colour="white", size=2.5) 
    }
    if(xrot){
        g <- g + theme(axis.text.x=element_text(angle=60, hjust=1));
    }
    if(diagx){
        g <- g + geom_text(data=dt.jtable[Var1==Var2], label="X", colour="grey", size=diagx.size)
    }
    if(length(grouping.lines)>0){
        grouping <- grouping.lines+0.5
        g <- g + geom_hline(yintercept=grouping) + geom_vline(xintercept=grouping);
    }
    return(g + coord_fixed())
}

f$square_chart_area <- function(gg_with_legend, np=T){
    # plots a ggplot with a legend so that the chart area is still in proportion
    # specifically for jaccard heat map for now
    ggjac <- gg_with_legend
    tmp <- ggplot_gtable(ggplot_build(ggjac))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]

    if(np)
        grid.newpage()
    col1 <- 0.85
    pushViewport( viewport( layout = grid.layout(1, 2, 
                            widths = unit(c( col1, 0.1 ), "npc" ), 
                            heights = unit( c(col1, 1), "npc"), 
                            respect = matrix(rep(1,2),1 ) ) ) )
    print(ggjac + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1 ) )
    upViewport(0)
    vp3 <- viewport( width = unit(0.15, "npc"), x=0.9, y=0.5)
    pushViewport(vp3)
    grid.draw( legend )
    popViewport()
}


f$jaccard_heat_map <- function(jtable, title, lowcol="white", hicol="black", text=FALSE, xrot=FALSE, limits=c(0,0.99), grouping=c()){
    jac.tab.melt <- melt(jtable,id.vars="values");
    jac.tab.melt$values <- factor(jac.tab.melt$values,levels=c(as.character(jtable$values)));
    jplot <- ggplot(jac.tab.melt, aes(x=variable,y=values)) + geom_tile(aes(fill=value), colour = "white") +
        scale_fill_gradient(low = lowcol, high = hicol, limits=limits) + 
        scale_x_discrete("Samples") + scale_y_discrete("Samples") + labs(title=title, fill="Jaccard index");
    if(text){
        jplot <- jplot + geom_text(aes(label=sprintf("%.1f",value)), colour="white", size=2.5); 
    }
    if(xrot){
        jplot <- jplot + theme(axis.text.x=element_text(angle=60, hjust=1));
    }
    if(length(grouping)>0){
        grouping <- grouping+0.5
        jplot <- jplot + geom_hline(aes_now(yintercept=grouping)) + geom_vline(aes_now(xintercept=grouping));
    }
    return(jplot)
}

f$aes_now <- function(...) {
  structure(list(...),  class = "uneval")
}

f$getRegion <- function(region,seqid=region$seqid,start=region$start,end=region$end, pat, counts){
    patreg <- pat[pat$hit == seqid & pat$start >= start & pat$end <= end,]
    return(merge(patreg, counts, by.x="query", by.y="seq"))
}

spBind <- function(..., spList=NULL, libname="lib"){
    # bind rows of multiple data.frames, using named arguments
    if(is.null(spList)){
        spList <- list(...)
    }
    names.sp <- names(spList)
    for (i in 1:length(spList)){
        spList[[i]][[libname]] <- names.sp[i]
    }
    return(do.call("rbind", spList))
}

####
# GGPLOT MULTIPLOT FUNCTION
####

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, by.row=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
f$multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
    }
#    print(layout)

    if (numPlots==1) {
        print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col))
        }
    }
}

ma.rep.counts <- function(counts, sample.names, treatment.names, sample.reg, rep.reg){
    counts.m <- melt_counts(counts, sample.names)
    counts.m[, treatment := gsub(sample.reg, "\\1", sample)]
    counts.m[, replicate := paste0("R",gsub(rep.reg, "\\1", sample))]
    counts.r <- dcast.data.table(counts.m, read + treatment ~ replicate, value.var="count")
    counts.r[, treatment := factor(treatment, levels=treatment.names)]
    ma.list <- list()
    reps <- names(counts.r)[-c(1,2)]
    nzcounts.r <- nzcounts(counts.r, reps)
    for(i in 1:length(reps)){
        for(j in i:length(reps)){
            if(i!=j){
                rcmp <- reps[c(i,j)]
                ma.list[[paste(rcmp, collapse=",")]] <- data.table(read=nzcounts.r[,read], treatment=nzcounts.r[,treatment],
                    ref=nzcounts.r[[reps[[i]]]], obs=nzcounts.r[[reps[[j]]]], getMA(nzcounts.r[,rcmp,with=F], base=2,offset=0))
            }
        }
    }
    return(spBind(spList=ma.list, libnames="comparison"))
}




f$docsv <- function(x,file){write.csv(x, file=file, row.names=F, quote=F)}
f$dopng <- function(filename, width=6, height=6){png(filename, units="in", res=300)}

f$total.dt <- function(counts, r=T){
    t.fun <- sum
    if(!r){
        t.fun <- function(x) {sum(x>0)}
    }
    total.reads <- apply(counts[,-1,with=F], 2, t.fun)
    return(total.reads)
}

f$mapping.summary <- function(all.counts, genome.counts, r=T){
    all.total <- total.dt(all.counts, r)
    genome.total <- total.dt(genome.counts, r)
    summary <- data.table(samples=colnames(all.counts)[-1], total=all.total, mapped=genome.total, percentage=100*genome.total/all.total)
    return(summary)
}

# sizeclass distributions
getScd <- function(unmapped, mapped, sample.names, t.pat, r.pat){
    counts.scd <- melt_counts2(mapped, sample.names, t.pat, r.pat)[,list(nr=sum(count > 0), 
        r=sum(count)), by=list(size=nchar(read), sample, treatment, replicate)]
    counts.scd[, comp := nr/r]

    counts.um.scd <- melt_counts2(unmapped, sample.names,t.pat,r.pat)[,list(nr=sum(count > 0),
        r=sum(count)), by=list(size=nchar(read),sample, treatment, replicate)]
    counts.um.scd[, comp := nr/r]

    counts.scd.all <- spBind(m=counts.scd, um=counts.um.scd, libname="mapped")

    counts.prop.mapped.scd <- counts.scd.all[,list(r=r[1]/r[2], nr=nr[1]/nr[2]),by=list(size,sample,treatment,replicate)]
    return(list(counts=counts.scd.all, prop=counts.prop.mapped.scd))
}
gg_scd <- function(scd, y, ylab, guide=F){
    g <- ggplot(scd, aes(x=factor(size), group=sample, colour=treatment)) + geom_line(size=0.5, aes_string(y=y)) +
        theme_nar + theme(panel.grid=element_line(colour="grey")) + labs(x="Read length", y=ylab)

    if(!guide)
        g <- g + scale_colour_discrete(guide=F)
    return(g)
}

# if too stringent, attempting group using proximity (200nt)
# calculate a lag difference - n[,end] - n-1[,start]
de_loci <- function(aligned.reads, proximity=0, othercols = c()){
    aligned.reads <- aligned.reads[order(seqid,start,end)]
    region.gaps <- data.table(seqidA=aligned.reads[-1,seqid], seqidB=aligned.reads[-nrow(aligned.reads),seqid], 
        idA=2:nrow(aligned.reads), idB=1:(nrow(aligned.reads)-1),
        gap=aligned.reads[-1,start] -  aligned.reads[-nrow(aligned.reads),end])

    region.gaps[,grouped := seqidA==seqidB & gap < proximity]
    print(region.gaps,20)

    # run-length encoding of gaps between sequences being < 200
    # T if gap is < 200
    # F if gap is > 200
    region.rle <- rle(region.gaps[,grouped])
    regions.grouped <- data.table(seqid=character(0), start=numeric(0), end=numeric(0), strand=character(0), type=character(0),  name=character(0), 
                                  arm=character(0),
                                    N=numeric(0), regionid=numeric(0))

    thisid <- 1
    gid <- 1
    # iterate over runs of T or F
    # T means that all ids involved are linked, including with the previous id
    # F means that all ids involved are seperate regions each, excluding the very last id
    for(i in 1:length(region.rle$values)){
        # thislen is the size of the run of T or F
        thislen <- region.rle$length[i]

        # if this is a run of T (seqs with gaps below threshold and should be grouped)
        if(region.rle$values[i]==T){

            # special case:
            # if(thisid+thislen > nrow(region.gaps))
            #    thislen <- thislen-1

            thisgroup <- region.gaps[thisid:(thisid+(thislen-1)), list(first=min(idB), last=max(idA))]
            #aligned.reads[thisgroup[,first]:thisgroup[,last], regionid := gid]

            thisregions <- aligned.reads[thisgroup[,first]:thisgroup[,last], 
                list(seqid=unique(seqid),start=min(start), end=max(end), strand=mymode(strand), type=mymode(type), name=mymode(name), arm=mymode(arm), N=.N) ]
            thisregions[is.null(arm),arm:=NA]
                print(thisregions)

            thisregions[,regionid := gid:(gid+nrow(thisregions)-1)]
            regions.grouped <- rbind(regions.grouped, thisregions)
            gid <- max(thisregions[,regionid])+1;
            print(c("T", thisregions[1,name], thisid, thisid+thislen-1))
            thisid <- thisid + thislen
        }
        else{
            # this is a run of single seqs
            if(thislen > 1){
                # for each id except the last id here, make it a seperate region
                Fids <- thisid:(thisid+(thislen - 2))
                for(j in Fids ){
                    regions.grouped <- rbind(regions.grouped, aligned.reads[region.gaps[j, idA],
                                            list(seqid, start, end, strand, type=type, name, arm=arm, N=1, regionid=gid)])
                    print(c("F", aligned.reads[j, name], j))
                    gid <- gid + 1 
                }
            }
            thisid <- thisid + (thislen)
        }
    }
    return(regions.grouped)
}

get_region_counts <- function(regions, bed, counts){
    setkey(counts,read)

    region.counts <- c() 
    for(i in 1:nrow(regions)){
        thisregion <- regions[i,]
        region.reads <- prof$getBedRegion(bed, thisregion)
        setkey(region.reads, read)

        read.counts <- counts[region.reads,nomatch=0]
        total.counts <- read.counts[,lapply(.SD, sum),.SDcols=colnames(counts[,!"read",with=F])]
        
        regionWithCounts <- data.table(thisregion, total.counts)
        region.counts <- rbind(region.counts, regionWithCounts)
    }
    return(region.counts)

}

get_all_region_reads <- function(region, bed, counts){
    setkey(counts,read)
    regions <- c() 
    
    for(i in 1:nrow(region)){
        thisregion <- region[i,]
        region.reads <- prof$getBedRegion(bed, thisregion)
        setkey(region.reads,read)
        read.counts <- counts[region.reads,nomatch=0]
        regions <- rbind(regions, cbind(thisregion, read.counts))
    }
    return(regions)
}

rep_ma <- function(counts.mapped, sample.names, t.pat, r.pat){
    counts.mapped.m <- melt_counts2(counts.mapped, sample.names, t.pat, r.pat)
    rep.names <- unique(counts.mapped.m[,replicate])
    treatment.names <- unique(counts.mapped.m[,treatment])
    counts.mapped.r <- dcast.data.table(counts.mapped.m, read + treatment ~ replicate, value.var="count")
    counts.mapped.r <- cleancounts(counts.mapped.r, count.cols=rep.names)
    counts.mapped.r[, treatment := factor(treatment, levels=treatment.names)]
    rep.pairs <- combn(as.character(rep.names),2)
    return(ma_from_reps(counts.mapped.r, c("read", "treatment"), rep.pairs))
}

ma_from_reps <- function(counts.rep, read.cols, rep.pairs){
    ma.list <- list()
    for(r in 1:ncol(rep.pairs)){
       thispair <- paste(rep.pairs[,r],collapse=",") 
       ma.list[[thispair]] <- data.table(counts.rep[,c(read.cols, rep.pairs[,r]),with=F], 
                              getMA(counts.rep[,rep.pairs[,r],with=F], o=1, b=2))
       setnames(ma.list[[thispair]], rep.pairs[,r], c("ref", "obs"))
    }
    rep.mas <- spBind(spList=ma.list, libname="rep.pair")
    rep.mas <- rep.mas[ref!=0|obs!=0]
    rm(ma.list)
    return(rep.mas)
}

# detach functions and then reattach to update
while("f" %in% search())
    detach("f")
attach(f)

mymode <- function(x){
    names(sort(-table(x)))[1]
}
