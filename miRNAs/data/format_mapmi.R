library(data.table)
# R was used to format and filter the table output by mapmi
mapmi.output.file <- "mapmi_output.csv"
mapmi <- as.data.table(read.csv(mapmi.output.file), header=F)
setnames(mapmi, c("name", "query", "species", "extension", "mismatches", "seqid", "strand", 
                  "mature.start", "mature.end", "mature.length",
    "pre.start", "pre.end", "pre.length", "structure.mismatches", "mfe", "score", "small.extension", 
    "long.extension", "mismatch.penalty", "max.mismatches", "pre.seq", "dust.fraction"))
mapmi2 <- mapmi[,list(name, query, mismatches, seqid, strand, mature.start, mature.end, mature.length, 
                      pre.start, pre.end, pre.length, mfe, score, pre.seq, dust.fraction)]
mapmi2[,id := gsub("([^\\s])\\s.+", "\\1", name)]
mapmi2[, sp.id := gsub("(\\w+)-.+", "\\1", id)]
mapmi2[, mir.id := gsub("\\w+-(.+)", "\\1", id)]

# Allow slight mature variant indicators (letters after mir number)
mapmi2[,mir.name := str_extract(mir.id,perl("((miR|let)(-?\\d+\\w?|-iab)|bantam)"))]
mapmi2[, sp.name := gsub("\\S+\\s\\S+\\s(\\S+\\s\\S+)\\s\\S+", "\\1", name)]
setkey(mapmi2, query)
mapmi2[id=="asu-miR-79-5p", mir.name:="miR-9"]

# extract only hexapoda mapmi
hexapoda.sp <- c("ame", "aae","aga","api","bmo","cqu","dan","der","dgr","dme","dmo","dpe","dps","dse",
                 "dsi","dvi","dwi","dya","hme","lmi","mse","ngi","nlo","nvi","pxy","tca")
mapmi2.hx <- mapmi2[sp.id %in% hexapoda.sp]

docsv(mapmi2.hx, "mapmi_formatted.csv")
