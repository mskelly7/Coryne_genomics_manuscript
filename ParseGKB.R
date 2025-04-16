
# Parse antismash outputs

# https://thackl.github.io/gggenomes/index.html
require(gggenomes)

require(read.gb)

FindFeature = function(start, stop, type, gblist) {
    # in a file imported using read.gb, find the feature(s) matching the start or end and type provided in the function call arguments
    recordid <- c()
    # position string to match, include + and - strand orientations.
    # assume there won't be a by change complement match
    posStringA <- paste(c(start), "..", c(stop), sep = "")
    posStringB <- paste("complement(", posStringA, ")", sep = "")
    for ( i in 1:length(gblist[[1]][["FEATURES"]]) ) {
        if( (any(gblist[[1]][["FEATURES"]][[i]][["Qualifier"]] == posStringA) | any(gblist[[1]][["FEATURES"]][[i]][["Qualifier"]] == posStringB) ) & any(gblist[[1]][["FEATURES"]][[i]][["Location"]] == type)) { 
            recordid <- c(recordid, i)
        } 
    }
    # if a record was found, get the protein sequence. protein sequence should be the first translation in the series.  
    # all subsequent translations correspond to other feature types not correctly imported by the read.gb function
    proteinSeq <- ifelse(length(recordid)==1,
        gblist[[1]][["FEATURES"]][[recordid]][["Qualifier"]][which(gblist[[1]][["FEATURES"]][[recordid]][["Location"]] == "translation")[1]],
        "not found"
    )
    return(proteinSeq)
}

ExctractProtoCore = function( gblistA ) {
    pc <- gblistA[which(gblistA$type == "proto_core"),, drop = FALSE]
    return(pc)
}

ExtractFeatureCDS = function(pc, gblistA, gblistB) {
    pc <- data.frame(pc)
    # check exact match for start or end
    pcCDS <- gblistA[which( (gblistA$start == pc$start | gblistA$end == pc$end) & gblistA$type == "CDS"),, drop = FALSE]
    #pcCDS = pcCDS[, which(colSums(is.na(pcCDS)) < nrow(pcCDS))]

    # it is possible to  have more than one CDS feature associated with a single pc
    if(nrow(pcCDS) > 0) { 
        pcCDS$translation <- sapply(1:nrow(pcCDS), function(i) FindFeature(pcCDS$start[i], pcCDS$end[i],pcCDS$type[i], gblistB) )
    }
    return(pcCDS)
}

ParseSampleXRegion = function(fh) {
    # subject and region from fh
    subjectid <- gsub( "\\/.*$", "", gsub("^.*prokka\\.", "", fh)) 
    regionName <- gsub(".gbk", "", basename(fh))

    gblistA <- read_gbk(fh)
    gblistB <- read.gb(fh)

    # get the proto core feature, might be more than 1 per region
    proto_core <- ExctractProtoCore(gblistA )
    # for each proto core feature, grab all CDS features and the corresponding protein sequence
    proto_core_CDS <- do.call(rbind, lapply(1:nrow(proto_core), function(i) {
        out <- ExtractFeatureCDS(proto_core[i,], gblistA ,gblistB)
        if(nrow(out) > 0) { 
            out <- merge(proto_core[i,], out, by = "seq_id", suffix = c("_protoCore", "_CDS"), all.x = TRUE)
        } else {
            out <- proto_core[i,]
        }
        if(!"sec_met_domain" %in% colnames(out)) {
            out$sec_met_domain <- NA
        } 
        if(! "ec_number" %in% colnames(out)) {
            out$ec_number <- NA
        }
        out
    }))

    proto_core_CDS$subjectid <- subjectid
    proto_core_CDS$regionName <- regionName

    return( proto_core_CDS)
}

# potential future work flow - 
    # get the protein translations for all CDS features listed in region 
    # all_CDS = gblistA[which(gblistA$type == "CDS"),]
    # all_CDS$translation = sapply(1:nrow(all_CDS), 
    #   function(i) FindFeature(all_CDS$start[i], all_CDS$end[i], all_CDS$type[i], gblistB) )



# test as single file 
# fh = "~/Documents/2022_Coryne_MKelly/Data/WGS/processed/antismashOutputs/prokka.05-003/GEGLKPGJ_30.region001.gbk"
# test <- ParseSampleXRegion(fh)

# test a whole sample
# files <- list.files("~/Documents/2022_Coryne_MKelly/Data/WGS/processed/antismashOutputs/prokka.05-002/", 
#        pattern = "region.*gbk",
#        full.names = TRUE)
#test <- lapply(files, ParseSampleXRegion)
#names(test) <-gsub(".gkb", "", basename(files))

# directory to the directory containing the antismash output folders / files. Assumes one smaple per antismash folder. 
antismashDir = "Data/WGS/processed/antismashOutputs"

# get full paths for all antismash region.gbk files
files <- list.files(antismashDir , 
        pattern = "region.*gbk",
        recursive = TRUE, 
        full.names = TRUE)
# parse the proto_core and associated protien sequence(s) for reach region file.
# some regions have mutiple proto_cores
# some proto_cores have multikple protiens sequences 
# records are expanded to the unique protein sequence for the proto_core 'hit'
test <- lapply(files, ParseSampleXRegion)
names(test) <-gsub(".gkb", "", basename(files))

# master set of columns
cols = unique(unlist(lapply(test, names)))
# add any missing columns
for (i in 1:length(test)) {
    test[[i]][, cols[which(!cols %in% colnames(test[[i]]))]] = NA 
}
# convert from list to data frame format
MasterOut = do.call(rbind, test)

#remove columns NA in all rows
MasterOut[, which(colSums(is.na(MasterOut)) == nrow(MasterOut))] = NULL
MasterOut_v2 <- apply(MasterOut,2,as.character)
# write.table(data.frame(MasterOut_v2), "All_antismashProto_cores.txt", row.names = FALSE, sep = "\t")

# verify that the data can be imported and match 
testIn = read.table("All_antismashProto_cores.txt", sep = "\t", header= TRUE)
table(sapply(colnames(testIn), function(i) all(testIn[,i] == MasterOut_v2[,i], na.rm = TRUE) ))