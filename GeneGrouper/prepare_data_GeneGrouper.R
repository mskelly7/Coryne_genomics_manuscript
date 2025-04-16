library(geneviewer)
library(GenomicFeatures)
library(Biostrings)

# Prepare data for GeneGrouper

# Copy gbk files
rm(list= ls(all= TRUE))

bgc_type="terpene"

dat <- data.table::fread("UniqueSequences.csv")

head(dat)

table(dat$product_protoCore)

tmp <- dat[which(dat$product_protoCore == bgc_type & dat$contig_edge == F), ]

head(tmp)


gbk_dir <- paste0("gene_grouper/", bgc_type, "/genomes_gbk/")

if(!dir.exists(gbk_dir)) {
  dir.create(gbk_dir)
}

for(i in 1:nrow(tmp)) {
  # Copy genome gbk file from the shared drive:
  com <- paste("cp", tmp$gbk_filename[i], paste0(gbk_dir, "/", tmp$UniqueClusterID[i], ".gbk"))
  system(com)
}


# Convert gbk to gbff (use python script for this):
# python convert_to_gbff.py

# Extract query genes
table(tmp$name_CDS)
query_genes <- names(table(tmp$name_CDS))



genes <- list()
ngenes <- 0
for(i in 1:nrow(tmp)) { # for(i in 2:nrow(tmp)) {
  tryCatch({
      gb <- read_gbk(tmp$gbk_filename[i])
      ans <- gbk_features_to_df(gb)
      for(g in query_genes) {
        if(g %in% ans$gene) {
          genes[[g]] <- ans[which(ans$gene == g), ]$translation
          ngenes <- ngenes + 1
          if(ngenes == length(query_genes)) {
            break
          }
        }      
      }
      if(ngenes == length(query_genes)) {
        break
      }
    }, 
    error=function(e) {
      message(e)
    },
    finally = {
      next
    })
}

ans <- lapply(1:length(genes), function(i) AAString(genes[[i]]))
names(ans) <- names(genes)
ans <- AAStringSet(ans)
#res <- getSeq(gb@sequence, setNames(GENES, GENES$gene_id))
writeXStringSet(ans, paste0("gene_grouper/", bgc_type, "/query_genes.fna"))


# Convert gbk to gbff:
gbff_dir <- paste0("gene_grouper/", bgc_type, "/genomes_gbff/")
if(!dir.exists(gbff_dir)) {
  dir.create(gbff_dir, recursive = T)
}

com = paste("python", 
            "convert_to_gbff.py", 
            paste("-i", gbk_dir),
            paste("-o", gbff_dir))
system(com)

# Build gene grouper database



