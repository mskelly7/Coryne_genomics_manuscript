### Potential to encode siderophore among respiratory Corynebacterium species
library(dplyr)
library(ggplot2)

dat <- read.csv("/Users/iz12/Library/CloudStorage/Box-Box/Coryne Genomics/Results/protocores_respiratory.csv")
dat <- dat[which(!duplicated(dat$strain_id)),]

species <- c("propinquum", "pseudodiphtheriticum", "accolens")

dat$final_species2 <- ifelse(dat$final_species %in% species, species, "Other")
dat$Cluster2 <- ifelse(dat$product_protoCore == "siderophore", "siderophore", "Other")
dat$final_species2 <- ifelse(dat$final_species2 == "propinquum", "C. propinquum", dat$final_species2)
dat$final_species2 <- ifelse(dat$final_species2 == "pseudodiphtheriticum", "C. pseudodiphtheriticum", dat$final_species2)
dat$final_species2 <- ifelse(dat$final_species2 == "accolens", "C. accolens", dat$final_species2)


ans <- as.data.frame(dat[, c("Cluster2", "final_species2", "product_protoCore")])
colnames(ans) <- c("Cluster", "Species", "product_protoCore")
ans <- ans %>% group_by(Species, Cluster, product_protoCore) %>% summarise(n=n())
ans <- ans[which(ans$n >= 1 & ans$Species != "Other"), ]
ans <- ans %>% group_by(Species, Cluster, n)

p <- ggplot(data=ans, aes(x=Species, y=n, fill=Cluster)) + geom_bar(stat="identity") + theme_bw() 
p <- p + ylab("# of species encoding siderophore BGCs") 
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
               panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_fill_manual(name="", values = c("siderophore" = "lightgreen", "Other" = "deepskyblue"), labels = c("Absence", "Presence"))
p
