#### Fine mapping AS association on 2p15

library(data.table)
library(ggplot2)
library(corrcoverage)
library(cowplot)

x <- fread("zcat hg19_gwas_ic_as_igas_4_19_1.tab.gz")

tagSNP <- x[which(x$Marker=="rs6759298"),] 

snpChr <- tagSNP$Chr

snpPos <- tagSNP$Position

X <- readRDS(paste0("geno-",snpChr,".RDS"))

lowr <- snpPos-50000

uppr <- snpPos+50000

# Find all snps in this region

head(X@snps)

snp_indices <- which(X@snps$chromosome==snpChr & dplyr::between(X@snps$position, lowr, uppr))

snpdata <- X@snps[snp_indices,]

# Match SNPs by positon and drop those that cannot be matched in immunochip data

tokeep <- snpdata$position %in% x$Position

snpdata_matched <- snpdata[tokeep,]

snp_indices_matched <- snp_indices[tokeep]

# all(snpdata_matched$position %in% x$Position)

# Find SNP correlations

Xsub <- X[,snp_indices_matched]

Sigma <- ld(sm(Xsub), sm(Xsub), stat="R") # note this is r not r^2

colnames(Sigma) <- snpdata_matched$snp.name
rownames(Sigma) <- snpdata_matched$snp.name

# Make data frame with summary stats

snpdata_matched$Pvalue <- x$PValue[match(snpdata_matched$position, x$Position)]

snpdata_matched$OR <- x[match(snpdata_matched$position, x$Position), 5]

snpdata_matched$LowerOR <- x$LowerOR[match(snpdata_matched$position, x$Position)]

snpdata_matched$UpperOR <- x$UpperOR[match(snpdata_matched$position, x$Position)]

# Find SNPs with r2>=0.8 with lead SNP

snpdata_matched[which(snpdata_matched$position==snpPos),]

which(Sigma[which(snpdata_matched$position==snpPos),]^2>=0.8)

####### FINE-MAPPING

# Need estimated effect sizes (betas) and SEs

# Convert P vals to Z scores

snpdata_matched$Z <- qnorm(snpdata_matched$Pvalue/2,lower.tail=FALSE)

snpdata_matched[,8] <- as.numeric(unlist(snpdata_matched[,8]))

snpdata_matched$beta <- log(snpdata_matched[,8])

snpdata_matched$SE <- snpdata_matched$beta/snpdata_matched$Z

snpdata_matched$PP <- ppfunc(snpdata_matched$Z, snpdata_matched$SE^2)

snpdata_matched$r2 <- Sigma[which(snpdata_matched$position==snpPos),]^2

# write.csv(snpdata_matched, "AS_2p15.csv")

betas <- snpdata_matched$beta

names(betas) <- snpdata_matched$snp.name

orig_99 <- credset(snpdata_matched$PP, thr = 0.99)

corrected_99cs <- corrected_cs_bhat(betas, snpdata_matched$SE^2, N0 = 15145, N1 = 10619, Sigma = Sigma, desired.cov = 0.99)

which( snp_matched$snp.name == "rs6759298") # lead SNP
which( snp_matched$snp.name == "rs13001372") # SNP mentioned in manuscript

p.plot <- ggplot(snpdata_matched, aes(x=position, y=-log10(Pvalue))) +
  geom_point( aes(color="black"), alpha=0.2, size=2.5) + theme_bw() +
  scale_color_manual(values = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+labs(x = "Position",y = "-log10(P)") +theme(plot.title = element_text(color="black", size=10, face="bold", hjust=0.5))+theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.line.x = element_line(color="black", size = 0.4), axis.line.y = element_line(color="black", size = 0.4)) + ggtitle("Manhattan Plot") + geom_point(data = snpdata_matched[which(snpdata_matched$snp.name=="rs6759298"),], aes(x=position, y=-log10(Pvalue)), color="red", size=2.5) + geom_point(data = snpdata_matched[which(snpdata_matched$snp.name=="rs13001372"),], aes(x=position, y=-log10(Pvalue)), color="blue", size=2.5)

pp.plot <- ggplot(snpdata_matched, aes(x=position, y=PP)) +
  geom_point( aes(color="black"), alpha=0.2, size=2.5) + theme_bw() +
  scale_color_manual(values = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+labs(x = "Position",y = "Posterior Probability") +theme(plot.title = element_text(color="black", size=10, face="bold", hjust=0.5))+theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.line.x = element_line(color="black", size = 0.4),axis.line.y = element_line(color="black", size = 0.4)) + ggtitle("Posterior Probability Plot") + geom_point(data = snpdata_matched[which(snpdata_matched$snp.name=="rs6759298"),], aes(x=position, y=PP), color="red", size=2.5) + geom_point(data = snpdata_matched[which(snpdata_matched$snp.name=="rs13001372"),], aes(x=position, y=PP), color="blue", size=2.5)

plot_grid(p.plot, pp.plot, nrow = 1)