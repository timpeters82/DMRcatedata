### =========================================================================
### DMRcatedata metadata 
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("crosshyb", "snpsall", "XY.probes", "hg19.generanges", "hg19.grt",
            "hg38.generanges", "hg38.grt", "mm10.generanges", "mm10.grt"),
  Description = c(paste0("Character vector of Illumina probes whose ", 
                         "probe sequence promiscuously aligns to non-target ",
                          "sections of the genome matching 47bp or higher."), 
                  paste0("Data.frame of Illumina probes whose target ",
                         "CpG lies on or near a known SNP. SNP ID, distance to ",
                         "CpG and minor allele frequency are all reported here."), 
                  paste0("Character vector of Illumina probes whose target CpG ",
                         "site is on a sex chromosome."),
                  paste0("GRanges object giving the genomic intervals of all gene ",
                         "regions in the Ensembl Release 75 of hg19."), 
                  "GeneRegionTrack formulated from TxDb.Hsapiens.UCSC.hg19.knownGene.",
                  paste0("GRanges object giving the genomic intervals of all gene ",
                         "regions in the Ensembl Release 102 of hg38."),
                  "GeneRegionTrack formulated from TxDb.Hsapiens.UCSC.hg38.knownGene.",
                  paste0("GRanges object giving the genomic intervals of all gene ",
                         "regions in the Ensembl Release 102 of mm10."),
                  "GeneRegionTrack formulated from TxDb.Mmusculus.UCSC.mm10.knownGene."),
  BiocVersion = rep("3.10", 9),
  Genome = c(rep("hg19", 5), rep("hg38", 2), rep("mm10", 2)),
  SourceType = c(rep("CSV", 3), rep("GTF", 6)), 
  SourceUrl = c(paste0("https://tinyurl.com/y339e2lh, ", "https://tinyurl.com/yxm6c7g6, ", "https://tinyurl.com/y33fdy2w"),
                paste0("https://tinyurl.com/yxfco4sz, ", "https://tinyurl.com/y6lwcm9v, ", "https://tinyurl.com/y5dykrzh, ", "https://tinyurl.com/y54dfss8"),
                paste0("https://bitbucket.org/hansenlab/illuminahumanmethylationepicanno.ilm10b4.hg19/src/master/, ", 
                       "https://bitbucket.org/hansenlab/illuminahumanmethylation450kanno.ilmn12.hg19/src/master/"),
                rep("ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz", 2),
                rep("ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr.gtf.gz", 2),
                rep("ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.chr.gtf.gz", 2)),
  SourceVersion = rep("Feb 04 2021", 9),
  Species = c(rep("Homo sapiens", 7), rep("Mus musculus", 2)),
  TaxonomyId = c(rep(9606, 7), rep(10090, 2)),
  Coordinate_1_based = TRUE,
  DataProvider = c("SickKids, Springer", "Springer", "Bioconductor", rep("Ensembl", 6)),
  Maintainer = "Tim Peters <t.peters@garvan.org.au>",
  RDataClass = c("character", "data.frame", "character", rep(c("GRanges", "GeneRegionTrack"), times=3)),
  DispatchClass = c(rep("Rda",9)),
  RDataPath = paste("DMRcatedata", c("crosshyb.Rda", "snpsall.Rda", "XY.probes.Rda", "hg19.generanges.Rda", "hg19.grt.Rda", 
                "hg38.generanges.Rda", "hg38.grt.Rda", "mm10.generanges.Rda", "mm10.grt.Rda"), sep="/")
)
write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
