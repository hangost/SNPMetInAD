library(VariantAnnotation)
library(GenomicFeatures)
library(doParallel)
library(BiocParallel)
library(data.table)

change.nm <- function(temp.mat,te.mat,id.nm,temp.nm,te.nm,te1="[.]"){
    rownames(temp.mat) <- temp.mat[,temp.nm]
    colnames(te.mat) <- gsub(te1,"-",colnames(te.mat))
    over.nm <- intersect(colnames(te.mat)[te.nm],temp.mat[,temp.nm])
    cn <- c(id.nm,over.nm)
    if (class(te.mat)[1] == "data.table"){
        f.mat <- te.mat[,..cn]
    }
    else {
        f.mat <- te.mat[,cn]
    }
    colnames(f.mat)[colnames(f.mat) != id.nm] <- temp.mat[over.nm,"wgs_id"]
    return (f.mat)
}

er.te <- function(te.mat,int.te=NULL){
    if (length(grep("Error",te.mat)))    return ("NA")
    if (!length(int.te))    return (te.mat)
    te.p <- try(summary(te.mat)$coefficient[int.te,4],silent=T)
    if (length(grep("Error",te.p)))    return (NA)
    return (te.p)
}

rm.out <- function(te.mat){
    ea.mat <- te.mat
    fi.q <- summary(as.double(ea.mat))[2]
    se.q <- summary(as.double(ea.mat))[5]
    int.q <- se.q - fi.q
    out.v <- c(fi.q - (3 * int.q),se.q + (3 * int.q))
    re.nms <- c(which(as.double(ea.mat) > out.v[1] & as.double(ea.mat) < out.v[2]))
    return (re.nms)
}

vcf.names <- c("NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_",".recalibrated_variants.Broad_Rush.vcf.gz")

exp.data <- gsub(" ","",as.matrix(read.table("./4.rsem_result/ROSMAP_exp",sep='\t',header=T)))
medi.exp <- apply(exp.data[,-1],1,function(x)    median(as.double(x)))
exp.data <- exp.data[which(medi.exp > 0.1),]

methyl.exp <- fread("./ROSMAP_arrayMethylation_imputed.tsv")
methyl.info <- gsub(" ","",as.matrix(read.csv("./ROSMAP_arrayMethylation_metaData.tsv",sep='\t',header=T,comment.char="")))
methyl.info <- methyl.info[methyl.info[,"CHR"] != "" & methyl.info[,"MAPINFO"] != "",]
methyl.ran <- GRanges(Rle(methyl.info[,"CHR"]),IRanges(as.integer(methyl.info[,"MAPINFO"]),as.integer(methyl.info[,"MAPINFO"])))

colnames(exp.data) <- gsub("ROSMAP_","",colnames(exp.data))
ROSMAP.cli <- gsub(" ","",as.matrix(read.table("./ROSMAP_clinical.csv",sep=',',header=T)))
ROSMAP.ID <- gsub(" ","",as.matrix(read.table("./ROSMAP_IDkey.csv",sep=',',header=T)))
ROSMAP.mer.cli <- as.matrix(merge(ROSMAP.cli,ROSMAP.ID,by.x="projid",by.y="projid"))
ROSMAP.mer.cli <- ROSMAP.mer.cli[ROSMAP.mer.cli[,"rnaseq_id"] != "" & ROSMAP.mer.cli[,"wgs_id"] != "" & ROSMAP.mer.cli[,"mwas_id"] != "",]
f.cli.mat <- ROSMAP.mer.cli

AD.sam <- f.cli.mat[as.integer(f.cli.mat[,"ceradsc"]) < 4,"wgs_id"]
CN.sam <- f.cli.mat[as.integer(f.cli.mat[,"ceradsc"]) == 4,"wgs_id"]

# change sample ids for expression and methylation into WGS ids.
te.exp.data <- change.nm(f.cli.mat,exp.data,"TXID","rnaseq_id",c(2:637))
te.methyl.exp <- change.nm(f.cli.mat,methyl.exp,"TargetID","mwas_id",c(2:741))

GTFdb <- loadDb("~/txDB_GTF75")
promo.re <- promoters(GTFdb,2000,200)
promo.re <- promo.re[exp.data[,"TXID"],]
promo.mat <- cbind(as.matrix(elementMetadata(promo.re))[,2],as.character(seqnames(promo.re)),start(promo.re),end(promo.re))
over.re <- as.matrix(findOverlaps(promo.re,methyl.ran))
u.over.re <- unique(over.re[,1])


paral.re <- foreach(i=1:length(u.over.re),.combine=rbind ,.errorhandling = "pass") %dopar% {
    ea.re <- rbind(over.re[over.re[,1] == u.over.re[i],])
    tx.id <- promo.mat[as.integer(u.over.re[i]),1]
    methyl.id <- methyl.info[as.integer(ea.re[,2]),"TargetID"]

    tx.exp <- rbind(te.exp.data[te.exp.data[,1] == tx.id,])
    ea.pro.re <- promo.mat[promo.mat[,1] == tx.id,]

    params <- ScanVcfParam(which=GRanges(Rle(ea.pro.re[2]),IRanges(as.integer(ea.pro.re[3]),as.integer(ea.pro.re[4]))))
    vcf.file.path <- paste("./6.WGS/",vcf.names[1],ea.pro.re[2],vcf.names[2],sep="")
    geno.data <- rbind(geno(readVcf(TabixFile(vcf.file.path), "hg19", params))$GT)
    te.geno.data <- rbind(geno.data[,is.element(colnames(geno.data),intersect(colnames(te.methyl.exp),colnames(tx.exp)))])
    AD.geno.dt <- rbind(te.geno.data[,is.element(colnames(te.geno.data),AD.sam)])
    maf.f <- apply(AD.geno.dt,1,function(x){
        t.x <- table(unlist(strsplit(x,"/")))
        sum(na.omit(t.x[c("1","2")]))/(length(x[x != "./."]) * 2)
    })
    com.mu <- names(maf.f[which(maf.f > 0.05 & maf.f < 0.95)])
    te.snp.nm <- grep("_[A,C,G,T]/[A,C,G,T]$",com.mu)








registerDoParallel(cores=12)



me.sams <- colnames(te.methyl.exp)








