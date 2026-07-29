# Bambu annotation for Fullscope long-read alignments.
suppressPackageStartupMessages({
    library(bambu)
    library(qs)
})
OuttablefromBambu = function(readtrans,sample){
    readtransf = readtrans[(lengths(readtrans$compatibleMatches) + lengths(readtrans$equalMatches)) > 0,]
    nrow(readtransf)
    rtnume = lengths(readtransf$equalMatches)
    rtnumc = lengths(readtransf$compatibleMatches)
    readid = rep(readtransf$readId,rtnume)
    tranid = unlist(readtransf$equalMatches)

    readtransf.dfe = data.frame("readid" = readid,
                              "tranid" = tranid,
                               "type" = "equal")
    readid = rep(readtransf$readId,rtnumc)
    tranid = unlist(readtransf$compatibleMatches)

    readtransf.dfc = data.frame("readid" = readid,
                              "tranid" = tranid,
                               "type" = "compatible")
    readtransf.df = rbind(readtransf.dfe,readtransf.dfc)
    readtransf.df = readtransf.df[!duplicated(readtransf.df[,1:3]),]
    readtransf.df$sample = sample
    return(readtransf.df)
}

BambuMatrixBuild = function(bamls,bambuAnnotations,genome,outfile,samples){
    se <- bambu(reads = c(bamls),
                annotations = bambuAnnotations,
                genome = genome,discovery = FALSE,
                NDR = 0,
               trackReads = TRUE)
    se
    qsave(se,file = paste0(outfile,"_se.qs"))

    readtransf.df = NULL
    for(i in 1:length(metadata(se)$readToTranscriptMaps)){
        readtransf.dfi = OuttablefromBambu(metadata(se)$readToTranscriptMaps[[i]],samples[i])
        readtransf.df = rbind(readtransf.df,readtransf.dfi)
    }
    sean = as.data.frame(rowData(se))
    sean = sean[,c("TXNAME","GENEID","txid")]
    colnames(sean) = c("transcript_id","gene_id","txid")
    matchid = match(readtransf.df$tranid,sean$txid)
    readtransf.df$transcript_id = sean$transcript_id[matchid]
    readtransf.df$gene_id = sean$gene_id[matchid]
    print(length(unique(readtransf.df$readid)))
    qsave(readtransf.df,file = paste0(outfile,"_trans_total_anno.qs"))
    return(readtransf.df)
}
print("start to run")

# 从命令行接收参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) < 5) {
  stop("Usage: Rscript bambu_process.R annotations.gtf genome.fa input.bam output_prefix sample")
}

# 解析参数
gtf <- args[1]
genome <- args[2]
bamls <- args[3]
outfile <- args[4]
samples <- args[5]  # 将逗号分隔的字符串转换为向量

bambuAnnotations <- prepareAnnotations(gtf)
readtransf.df = BambuMatrixBuild(bamls,bambuAnnotations,genome,outfile,samples)
