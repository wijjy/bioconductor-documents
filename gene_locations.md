## Getting gene locations

```
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
sym <- as.list(org.Hs.egSYMBOL2EG)
id2symbol <- rep(names(sym),sapply(sym,length))
names(id2symbol) <- unlist(sym)

## we need to find the transcripts of these genes
allgenes <- transcriptsBy(txdb, "gene")
allgenes_df <- as.data.frame(reduce(allgenes))

row_match <- match(allgenes_df$group_name, sym)
allgenes_df <- cbind(allgenes_df, name = names(sym)[row_match])
genesGR <- GRanges(seqnames=Rle(allgenes_df$seqnames),
                   ranges = IRanges(allgenes_df$start, allgenes_df$end,
                                    names=allgenes_df$name), strand=allgenes_df$strand)

rm(allgenes,allgenes_df,sym, txdb, row_match, id2symbol)
genome(genesGR) <- "hg19"
save(genesGR, file="genesGR.rda")
```
