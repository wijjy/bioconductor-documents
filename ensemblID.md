Problem
=====

Introduction
----------------

I have data associated with a set of Ensembl IDs like this


Ensembl ID	| baseMean	| FoldChange|	gene_symbol
-------------------|------------------|-----------------|	
ENSMUSG00000083414	|251.9660149612	| 0.20|	Gm14840	
ENSMUSG00000053024	|321.8617500455	| 0.21| Cntn2	
ENSMUSG00000058625	|288.034261079	|	0.2420134298|Gm17383	
ENSMUSG00000100720	|47.9402007709|	0.1550504607	|Gm10075	

I also have a set of regions in mm9.

chromosome|start|end
--|--|--
chr1|	9956201	|9959400|
chr1|	12823301	|12827600
chr1	|34184301|	34188400
chr1	|38493701|	38497600


I want to know something about the relationship between the quantitative data and the locations.  (There is other data and other information, the problem is not so open ended).  

## Question 1: How do I get locations for the genes in mm9 from the ensembl ID?

### Packages?

The `org.Mm.eg.db` package has a mapping between Ensembl IDs and other sorts of ID.

This package is describes as

> Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers.

This package allows you to look up some information based on the Ensembl ID.

All this stuff seems to be included in the `Mus.musculus` package.  This is a package which
abstracts a number of _Mus musculus_ packages within bioconductor.  Information about
functionality is in the [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) documentation, which includes a 
[video](https://www.youtube.com/watch?v=8qvGNTVz3Ik). 

The `Mus.musculus` package removes some of the work done at the end of this file as it seems to automate the 
work that I was doing by hand.  


```
require(Mus.musculus)
columns(Mus.musculus)
```

Output:

```
 [1] "ACCNUM"       "ALIAS"        "CDSCHROM"     "CDSEND"       "CDSID"        "CDSNAME"      "CDSSTART"     "CDSSTRAND"   
 [9] "DEFINITION"   "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
[17] "EXONCHROM"    "EXONEND"      "EXONID"       "EXONNAME"     "EXONRANK"     "EXONSTART"    "EXONSTRAND"   "GENEID"      
[25] "GENENAME"     "GO"           "GOALL"        "GOID"         "IPI"          "MGI"          "ONTOLOGY"     "ONTOLOGYALL" 
[33] "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "TERM"         "TXCHROM"     
[41] "TXEND"        "TXID"         "TXNAME"       "TXSTART"      "TXSTRAND"     "TXTYPE"       "UNIGENE"      "UNIPROT"  
```

And we can find out those that can be used as keys

```
keytypes(Mus.musculus)
 [1] "ACCNUM"       "ALIAS"        "CDSID"        "CDSNAME"      "DEFINITION"   "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
 [9] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "EXONID"       "EXONNAME"     "GENEID"       "GENENAME"    
[17] "GO"           "GOALL"        "GOID"         "IPI"          "MGI"          "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
[25] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "TERM"         "TXID"         "TXNAME"      
[33] "UNIGENE"      "UNIPROT"     
```

This is clearly the way to go.

So my approach:

### Select the genes using SYMBOL (can also be sone with ensembl ID)

```{r}
res.down <- select(Mus.musculus, keys=paste(down$gene_symbol), columns=c("CDSCHROM", "CDSSTART", "CDSEND", "CDSSTRAND","REFSEQ","GENEID", "GENENAME"), keytype="SYMBOL")
```

### Turn these to GRanges

```{r}
mus.down <- GRanges(seqnames=res.down$CDSCHROM, IRanges(start=res.down$CDSSTART, end=res.down$CDSEND, name=res.down$SYMBOL), strand=res.down$CDSSTRAND, geneid=res.down$GENEID)
```

### Find overlaps between these and the target

```{r}
overDown <- findOverlaps(promoters(mus.down), regions10)
overSYMBOL <- names(mus.down)[queryHits(overDown)]
```

### And get the results for those keys

```
target <- select(Mus.musculus, 
                 keys=paste(overSYMBOL), 
                 columns=c("GENEID", "GENENAME", "GO", "CDSCHROM", "CDSSTART", "CDSEND", "ENSEMBL", "PATH"), 
                 keytype="SYMBOL")
```


***



Older Libraries
================


The information we have from the experiment in terms of selecting genes
are the ENSEMBL ID's and the and the gene symbol.  Now both of these are available in
the `org.Mm.eg.db` library.   

### Transcript libraries

```
require("TxDb.Mmusculus.UCSC.mm9.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
txdb
keytypes(txdb)
```




```
columns(txdb)
keytypes(txdb)
```

The key column is going to be GENEID, which is the ENTREZ gene ID, which is also available in the `org.Mm.eg.db` library.

### Gene information

```
require(org.Mm.eg.db)
ls(2)
columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
ky <- keys(org.Mm.eg.db, keytype="ENTREZID")
```


```
s <- select(org.Mm.eg.db, 
       keys=paste(down$Ensembl.ID),
       columns=c("ENTREZID", "GENENAME"),
       keytype="ENSEMBL")
s <- s[!is.na(s$ENTREZID),]

```

Now we need to be able to merge tables

```
GR <- transcripts(txdb)
s2 <- select(txdb,
             keys=s$ENTREZID,
             columns=c("TXCHROM", "TXSTART", "TXEND","TXSTRAND"),
             keytype="GENEID")

s3 <- merge(s,s2,by.x="ENTREZID", by.y="GENEID")

## Now turn into sequences

```
s3 <- s3[complete.cases(s3),]
mus.down <- GRanges(seqnames=s3$TXCHROM, IRanges(start=s3$TXSTART, end=s3$TXEND, name=s3$ENSEMBL), strand=s3$TXSTRAND, geneid=s3$ENTREZID)


overDown <- findOverlaps(promoters(mus.down), regions)
overSYMBOL <- names(mus.down)[queryHits(overDown)]
overRegions <- regions[subjectHits(overDown)]

```



Now write the results out

```{r}
target <- select(Mus.musculus, 
                 keys=paste(overSYMBOL), 
                 columns=c("GENEID", "GENENAME", "GO", "TXCHROM", "TXSTART", "TXEND", "ENSEMBL", "PATH"), 
                 keytype="SYMBOL")
```


## Old Code That may Still be useful


```{r}
require(org.Mm.eg.db)
yy <- as.list(org.Mm.egENSEMBL2EG)
yy[["ENSMUSG00000053024"]]
```

Gives "21367", which is the entrez ID of this gene.  We can see this 
using

```{r}
x<- org.Mm.egGENENAME
x[["21367"]]
```
gives "contactin 2".

Now the `TxDb.Mmusculus.UCSC.mm9.knownGene` contains transcript information, which should allow me to get information about promoters as well as locations.

```{r}
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
ts <- transcripts(txdb, columns="gene_id")
```


