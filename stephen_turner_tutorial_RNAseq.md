
Test from [http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html](http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html)




#Tutorial: RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview

## BACKGROUND

This tutorial shows an example of RNA-seq data analysis with DESeq2,
followed by KEGG pathway analysis using GAGE. Using data from
GSE37704, with processed data available on Figshare
DOI: 10.6084/m9.figshare.1601975. This dataset has six samples from
GSE37704, where expression was quantified by either: (A) mapping to to
GRCh38 using STAR then counting reads mapped to genes with
featureCounts under the union-intersection model, or (B)
alignment-free quantification using Sailfish, summarized at the gene
level using the GRCh38 GTF file. Both datasets are restricted to
protein-coding genes only. Here I’ll use the Sailfish gene-level
estimated counts.

## DIFFERENTIAL EXPRESSION ANALYSIS

First, import the countdata and metadata directly from the web. Set up
the DESeqDataSet, run the DESeq2 pipeline.

```

# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)

# Which data do you want to use? Let's use the sailfish counts.
# browseURL("http://dx.doi.org/10.6084/m9.figshare.1601975")
# countDataURL = "http://files.figshare.com/2439061/GSE37704_featurecounts.csv"
countDataURL = "http://files.figshare.com/2600373/GSE37704_sailfish_genecounts.csv"

# Import countdata
countData = read.csv(countDataURL, row.names=1) %>% 
  dplyr::select(-length) %>% 
  as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)
##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
## ENSG00000198888     17528     23007     30241     24418     29152
## ENSG00000198763     21264     26720     35550     28878     32416
## ENSG00000198804    130975    151207    195514    178130    196727
## ENSG00000198712     49769     61906     78608     66478     69758
## ENSG00000228253      9304     11160     12830     12608     13041
## ENSG00000198899     45401     51260     66851     63433     66123
##                 SRR493371
## ENSG00000198888     34416
## ENSG00000198763     38422
## ENSG00000198804    244670
## ENSG00000198712     86808
## ENSG00000228253     16063
## ENSG00000198899     79215
# Import metadata
colData = read.csv("http://files.figshare.com/2439060/GSE37704_metadata.csv", row.names=1)
colData
##               condition
## SRR493366 control_sirna
## SRR493367 control_sirna
## SRR493368 control_sirna
## SRR493369      hoxa1_kd
## SRR493370      hoxa1_kd
## SRR493371      hoxa1_kd
# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~condition)
dds = DESeq(dds)
dds
## class: DESeqDataSet 
## dim: 16755 6 
## metadata(0):
## assays(3): counts mu cooks
## rownames(16755): ENSG00000198888 ENSG00000198763 ...
##   ENSG00000267795 ENSG00000165795
## rowRanges metadata column names(27): baseMean baseVar ... deviance
##   maxCooks
## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
## colData names(2): condition sizeFactor
```

Next, get results for the HoxA1 knockdown versus control siRNA, and
reorder them by p-value. Call summary on the results object to get a
sense of how many genes are up or down-regulated at FDR 0.1.

```
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)
## 
## out of 16755 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4193, 25% 
## LFC < 0 (down)   : 4286, 26% 
## outliers [1]     : 22, 0.13% 
## low counts [2]   : 1299, 7.8% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Since we mapped and counted against the Ensembl annotation, our
results only have information about Ensembl gene IDs. But, our pathway
analysis downstream will use KEGG pathways, and genes in KEGG pathways
are annotated with Entrez gene IDs. I wrote an R package for doing
this offline the dplyr way
[https://github.com/stephenturner/annotables](https://github.com/stephenturner/annotables),
but the canonical Bioconductor way to do it is with the AnnotationDbi
and organism annotation packages. Here we’re using the organism
package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi
database package (“db”) using Entrez Gene IDs (“eg”) as primary
keys. To see what all the keys are, use the columns function.

```
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
## [1] "ACCNUM" "ALIAS" "ENSEMBL" "ENSEMBLPROT"

##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
```

Let’s use the mapIds function to add more columns to the results. The
row.names of our results table has the Ensembl gene ID (our key), so
we need to specify keytype=ENSEMBL. The column argument tells the
mapIds function which information we want, and the multiVals argument
tells the function what to do if there are multiple possible values
for a single input value. Here we ask to just give us back the first
one that occurs in the database. Let’s get the Entrez IDs, gene
symbols, and full gene names.

```
res$symbol = mapIds(org.Hs.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL"
, multiVals="first")

res$entrez = mapIds(org.Hs.eg.db, keys=row.names(res)
, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

res$name = mapIds(org.Hs.eg.db, keys=row.names(res), column="GENENAME", keytype="ENSEMBL",
multiVals="first")

head(res, 10)
## log2 fold change (MAP): condition hoxa1_kd vs control_sirna 
## Wald test p-value: condition hoxa1_kd vs control_sirna 
## DataFrame with 10 rows and 9 columns
##                  baseMean log2FoldChange      lfcSE      stat    pvalue
##                           
## ENSG00000148773  1885.344      -3.172502 0.07868572 -40.31865         0
## ENSG00000138623  2939.936      -2.418238 0.05889229 -41.06205         0
## ENSG00000104368 13601.963       2.016802 0.05249643  38.41789         0
## ENSG00000124766  2692.200       2.379545 0.06193654  38.41908         0
## ENSG00000122861 35889.413       2.224779 0.05258658  42.30697         0
## ENSG00000116016  4558.157      -1.885339 0.04258766 -44.26961         0
## ENSG00000164251  2404.103       3.325196 0.07021236  47.35912         0
## ENSG00000125257  6187.386       1.943762 0.04259189  45.63692         0
## ENSG00000104321  9334.555       3.186856 0.06227530  51.17367         0
## ENSG00000183508  2110.345       3.190612 0.07488305  42.60794         0
##                      padj      symbol      entrez
##                   
## ENSG00000148773         0       MKI67        4288
## ENSG00000138623         0      SEMA7A        8482
## ENSG00000104368         0        PLAT        5327
## ENSG00000124766         0        SOX4        6659
## ENSG00000122861         0        PLAU        5328
## ENSG00000116016         0       EPAS1        2034
## ENSG00000164251         0       F2RL1        2150
## ENSG00000125257         0       ABCC4       10257
## ENSG00000104321         0       TRPA1        8989
## ENSG00000183508         0      FAM46C       54855
##                                                                               name
##                                                                        
## ENSG00000148773                                      marker of proliferation Ki-67
## ENSG00000138623 semaphorin 7A, GPI membrane anchor (John Milton Hagen blood group)
## ENSG00000104368                                      plasminogen activator, tissue
## ENSG00000124766                               SRY (sex determining region Y)-box 4
## ENSG00000122861                                   plasminogen activator, urokinase
## ENSG00000116016                                   endothelial PAS domain protein 1
## ENSG00000164251                   coagulation factor II (thrombin) receptor-like 1
## ENSG00000125257            ATP-binding cassette, sub-family C (CFTR/MRP), member 4
## ENSG00000104321 transient receptor potential cation channel, subfamily A, member 1
## ENSG00000183508                       family with sequence similarity 46, member C
```



## PATHWAY ANALYSIS

We’re going to use the gage package (Generally Applicable Gene-set
Enrichment for Pathway Analysis) for pathway analysis. See also the
gage package workflow vignette for RNA-seq pathway analysis. Once we
have a list of enriched pathways, we’re going to use the pathview
package to draw pathway diagrams, shading the molecules in the pathway
by their degree of up/down-regulation.  KEGG pathways

The gageData package has pre-compiled databases mapping genes to KEGG
pathways and GO terms for common organisms. kegg.sets.hs is a named
list of 229 elements. Each element is a character vector of member
gene Entrez IDs for a single KEGG pathway. (See also
go.sets.hs). sigmet.idx.hs is an index of numbers of sinaling and
metabolic pathways in kegg.set.gs. In other words, KEGG pathway
include other types of pathway definitions, like “Global Map” and
“Human Diseases”, which may be undesirable in pathway
analysis. Therefore, kegg.sets.hs[sigmet.idx.hs] gives you the
“cleaner” gene sets of sinaling and metabolic pathways only.

```
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
## $`hsa00232 Caffeine metabolism`
## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
## 
## $`hsa00983 Drug metabolism - other enzymes`
##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
## [50] "8833"   "9"      "978"   
## 
## $`hsa00230 Purine metabolism`
##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
## [162] "9615"
```

The gage() function requires a named vector of fold changes, where the
names of the values are the Entrez gene IDs.

```
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
##      4288      8482      5327      6659      5328      2034 
## -3.172502 -2.418238  2.016802  2.379545  2.224779 -1.885339
```


Now, let’s run the pathway analysis. See help on the gage function
with ?gage. Specifically, you might want to try changing the value of
s`ame.dir. This value determins whether to test for changes in a gene
set toward a single direction (all genes up or down regulated) or
changes towards both directions simultaneously (any genes in the
pathway dysregulated).

For experimentally derived gene sets, GO term groups, etc,
coregulation is commonly the case, hence same.dir = TRUE (default); In
KEGG, BioCarta pathways, genes frequently are not coregulated, hence
it could be informative to let same.dir = FALSE. Although same.dir =
TRUE could also be interesting for pathways.

Here, we’re using same.dir = TRUE, which will give us separate lists
for pathways that are upregulated versus pathways that are
downregulated. Let’s look at the first few results from each.

```
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)
## $greater
##                                          p.geomean stat.mean        p.val
## hsa04142 Lysosome                     0.0002630657  3.517890 0.0002630657
## hsa04640 Hematopoietic cell lineage   0.0017919390  2.976432 0.0017919390
## hsa04630 Jak-STAT signaling pathway   0.0048980977  2.604390 0.0048980977
## hsa00140 Steroid hormone biosynthesis 0.0051115493  2.636206 0.0051115493
## hsa04062 Chemokine signaling pathway  0.0125582961  2.250765 0.0125582961
## hsa00511 Other glycan degradation     0.0223819919  2.104311 0.0223819919
##                                            q.val set.size         exp1
## hsa04142 Lysosome                     0.04261664      116 0.0002630657
## hsa04640 Hematopoietic cell lineage   0.14514706       61 0.0017919390
## hsa04630 Jak-STAT signaling pathway   0.20701775      119 0.0048980977
## hsa00140 Steroid hormone biosynthesis 0.20701775       39 0.0051115493
## hsa04062 Chemokine signaling pathway  0.40688879      156 0.0125582961
## hsa00511 Other glycan degradation     0.49956506       15 0.0223819919
## 
## $less
##                                      p.geomean stat.mean        p.val
## hsa04110 Cell cycle               2.165725e-06 -4.722301 2.165725e-06
## hsa03030 DNA replication          3.807440e-06 -4.835336 3.807440e-06
## hsa04114 Oocyte meiosis           1.109869e-04 -3.767561 1.109869e-04
## hsa03013 RNA transport            1.181787e-03 -3.071947 1.181787e-03
## hsa03440 Homologous recombination 1.197124e-03 -3.190747 1.197124e-03
## hsa00240 Pyrimidine metabolism    1.570318e-03 -2.992059 1.570318e-03
##                                          q.val set.size         exp1
## hsa04110 Cell cycle               0.0003084027      121 2.165725e-06
## hsa03030 DNA replication          0.0003084027       36 3.807440e-06
## hsa04114 Oocyte meiosis           0.0059932916      101 1.109869e-04
## hsa03013 RNA transport            0.0387868193      145 1.181787e-03
## hsa03440 Homologous recombination 0.0387868193       28 1.197124e-03
## hsa00240 Pyrimidine metabolism    0.0423985796       96 1.570318e-03
## 
## $stats
##                                       stat.mean     exp1
## hsa04142 Lysosome                      3.517890 3.517890
## hsa04640 Hematopoietic cell lineage    2.976432 2.976432
## hsa04630 Jak-STAT signaling pathway    2.604390 2.604390
## hsa00140 Steroid hormone biosynthesis  2.636206 2.636206
## hsa04062 Chemokine signaling pathway   2.250765 2.250765
## hsa00511 Other glycan degradation      2.104311 2.104311
```

Now, let’s process the results to pull out the top 5 upregulated
pathways, then further process that just to get the IDs. We’ll use
these KEGG pathway IDs downstream for plotting.

```
# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways
## [1] "hsa04142 Lysosome"                    
## [2] "hsa04640 Hematopoietic cell lineage"  
## [3] "hsa04630 Jak-STAT signaling pathway"  
## [4] "hsa00140 Steroid hormone biosynthesis"
## [5] "hsa04062 Chemokine signaling pathway"
# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
## [1] "hsa04142" "hsa04640" "hsa04630" "hsa00140" "hsa04062"
Finally, the pathview() function in the pathview package makes the plots. Let’s write a function so we can loop through and draw plots for the top 5 pathways we created above.
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
```

Here are the plots:





## Gene Ontology (GO)

We can also do a similar procedure with gene ontology. Similar to
above, go.sets.hs has all GO terms. go.subs.hs is a named list
containing indexes for the BP, CC, and MF ontologies. Let’s only do
Biological Process.

```
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)
## $greater
##                                                             p.geomean
## GO:0007156 homophilic cell adhesion                      3.914568e-05
## GO:0008285 negative regulation of cell proliferation     2.907332e-04
## GO:0016339 calcium-dependent cell-cell adhesion          4.218753e-04
## GO:0016337 cell-cell adhesion                            6.170551e-04
## GO:0048729 tissue morphogenesis                          6.581460e-04
## GO:1901617 organic hydroxy compound biosynthetic process 8.876161e-04
##                                                          stat.mean
## GO:0007156 homophilic cell adhesion                       4.017207
## GO:0008285 negative regulation of cell proliferation      3.453345
## GO:0016339 calcium-dependent cell-cell adhesion           3.543891
## GO:0016337 cell-cell adhesion                             3.244296
## GO:0048729 tissue morphogenesis                           3.223979
## GO:1901617 organic hydroxy compound biosynthetic process  3.157421
##                                                                 p.val
## GO:0007156 homophilic cell adhesion                      3.914568e-05
## GO:0008285 negative regulation of cell proliferation     2.907332e-04
## GO:0016339 calcium-dependent cell-cell adhesion          4.218753e-04
## GO:0016337 cell-cell adhesion                            6.170551e-04
## GO:0048729 tissue morphogenesis                          6.581460e-04
## GO:1901617 organic hydroxy compound biosynthetic process 8.876161e-04
##                                                              q.val
## GO:0007156 homophilic cell adhesion                      0.1613977
## GO:0008285 negative regulation of cell proliferation     0.4720349
## GO:0016339 calcium-dependent cell-cell adhesion          0.4720349
## GO:0016337 cell-cell adhesion                            0.4720349
## GO:0048729 tissue morphogenesis                          0.4720349
## GO:1901617 organic hydroxy compound biosynthetic process 0.4720349
##                                                          set.size
## GO:0007156 homophilic cell adhesion                           124
## GO:0008285 negative regulation of cell proliferation          458
## GO:0016339 calcium-dependent cell-cell adhesion                27
## GO:0016337 cell-cell adhesion                                 355
## GO:0048729 tissue morphogenesis                               429
## GO:1901617 organic hydroxy compound biosynthetic process      141
##                                                                  exp1
## GO:0007156 homophilic cell adhesion                      3.914568e-05
## GO:0008285 negative regulation of cell proliferation     2.907332e-04
## GO:0016339 calcium-dependent cell-cell adhesion          4.218753e-04
## GO:0016337 cell-cell adhesion                            6.170551e-04
## GO:0048729 tissue morphogenesis                          6.581460e-04
## GO:1901617 organic hydroxy compound biosynthetic process 8.876161e-04
## 
## $less
##                                             p.geomean stat.mean
## GO:0048285 organelle fission             4.411540e-18 -8.850004
## GO:0000280 nuclear division              7.459684e-18 -8.805564
## GO:0007067 mitosis                       7.459684e-18 -8.805564
## GO:0000087 M phase of mitotic cell cycle 2.286444e-17 -8.655644
## GO:0007059 chromosome segregation        1.872901e-13 -7.686883
## GO:0051301 cell division                 5.841375e-12 -6.887763
##                                                 p.val        q.val
## GO:0048285 organelle fission             4.411540e-18 1.025209e-14
## GO:0000280 nuclear division              7.459684e-18 1.025209e-14
## GO:0007067 mitosis                       7.459684e-18 1.025209e-14
## GO:0000087 M phase of mitotic cell cycle 2.286444e-17 2.356752e-14
## GO:0007059 chromosome segregation        1.872901e-13 1.544394e-10
## GO:0051301 cell division                 5.841375e-12 4.013998e-09
##                                          set.size         exp1
## GO:0048285 organelle fission                  376 4.411540e-18
## GO:0000280 nuclear division                   352 7.459684e-18
## GO:0007067 mitosis                            352 7.459684e-18
## GO:0000087 M phase of mitotic cell cycle      362 2.286444e-17
## GO:0007059 chromosome segregation             141 1.872901e-13
## GO:0051301 cell division                      462 5.841375e-12
## 
## $stats
##                                                          stat.mean
## GO:0007156 homophilic cell adhesion                       4.017207
## GO:0008285 negative regulation of cell proliferation      3.453345
## GO:0016339 calcium-dependent cell-cell adhesion           3.543891
## GO:0016337 cell-cell adhesion                             3.244296
## GO:0048729 tissue morphogenesis                           3.223979
## GO:1901617 organic hydroxy compound biosynthetic process  3.157421
##                                                              exp1
## GO:0007156 homophilic cell adhesion                      4.017207
## GO:0008285 negative regulation of cell proliferation     3.453345
## GO:0016339 calcium-dependent cell-cell adhesion          3.543891
## GO:0016337 cell-cell adhesion                            3.244296
## GO:0048729 tissue morphogenesis                          3.223979
## GO:1901617 organic hydroxy compound biosynthetic process 3.157421
```
