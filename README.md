# RNAseq
Analysis of RNA seq data. From the form:

This run was for RNA isolated from new muscle biopsies of human breast cancer patients.  I want to have these included the group and have everything analyzed together as a whole.  We can simply run non-cancer controls versus all cancer, including the previous data you generated.


## Steps:

### Quality Control
In the reads folder, create a fastqc folder. Then run:

```
for f in *.gz; do echo $f; fastq -o fastqc -t 40 $f; done
```

This makes QC reports for each file. To get a single report, cd to the fastqc folder and use:
```
multiqc -o multiqc .
```

Quality looks good. There are two files called "test" that are different. (They are Ryan testing something else, and can safely be ignored.) I will exclude them from further analysis.


### Quantification

Get the reference from Ensembl, release 94. Got cDNA_all and ncRNA, cat them together and build salmon index.

```
conda activate salmon

# Build index
salmon index -t CDNA_and_NCrna.fa -i hg38.94 --type quasi -k 31

# Do quantification
for f in P*R1*; do echo $f; salmon quant -i ../Data/hg38.94 --libType A --gcBias -p 20 --numBootstraps 50 -o ../salmon/${f%_S*} -1 $f -2 ${f/R1/R2}; done

salmon version
# salmon v0.11.0

conda deativate
```


In the reads folder, use salmon to quantify each read

Adding in previously sequenced samples. Redoing quantification, to make sure everything is done equally. From Apr 2017 folder, run salmon again, writting to the current salmon directory.




### Differential expression

Use the R script doDESeq.R to find differentially expressed genes, and produce some figures

<<<<<<< HEAD
The set up for different prefixes:

##### All
```
PCA_Group <- 'Cancer'
design =~ Run + Surgeon + Cancer
contrast <- c('Cancer','Cancer','Control')
meta <- metadata
```

##### ERPR
```
outPrefix <- 'ERPR'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','ERPR','Control')
meta <- metadata %>% dplyr::filter(Group=='ERPR' | Group=='Control')
```

##### HER
```
outPrefix <- 'HER'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','HER','Control')
meta <- metadata %>% dplyr::filter(Group=='HER' | Group=='Control')
```

##### TN
```
outPrefix <- 'TN'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','TN','Control')
meta <- metadata %>% dplyr::filter(Group=='TN' | Group=='Control')
```

##### TP
```
outPrefix <- 'TP'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','TP','Control')
outPrefix <- 'TP'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','TP','Control')
```

### Emperor

Use the R script makeEmperor.R to output PCA and meta data files. Headers for both files need to be manually altered; "pc vector number" as a column header in pca.dat, and a #sampleID to meta.dat. Then use:
```
conda activate emperor
make_emperor.py -i pca.dat -m meta.dat
conda deactivate
```

=======
>>>>>>> 2fa6b21439b8173f68a9a7bd1dcfa752967ab92b
### GO Analysis

### Pathway Analysis
