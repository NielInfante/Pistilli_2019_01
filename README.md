# RNAseq
Analysis of RNA seq data. From the form:

This run was for RNA isolated from new muscle biopsies of human breast cancer patients.  I want to have these included the group and have everything analyzed together as a whole.  We can simply run non-cancer controls versus all cancer, including the previous data you generated.


## Steps:

###  :white_check_mark: Quality Control
In the reads folder, create a fastqc folder. Then run:

```
for f in *.gz; do echo $f; fastq -o fastqc -t 40 $f; done
```

This makes QC reports for each file. To get a single report, cd to the fastqc folder and use:
```
multiqc -o multiqc .
```

Quality looks good. There are two files called "test" that are different. (They are Ryan testing something else, and can safely be ignored.) I will exclude them from further analysis.


### :white_check_mark: Quantification

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




### :large_orange_diamond: <!--:white_check_mark:--> Differential expression

Use the R script doDESeq.R to find differentially expressed genes, and produce some figures

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
design =~ Run + Group
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

##### TP_fixed
Remove two outliers
```
outPrefix <- 'TP_fixed'
PCA_Group <- 'Group'
design =~ Run + Surgeon + Group
contrast <- c('Group','TP','Control')

meta <- metadata %>% dplyr::filter(Group=='TP' | Group=='Control') %>%
	dplyr::filter(!SampleID %in% c('P4','P310'))

```

##### All_fixed
Remove two outliers
```
outPrefix <- 'All_fixed'
PCA_Group <- 'Cancer'
design =~ Run + Surgeon + Cancer
contrast <- c('Cancer','Cancer','Control')
meta <- metadata %>% dplyr::filter(!SampleID %in% c('P4','P310'))
```

##### TP_vs_Out
Both outliers are in group TP. How do they differ from other members of this group?
```
outPrefix <- 'TP_vs_Out'
PCA_Group <- 'Out'
design =~ Run + Surgeon + Out
contrast <- c('Out','Good','Outlier')

meta <- metadata %>% dplyr::filter(Group=='TP')
meta$Out <- 'Good'
meta <- mutate(meta, Out=ifelse(SampleID %in% c('P4','P310'), 'Outlier','Good'))
```


### :white_check_mark: Emperor

Use the R script makeEmperor.R to output PCA and meta data files. Headers for both files need to be manually altered; "pc vector number" as a column header in pca.dat, and a #sampleID to meta.dat. Then use:
```
conda activate emperor
make_emperor.py -i pca.dat -m meta.dat
conda deactivate
```


### :x: <!--:large_orange_diamond: :white_check_mark:--> GO Analysis

### :x: <!--:large_orange_diamond: :white_check_mark:--> Pathway Analysis
