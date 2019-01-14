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

Quality looks good. There are two files called "test" that are different, I will exclude them from further analysis.


### Quantification

Gat the reference from Ensembl, release 94. Got cDNA_all and ncRNA, cat them together and build salmon index.

```
conda activate

# Build index
salmon index -t CDNA_and_NCrna.fa -i hg38.94 --type quasi -k 31

# Do quantification
for f in P*R1*; do echo $f; salmon quant -i ../Data/hg38.94 --libType A --gcBias -p 20 --numBootstraps 50 -o ../salmon/${f%_S*} -1 $f -2 ${f/R1/R2}; done

conda deativate
```


In the reads folder, use salmon to quantify each read

### Differential expression

Use the R script doDESeq.R to find differentially expressed genes, and produce some figures

### GO Analysis

### Pathway Analysis
