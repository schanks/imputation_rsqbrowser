# Genomic features 
Data sources and scripts used to perform analyses associating Rsq with local genomic features for Figures 4,S9,S10.  Includes scripts for testing association between genomic features/repeats and dichotomous imputation quality (logit.R) and continuous imputation quality (zoib.R)

# Download genomic features
GC content (chr20.gc):
```
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.wig.gz
```

Structural variants: UCSC Genome Browser: Variation -> DGV Struct Var -> dgvM
erged. Keep columns chromStart, chromEnd, varType (chr20.struct)

Segmental duplications: UCSC Genome Browser: Repeats -> Segmental Duplications -> genomicSuperDups. Keep columns chrom, chromStart, chromEnd (chr20.sd)

Repeats: UCSC Genome Browser: Repeats -> RepeatMasker -> rmsk. Keep columns chromosome, repStart, repEnd, repClass (chr20.repeat)

Recombination rate:
Downloaded from http://faculty.washington.edu/browning/beagle/beagle.html (chr20.recomb)

## Input files
Needed input files for each ancestry:
wgs.bi.snv.sites: Tab delimited file with CHR, POS, REF, ALT, AF for each sequenced biallelic SNV.  Can be created from ../../[ancestry].bi.snv.tab downloaded from RsqBrowser.
chr20.gccontent, chr20.struct, chr20.sd,chr20.repeat,chr20.recomb from Step 1

chr*.bi.snv.vcf.gz: WGS VCF of biallelic SNVs

wgs.bi.snv.sites: Tab delimited file with CHR, POS, REF, ALT, AF for each sequenced biallelic SNV

../../[ancestry].bi.snv.tab: Variant-level Rsq file, downloaded from [RsqBrowser](https://imputationserver.sph.umich.edu/rsq-browser/downloads)

## Scripts:
gc.R: Outputs mean GC content in 10kb window centered at each variant

sv.R: Outputs number of structural variants in 10kb window

segdup.R: Outputs presence/absence of segmental duplications in 10kb window

repeat.R: Outputs number of repeats (total) in 10kb window

repeat_class.R: Outputs number of repeats by repeat class in 10kb window

recombMean.R: Outputs mean recombination rate in 10kb window 

distance.R: Outputs distance of each variant to nearest genotyped marker

## Calculate aggregate statistics around 10kb window of each variant for each feature.  Perform for each ancestry.
```
Rscript gc.R -r wgs.bi.snv.sites -g chr20.gccontent -c chr20 -w 5000 -o chr20.gc
Rscript sv.R -r wgs.bi.snv.sites -s chr20.struct -c chr20 -w 5000 -o chr20.sv
Rscript segdup.R -r wgs.bi.snv.sites -s chr20.sd -c chr20 -w 5000 -o chr20.segdup
Rscript repeat.R -r wgs.bi.snv.sites -p chr20.repeat -c chr20 -w 5000 -o chr20.rep
Rscript repeat_class.R -r wgs.bi.snv.sites -p chr20.repeat -c chr20 -w 5000 -o chr20.class
Rscript recombMean.R -r wgs.bi.snv.sites -c chr20.recomb -c chr20 -w 5000 -o chr20.recomb.mean
Rscript distance.R -r wgs.bi.snv.sites -a ../../../arrays/omni25.bed -c chr20 -o chr20.omni25
```

Produces the following output: 

chr20.gc: Mean GC content in 10kb window centered at each variant

chr20.sv: Number of SVs in 10kb window

chr20.seg: Presence of segmental duplication in 10kb window

chr20.rep: Number of repeats (total) in 10kb window

chr20.class: Number of repeats by repeat class in 10kb window

chr20.recomb.mean: Mean recombination rate in 10kb window

chr20.omni25: Distance of each variant to nearest genotyped marker on Omni2.5M array

## Subset to independent SNPs
```
for i in {1..22}; do echo "plink2 --vcf chr$i.bi.snv.vcf.gz --set-all-var-ids @:# --indep-pairwise 50 5 0.2 --out chr$i"; done > command.sh
cat chr1.prune.in | sed 's/^/chr/g' | sed 's/:/\t/g' > all.prune.in
for i in {2..22}; do echo "cat chr$i.prune.in | sed 's/^/chr/g' | sed 's/:/\t/g' >> all.prune.in"; done > command.sh
```
Will need to merge all.prune.in with wgs.bi.snv.sites to create [ancestry].bi.snv.indep.sites

## Perform statistical analyses
a) Logistic regression 
```
Rscript logit.R
```
Produces output in the form [ancesty].logit.Omni25_TOPMed.coefs (for Figure 4) and [ancestry].represults.Omni25_TOPMed.coefs (for Figure S9).  Also calculates Nagelkerke Rsq for Figure 4B but does not output results. Example script is included for METSIM

b) Zero-one inflated beta regression
```
Rscript zoib.R
```
Produces output in form [ancestry].zoibresults.Omni25_TOPMed.coefs.  Example is included for METSIM







