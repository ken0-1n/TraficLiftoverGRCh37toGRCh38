# TraficLiftoverGRCh37toGRCh38
Lift over the Trafic result file from GRch37 to GRCh38.

## preparations
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

gunzip hg19ToHg38.over.chain.gz

chmod 750 liftOver

export PATH=$PATH:${path to liftOVer directory}
```

## run
```
python lift_over_trafic.py TraFiC.vcf TraFiC_hg38.vcf hg19ToHg38.over.chain 
```
