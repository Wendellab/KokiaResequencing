### From rawreads to final VCF

#### Reads trimming
```
module load  trimmomatic/0.39-zwxnnrx
trimmomatic PE -threads $thr $file1 $file2 $tDir/$name.R1.fq.gz $tDir/$name.U1.fq.gz $tDir/$name.R2.fq.gz $tDir/$name.U2.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:75
```

#### Reads mapping to one genome Kk (Kokau.chrONLY.fasta)
```
# specify the directory, num threads to use, and genome
DIR=/ptmp/LAS/jfw-lab/corrinne/redoKokia/1-filt/
thr=$SLURM_CPUS_PER_TASK
ref=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kokau.chrONLY.fasta

ml sentieon-genomics

export bwt_max_mem=275G
sentieon bwa mem -M -K 10000000 -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" -t $thr $ref $DIR/$baseName.R1.fq.gz $DIR/$baseName.R2.fq.gz | sentieon util sort -o $name.sort.bam -t $thr --sam2bam -i -

#remove duplicated reads
sentieon driver -t $thr -i $name.sort.bam --algo LocusCollector --fun score_info $name.score 
sentieon driver -t $thr -i $name.sort.bam --algo Dedup --rmdup --score_info $name.score --metrics $name.dedup.metric $name.dedup.bam

#realign with sentieon-genomics:
sentieon driver -t $thr -r $ref -i $name.dedup.bam --algo Realigner $name.realign.bam
sentieon driver -t $thr -r $ref -i $name.realign.bam --algo QualCal $name.recal_data.table

#calling gvcf
sentieon driver -t $thr -r $ref -i $name.realign.bam -q $name.recal_data.table --algo Haplotyper $name.gVCF --emit_mode gvcf 
```





#### Calling VCF for each population/species
```
cat AD1_MK_n25.txt | sentieon driver --interval Ah_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Ah_$seq.vcf -
mv $TMPDIR/$output1.Ah_$seq.vcf* $Dir/$output1/

cat AD1_MK_n25.txt | sentieon driver --interval Dh_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Dh_$seq.vcf -
mv $TMPDIR/$output1.Dh_$seq.vcf* $Dir/$output1/
```

#### Filtering output VCF using depth and maximum alleles
```
ml vcftools bcftools

vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Ah_$seq.variant
vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Ah_$seq.invariant

vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Dh_$seq.variant
vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Dh_$seq.invariant

module load parallel/20220522-sxcww47

parallel bgzip {} ::: $output1.*h_$seq.*variant.recode.vcf
parallel tabix {} ::: $output1.*h_$seq.*variant.recode.vcf.gz

bcftools concat --allow-overlaps --threads $thr $output1.Ah_$seq.variant.recode.vcf.gz $output1.Ah_$seq.invariant.recode.vcf.gz -Oz -o $output1.Ah_$seq.combined.vcf.gz
bcftools concat --allow-overlaps --threads $thr $output1.Dh_$seq.variant.recode.vcf.gz $output1.Dh_$seq.invariant.recode.vcf.gz -Oz -o $output1.Dh_$seq.combined.vcf.gz

parallel tabix {} ::: $output1.*h_$seq.*combined.vcf.gz
```

#### Merging all populations/species 
```
ml vcftools bcftools

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD5/AD5.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Ah_$seq.combined.vcf.gz \
-Oz -o MKAD2AD5AD4_n51.Ah_$seq.combined.vcf.gz

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD2/AD2.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD4/AD4.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD5/AD5.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/00_interval/AD1_MK_n25/AD1_MK_n25.Dh_$seq.combined.vcf.gz \
-Oz -o MKAD2AD5AD4_n51.Dh_$seq.combined.vcf.gz

bcftools view MKAD2AD5AD4_n51.Ah_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o MKAD2AD5AD4_n51.Ah_$seq.combined.bi.vcf.gz

bcftools view MKAD2AD5AD4_n51.Dh_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o MKAD2AD5AD4_n51.Dh_$seq.combined.bi.vcf.gz

module load parallel/20220522-sxcww47
parallel tabix {} ::: MKAD2AD5AD4_n51.*h_$seq.combined.bi.vcf.gz
```


#### Building final VCF with biallelics and adding annotation
```
module load picard/2.27.4

picard GatherVcfs \
$(for vcf in *.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename.MKAD2AD4AD5_n51.txt MKAD2AD5AD4_n51.AhDh.combined.bi.vcf.gz -o MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.vcf >  MKAD2AD5AD4_n51.AhDh.combined.bi.rehead.id.vcf
```


