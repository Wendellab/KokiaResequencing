### From rawreads to final VCF

#### Reads trimming
```
module load  trimmomatic/0.39-zwxnnrx
trimmomatic PE -threads $thr $file1 $file2 $tDir/$name.R1.fq.gz $tDir/$name.U1.fq.gz $tDir/$name.R2.fq.gz $tDir/$name.U2.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:75
```

#### All sample reads mapping to one genome Kk (Kokau.chrONLY.fasta)
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
##### Generating three VCF files: Kc.redo.self.vcf / Kd.yesHAVO.redo.self.vcf / Kk.redo.self.vcf
```
ml sentieon-genomics

set1=$(ls gVCF/Kc*vsKc.gVCF)
set2=$(ls gVCF/Kd*vsKd.gVCF)
set3=$(ls gVCF/Kd*vsKd.gVCF  | sed '/HV/d')
set4=$(ls gVCF/Kk*vsKk.gVCF)
set5=$(ls gVCF/Kc*vsKk.gVCF)
set6=$(ls gVCF/Kd*vsKk.gVCF)
set7=$(ls gVCF/Kd*vsKk.gVCF  | sed '/HV/d')

genome=(Kocoo.chrONLY.fasta Kodry.chrONLY.fasta Kodry.chrONLY.fasta Kokau.chrONLY.fasta Kokau.chrONLY.fasta Kokau.chrONLY.fasta Kokau.chrONLY.fasta)
sets=(set1 set2 set3 set4 set5 set6 set7)
analysis=(Kc.redo.self.vcf Kd.yesHAVO.redo.self.vcf Kd.noHAVO.redo.self.vcf Kk.redo.self.vcf Kc.redo.Kk.vcf Kd.yesHAVO.redo.Kk.vcf Kd.noHAVO.redo.Kk.vcf)

currentgenome="${genome[$SLURM_ARRAY_TASK_ID]}"
currentset="${sets[$SLURM_ARRAY_TASK_ID]}"
outfile="${analysis[$SLURM_ARRAY_TASK_ID]}"

infiles=$(eval "echo \${$currentset[@]}")

echo "my current genome is $currentgenome. I will write to $outfile. My infiles are $infiles"

sentieon driver -t 100 -r $currentgenome --algo GVCFtyper --emit_mode all $outfile $infiles
```

#### Filtering output VCF using depth and maximum alleles
```
#for all three VCFs below, each had been filtered via the pipeline below to produce the combined VCFs

input=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kc.redo.Kk.vcf
input=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kd.yesHAVO.redo.Kk.vcf
input=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kk.redo.self.vcf

###########
thr=30 #NUMBER_THREADS

outputdir=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/00_Kkref_n163
output=$(basename $input .vcf)

ml vcftools bcftools

vcftools --vcf $input --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out $outputdir/$output.variant
vcftools --vcf $input --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out $outputdir/$output.invariant

cd $outputdir

echo "invariant site =" $(wc -l $output.invariant*)
echo "variant site =" $(wc -l $output.variant*)

module load parallel/20220522-sxcww47

parallel bgzip {} ::: $output.*variant.recode.vcf
parallel tabix {} ::: $output.*variant.recode.vcf.gz

bcftools concat --allow-overlaps --threads $thr $output.variant.recode.vcf.gz $output.invariant.recode.vcf.gz -Oz -o $output.combined.vcf.gz

tabix $output.combined.vcf.gz
```

#### Merging all species from three VCFs
```
ml vcftools bcftools

bcftools merge --threads 10 \
/ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/00_Kkref_n163/Kc.redo.Kk.combined.vcf.gz \
/ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/00_Kkref_n163/Kd.yesHAVO.redo.Kk.combined.vcf.gz \
/ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/00_Kkref_n163/Kk.redo.self.combined.vcf.gz \
-Oz -o KcKdKk_Kkn163.combined.vcf.gz

tabix KcKdKk_Kkn163.combined.vcf.gz

bcftools view KcKdKk_Kkn163.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING=0' -q 0.001:minor -Oz -o KcKdKk_Kkn163.combined.bi.vcf.gz

echo "total merged sites for n163" $(zgrep -cv '#' KcKdKk_Kkn163.combined.vcf.gz)
echo "total merged biSNP sites for n163" $(zgrep -cv '#' KcKdKk_Kkn163.combined.bi.vcf.gz)
```


#### Building final VCF with biallelics and adding annotation
```
## Get a bed file that contains the genic regions
zcat /ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/04_gff/Kk.ncbi.gtf.gz |  awk '$3 == "gene" {print $1, $4-1, $5, $10}' OFS='\t' | grep -v "Unplaced" | grep -v "scaffold" | sed 's/"//g' | sed 's/;//g' > /ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/04_gff/Kk.ncbi.gtf.genic.bed

## Get the genic region SNPs
vcftools --gzvcf KcKdKk_Kkn163.combined.bi.vcf.gz --bed /ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/04_gff/Kk.ncbi.gtf.genic.bed --recode --out KcKdKk_Kkn163.combined.bi.genic

## This is for PLINK and LEA, and remove all fixed hetezygosity sites
bcftools view --exclude "F_PASS(GT='het')=1" KcKdKk_Kkn163.combined.bi.genic.recode.vcf -o KcKdKk_Kkn163.combined.bi.genic.nofixhet.vcf

## This is for Pixy to extract variant and invariant sites in genic regions only
vcftools --gzvcf KcKdKk_Kkn163.combined.vcf.gz --bed /ptmp/LAS/jfw-lab/corrinne/redoKokia/Weixuan/04_gff/Kk.ncbi.gtf.genic.bed --recode --out Pixy_genic_n160/KcKdKk_Kkn163.combined.genic
```


