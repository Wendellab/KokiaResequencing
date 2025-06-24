### Building intraspecific VCFs for Kc, Kd, and Kk
#
#### Using trimmed reads to build gVCFs for each species by applying different reference genomes
##### Reference genomes
```
ref=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kocoo.chrONLY.fasta
ref=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kodry.chrONLY.fasta
ref=/ptmp/LAS/jfw-lab/corrinne/redoKokia/Kokau.chrONLY.fasta
```

##### Kc, Kd, and Kk reads mapped to its own reference, HAVO samples mapped to Kc ref
```
# specify the directory, num threads to use, and genome
DIR=/ptmp/LAS/jfw-lab/corrinne/redoKokia/1-filt/
thr=$SLURM_CPUS_PER_TASK

# files to operate over; file2 is based on file1
file1=$(ls -1 $DIR/Kc*.R1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
baseName=$(basename $file1 .R1.fq.gz)
genome="Kc"
name=${baseName}.vs$genome

mkdir $baseName && cd $baseName

echo "the first file is " $file1
echo "the second file is " $file2
echo "the reference file is " $ref

ml sentieon-genomics

export bwt_max_mem=275G
sentieon bwa mem -M -K 10000000 -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" -t $thr $ref $DIR/$baseName.R1.fq.gz $DIR/$baseName.R2.fq.gz | sentieon util sort -o $name.sort.bam -t $thr --sam2bam -i -
sentieon driver -t $thr -i $name.sort.bam --algo LocusCollector --fun score_info $name.score 
sentieon driver -t $thr -i $name.sort.bam --algo Dedup --rmdup --score_info $name.score --metrics $name.dedup.metric $name.dedup.bam 
sentieon driver -t $thr -r $ref -i $name.dedup.bam --algo Realigner $name.realign.bam
sentieon driver -t $thr -r $ref -i $name.realign.bam --algo QualCal $name.recal_data.table

sentieon driver -t $thr -r $ref -i $name.realign.bam -q $name.recal_data.table --algo Haplotyper $name.gVCF --emit_mode gvcf
```

##### Building VCFs for each species: Kc.yesHAVO.redo.self.vcf // Kd.noHAVO.redo.self.vcf // Kk.redo.self.vcf 
```
sentieon driver -t $thr -r $ref --algo GVCFtyper --emit_mode all $analysis $infiles
```
#

#### Filtering VCFs
##### For each VCF, filtering by biallelic and depth for variant and invariant sites
```
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

##### Reference genomes
```

```

##### Reference genomes
```
```
