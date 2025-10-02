# Task
Multiplexing is done to reduce sequencing cost by combining multiple samples and sequencing them together all at once.
If the samples are from different individuals, we can utilize genetic differences between them to assign cells to the individuals.
Genotype information is usefull but not mandatory. 
Thankfully, (at the time of this writing), genotype information of 52 strains of mice are avalible from the mouse genome project (https://www.mousegenomes.org/snps-indels/). We can download the VCF file that lists all the SNPs (differences in genetic sequence) between these mice.  
If you haven't worked with VCF files before, I suggest take some time to understand what they are. They are in fact rather straightforward. 

# What you need
1. The bam file for your scRNA data. Can be found in the 10X output folder.
2. The VCF file. Contains the genotypic differences. (Check out the Mouse Genome Project)
3. Barcodes file. A list of barcodes that you want to demultiplex. The bam file will have data for all the barcodes, incliding empty droplets. You can provide a custom list of barcodes or use filtered_feature_bc_matrix/barcodes.tsv.gz in the 10x cellranger output.
4. Individuals file, just the names of strains that are supposed to be in the data. If you omit a starin, even if it is in the VCF file, it will be ignored while demultiplexing.
5. Libraries  
  a. bcftools, I used a singularity image from dockerhub  
  b. samtools, I used a singularity image from dockerhub  
  c. Demuxalot (or a similar tool). I installed demuxalot from pip. It also is present in Demuxify's singularity image the authors shared. I couldn't get it to run that way.  

# Steps
## 1. Filter the VCF
The VCF file has 52 strains, where we need only 3 of those. So we filter for what we need, using bcftools:

```
bcftools view -f PASS \
              --samples C57BL_6NJ,NZO_HlLtJ,CAST_EiJ \
              --trim-alt-alleles \
              --min-ac 1 \
              --output B6_NZO_CAST.vcf.gz \
              --output-type z \
              mgp_REL2021_snps.vcf.gz
```
This command does the following:  
The view subcommand of the bcftools program filters for all records (the genetic differences) that **PASS** the quality checks,  
grabs 3 "samples" namely **C57BL_6NJ, NZO_HlLtJ and CAST_EiJ** which are the strains we are interested in,  
**trim alternate alleles** NOT SEEN in our subset. **--min-ac 1** might be redundant here, it is supposed to remove everything that has no allele.  

You don't really have to remove the unnecessary strains from the VCF file. You will provide the "individuals" file to specifiy the strains of interest anyway. However shrinking the VCF file speeds things up.

We also need tp create an index file fo the new VCF. This provides faster access. And also, some tools demand you have the index file. We will keep doing this as we generate new vcf files.  
```
bcftools index B6_NZO_CAST.vcf.gz
```

## 2. Append "chr-" to chromosome names in the VCF
The demultiplexing tool requires the chromosome (called contigs in the VCF file) names to be the same. 
In the 10x generated bam file the chromosome names have "chr-" prefixes in them (like chr1, chr2...). Whereas the VCF file does not (just 1,2,3 etc). 
We fix those by tweaking the vcf file. You could also tweak the bam file instead. 

```
bcftools annotate \
  --rename-chrs chrom_map.txt \
  -Oz -o B6_NZO_CAST.chr.vcf.gz \
  B6_NZO_CAST.vcf.gz

bcftools index B6_NZO_CAST.chr.vcf.gz
```
chrom_map.txt is just a mapping between the old an new names. Like:   
```
1   chr1  
2   chr2  
3   chr3  
...
```
I included it in the repo.

## 3. Sort the chromosomes (contigs)
The tool demultiplexing too also requires the chromosomes to be sorted in the same order in the bam file and the VCF file.  
Unfortunately the 10X generated bam files have lexicographical ordering (1,10,11,...,2,3...) and the mouse genome project uses numerical (1,2,3...).

This is how we sort the records (not the header!):

```
bcftools view \
  --regions-file bam_regions.bed \
  -Oz -o B6_NZO_CAST.reordered.vcf.gz \
  B6_NZO_CAST.chr.vcf.gz

bcftools index B6_NZO_CAST.reordered.vcf.gz
```

Here, bam_regions.bed is in fact a tsv file with the ordering we mentioned. Below was the command I used to generate it, but you can manualy create it if you like. I also included it in the repo:
```
samtools idxstats \
    possorted_genome_bam.bam \
    | awk '$1 !~ /^\*/ {print $1 "\t1\t" $2}' > bam_regions.bed
```

So far, we fixed the chr- prefix issue and sorted the records, BUT did not fix the VCF header! I did that semi-manually. First extarct the current header into a txt file:

```
bcftools view -h B6_NZO_CAST.reordered.vcf.gz > header.txt
```

Then **manually edited** the txt file, sorted the contigs. Just check the file, you will understand. 
Now we can put the sorted header back in there:

```
bcftools reheader \
  --header header.txt \
  -o B6_NZO_CAST.fully.reordered.vcf.gz \
  B6_NZO_CAST.reordered.vcf.gz 

bcftool index B6_NZO_CAST.fully.reordered.vcf.gz
```

## 4. Demultiplex!
Finally we can run the tool itself with the VCF file that we prepared along with the bam file and others. I added the scripts to run demuxalot on an hpc cluster. Demuxalot is a python library, you can also use it as such.

