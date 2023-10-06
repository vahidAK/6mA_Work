# 6mA_Work
Detection of CpG methylation and chromatin accessibility simultaneously using nanopore sequencing.
To run the commands here you need these tools:  
[Megalodon](https://github.com/nanoporetech/megalodon) and ONT Guppy, [minimap2](https://github.com/lh3/minimap2), [Clair3](https://github.com/HKU-BAL/Clair3), [WhatsHap](https://github.com/whatshap/whatshap), [Samtools](https://github.com/samtools/samtools), [bgzip](http://www.htslib.org/doc/bgzip.html), [Tabix](http://www.htslib.org/doc/tabix.html), [Nanopolish](https://github.com/jts/nanopolish), [NanoMethPhase](https://github.com/vahidAK/NanoMethPhase), [DSS R package](https://bioconductor.org/packages/release/bioc/html/DSS.html)  

Table of Contents
=================
* **[6mA calling](https://github.com/vahidAK/6mA_Work/blob/master/README.md#6ma-calling)**
  * [Base and 6mA calling using megalodon](https://github.com/vahidAK/6mA_Work/blob/master/README.md#base-and-6ma-calling-using-megalodon)
  * [Detection of DNA accessibility peaks](https://github.com/vahidAK/6mA_Work/blob/master/README.md#detection-of-dna-accessibility-peaks)
  * [Phasing 6mA calls](https://github.com/vahidAK/6mA_Work/blob/master/README.md#phasing-6ma-calls)
* **[CpG methylation](https://github.com/vahidAK/6mA_Work/blob/master/README.md#cpg-methylation)**
  * [Methylation calling](https://github.com/vahidAK/6mA_Work/blob/master/README.md#methylation-calling)
  * [Phasing CpG methylation calls](https://github.com/vahidAK/6mA_Work/blob/master/README.md#phasing-cpg-methylation-calls) 

# 6mA calling
## Base and 6mA calling using megalodon 
To make the base and 6mA calling faster you can run this command on each fast5 file separately or on a bunch of fast5 files. Related models are provided in the [Models](https://github.com/vahidAK/6mA_Work/blob/main/Models) directory.
```
megalodon path/to/directory/containing/fast5/ \ 
  --guppy-server-path /path/to/guppy_basecall_server \
  --guppy-config dna_r9.4.1_450bps_sup_prom.cfg \
  --reference /path/to/reference.fa # Alternatively, use a minimap index of the reference (.mmi) for more efficient memory usage \
  --output-directory /path/to/output/directory \
  --guppy-params "-m /path/to/BaseCalingModel/from/this/repository" \ # Models/Guppy_BaseCallingModel_r9.4.1.json
  --device <device name, for example "cuda:0"> \
  --remora-model /path/to/RemoraModel/from/this/repository" \ # Models/Remora_6mACallingModel_r9.4.1.onnx
  --outputs basecalls per_read_mods mods \
  --mod-binary-threshold 0.75 \
  --guppy-timeout 300 \
  --processes <# of processes>
```
Output per-site results will be in the output directory. If you run the command for each fast5 file or a bunch of them separately, you can put all the per-site results from all the runs in a directory and use the [aggregate_sites.py](https://github.com/vahidAK/6mA_Work/blob/main/aggregate_sites.py) script to merge per-site results for each chromosome.
## Detection of DNA accessibility peaks
By doing a differential 6mA analysis of treated vs untreated sample using DSS R package we can detect 6mA DMRs/peaks. [HG002](https://zenodo.org/record/8408392) untreated 6mA frequency data can be used as the untreated control sample for data based on hg38. Using this command in R you can perform differential 6mA analysis:
```
library("DSS")
options(scipen = 15)
file1= read.table("/path/to/6mA_frequency_file/from/Treated_sample",header=TRUE,sep = '\t')
file1= file1[,c(1,2,5,6)] # Selecting for columns that represent chromosome, adenine position, coverage, number of reads showing 6mA 
colnames(file1)= c("chr","pos","N","X")

file2= read.table("/path/to/6mA_frequency_file/from/Untreated_sample",header=TRUE,sep = '\t')
file2= file2[,c(1,2,5,6)] # Selecting for columns that represent chromosome, adenine position, coverage, number of reads showing 6mA 
colnames(file2)= c("chr","pos","N","X")
DSObject<- makeBSseqData(list(file1,file2),c("C1", "N1"))

test<- DMLtest(DSObject, 
               group1=c("C1"), 
               group2=c("N1"), 
               equal.disp = FALSE,
               smoothing = TRUE,
               smoothing.span = 150,
               ncores=<# of cores>)

DM_region<- callDMR(test, delta=0.1, p.threshold=0.001, minlen=150, minCG=30, dis.merge=150, pct.sig=0.5)
write.table(DM_region,"/path/to/output_file", sep="\t", row.names=F, quote=F)
```
Afterward, you need to keep regions that have higher 6mA in the treated sample (diff.Methy > 0). You can also select more robust regions, for example, keeping regions with diff.Methy > 0.15.

## Phasing 6mA calls 
During the previous steps, we got unphased 6mA and peak data and here we are going to detect haplotype-specific 6mA calls and peaks
### Mapping, variant calling, and phasing 
If there are multiple fastq files, you need to concatenate all of them into a single file. Mapping and sorting can be done as follow:
```
minimap2 -ax map-ont --MD -L -t <# of threads> /path/to/reference.fa /path/to/input.fastq > file.sam
samtools sort -@ <# of threads> file.sam -o file.bam && rm file.sam
samtools index -@ <# of threads> file.bam
```
Now we can call variants from the bam using Clair3:
```
run_clair3.sh --bam_fn=file.bam \
 --ref_fn=/path/to/reference.fa \
 --output=/path/to/output_directory \
 --threads=<# of threads> \
 --platform=ont \
 --model_path=/path/to/clair3/model_directory
```
The result variants will be in the output directory. Then we should select the quality passed variants, bgzip and index the file:
```
gunzip -c /path/to/output_directory/merge_output.vcf.gz | awk -F'\t' '$1~/^#/ || $7=="PASS"' > clair3_passed.vcf
bgzip clair3_passed.vcf
tabix -p vcf clair3_passed.vcf.gz
```
Now variants can be phased using whatshap and then reads can also be tagged to haplotypes:
```
whatshap phase --ignore-read-groups \
 --reference /path/to/reference.fa \
 -o clair3_passed_whatshap_phased.vcf \
 clair3_passed.vcf.gz \
 file.bam

bgzip clair3_passed_whatshap_phased.vcf
tabix -p vcf clair3_passed_whatshap_phased.vcf.gz

# now haplotagging reads
whatshap haplotag --ignore-read-groups \
 --skip-missing-contigs \
 --reference /path/to/reference.fa \
 -o file_haplotagged.bam \
 --output-haplotag-list file_haplotagged_Reads.tsv \
 clair3_passed_whatshap_phased.vcf.gz \
 file.bam
```
The list of haplotype 1 and haplotype 2 reads will be in the file_haplotagged_Reads.tsv file. You need to put haplotype 1 and 2 read IDs in two separate files for the following step (1 read ID per line).
### Phasing 6mA calls and detecting allele-specific differentially accessible regions (aDARs)
Using the following command you can get the per-site 6mA data for each haplotype
```
megalodon_extras aggregate run \
 --outputs mods \
 --megalodon-directory /path/to/directory/containes/megalodon_results \
 --read-ids-filename Haplotype1_ReadIDs.tsv \
 --output-suffix out_put_suffix_haplotype_1 \
 --mod-binary-threshold 0.75 \
 --processes <# of processes>

megalodon_extras aggregate run \
 --outputs mods \
 --megalodon-directory /path/to/directory/containes/megalodon_results \
 --read-ids-filename Haplotype2_ReadIDs.tsv \
 --output-suffix out_put_suffix_haplotype_2 \
 --mod-binary-threshold 0.75 \
 --processes <# of processes>
```
Output per-site results will be in the megalodon results output directory. If you run the command for each fast5 file or a bunch of them separately, you can put all the per-site results from all the runs for each haplotype in a directory and use the [aggregate_sites.py](https://github.com/vahidAK/6mA_Work/blob/main/aggregate_sites.py) script to merge per-site results for each haplotype for each chromosome.  

After phasing 6mA, allele-specific accessibility peaks can be detected using DSS by performing differential 6mA analysis between haplotypes.  
```
library("DSS")

file1= read.table("/path/to/6mA_frequency_file/from/haplotype_1",header=TRUE,sep = '\t')
file1= file1[,c(1,2,5,6)] # Selecting for columns that represent chromosome, adenine position, coverage, number of reads showing 6mA 
colnames(file1)= c("chr","pos","N","X")

file2= read.table("/path/to/6mA_frequency_file/from/haplotype_2",header=TRUE,sep = '\t')
file2= file2[,c(1,2,5,6)] # Selecting for columns that represent chromosome, adenine position, coverage, number of reads showing 6mA 
colnames(file2)= c("chr","pos","N","X")

DSObject<- makeBSseqData(list(file1,file2),c("C1", "N1"))
test<- DMLtest(DSObject, 
               group1=c("C1"), 
               group2=c("N1"), 
               equal.disp = FALSE,
               smoothing = TRUE,
               smoothing.span = 150,
               ncores=<# of cores>)

DM_region<- callDMR(test, delta=0.15, p.threshold=0.001, minlen=150, minCG=30, dis.merge=150, pct.sig=0.5)
write.table(DM_region,"/path/to/output_file", sep="\t", row.names=F, quote=F)
```
Afterward, You can also select more robust regions. For example, keeping regions with diff.Methy > 0.2 or diff.Methy < -0.2 between haplotypes.  

# CpG methylation
## Methylation calling
Here we use nanopolish (You may use [f5c](https://github.com/hasindu2008/f5c) >=v0.7 which is an optimised re-implementation of nanopolish or other methylation callers). First, we need to index fastq file and then call per-read methylation.

```
nanopolish index -d /path/to/fast5s_directory input.fastq

nanopolish call-methylation \
  -t <number_of_threads> \
  -q cpg \
  -r input.fastq \
  -b file.bam \
  -g /path/to/reference.fa > CpG_MethylationCall.tsv
# For faster calling you can also run nanopolish call-methylation command on each chromosome or genomic intervals separately.
```
Now per-site unphased CpG methylation can be detected by running the ```calculate_methylation_frequency.py``` script from nanopolish:
```
calculate_methylation_frequency.py --call-threshold 1.5 -s CpG_MethylationCall.tsv > CpG_MethylationFrequency.tsv
```
## Phasing CpG methylation calls
Here we use NanoMethPhase. We first need to preprocess the methylation call file and then phase methylation data.
```
nanomethphase.py methyl_call_processor -mc CpG_MethylationCall.tsv -t <# of threads> | sort -k1,1 -k2,2n -k3,3n | bgzip > NanoMethPhase_CpG_MethylationCall.tsv.gz
tabix -p bed NanoMethPhase_CpG_MethylationCall.tsv.gz

nanomethphase.py phase -v /path/to/clair3_passed_whatshap_phased.vcf.gz \
  -mc NanoMethPhase_CpG_MethylationCall.tsv.gz \
  -o /path/to/output_directory_and_prefix \
  -of methylcall \
  -b /path/to/file.bam \
  -r /path/to/reference.fa \
  -t <# of threads>
```
These commands will give us haplotype 1 and 2 CpG methylation frequency files. 
### Detecting allele-specific differential methylated regions (aDMRs)
Now we can perform differential 5mC analysis between haplotypes to detect allele-specific DMRs.
```
nanomethphase.py dma -c 1,2,4,5,7 \
  -ca /path/to/5mC-Frequency_haplotype_1.tsv \
  -co /path/to/5mC-Frequency_haplotype_2.tsv \
  -o /path/to/output_directory \
  -op output_prefix
```
 The Output callDMR.txt file will include allelic DMRs. You can also select more robust regions. For example, keeping regions with diff.Methy > 0.2 or diff.Methy < -0.2 between haplotypes. 
