### We are sharing files from the lab directory. Create a symbolic link to data in the /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia directory and put it in your own directory.
#### Saves limited server space

### Example of sympbolic link: running in /local/workdir/RNA-Seq_workshop/eag
ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_C1_S1_R1_001.fastq SY_C1.fastq

### Run in shell script: 
ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_C2_S2_R1_001.fastq SY_C2.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_C3_S3_R1_001.fastq SY_C3.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_C4_S4_R1_001.fastq SY_C4.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_N1_S5_R1_001.fastq SY_N1.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_N2_S6_R1_001.fastq SY_N2.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_N3_S7_R1_001.fastq SY_N3.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_N4_2_S18_R1_001.fastq SY_N4.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_P1_2_S19_R1_001.fastq SY_P1.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_P2_S8_R1_001.fastq SY_P2.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_P3_S9_R1_001.fastq SY_P3.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_P4_2_S20_R1_001.fastq SY_P4.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_PN1_2_S21_R1_001.fastq SY_PN1.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_PN2_S10_R1_001.fastq SY_PN2.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_PN3_S11_R1_001.fastq SY_PN3.fastq

ln -s /local/workdir/RNA-Seq_workshop/Soybean_SCN_Pochonia/SY_PN4_2_S22_R1_001.fastq SY_PN4.fastq


### Run fastqc on file to inspect quality. To trim/remove adapters we will use Trimmamatic. No adapters were detected, but it is worth running the adapter removal step anyway. First, Headcrop will cut a specific amount of basepairs off the front, we will remove 13 bp. 
java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -threads 16 -phred33 SY_C1.fastq trimSY_C1.fastq HEADCROP:13

### Fastqc again, and we see that there are variable read lengths. Now we want to run Trimmomatic again to make sure that there is a minimum length. 
java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -threads 16 -phred33 trimSY_C1.fastq trim2SY_C1.fastq \ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:63 


### Bash script for running trimmomatic on all our files. 
#!/bin/bash

files=$(ls | grep '.fastq$')

for f in $files

do

	java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -phred33 $f crop_$f HEADCROP:13

	java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -phred33 crop_$f trim_crop_$f ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:63

done


### The data is of soybeans grown alone (Control-C), with potential nematode fungal pathogen Pochonia (P), with Soybean cyst nematodes (N), or with Soybean cyst nematodes and Pochonia (PN).
### We will use the soybean genome as the reference genome. We are interested to see how the soybean is responding with the fungi, nematodes, or both. 

### Download GCA_000004515.5_Glycine_max genome from NCBI. 

tar -xvf genome_assemblies_genome_fasta.tar

#Makes the directory: ncbi-genomes-2023-03-30



### We need to index of the reference genome. Remember to use the with updated seq length number of 63 as the overhang. 

export PATH=/programs/STAR-2.7.10b:$PATH

STAR --runMode genomeGenerate --runThreadN 24 --genomeDir indexgenome --genomeFastaFiles /local/workdir/RNA-Seq_workshop/eag/ncbi-genomes-2023-03-30/GCA_000004515.5_Glycine_max_v4.0_genomic.fasta  --sjdbGTFfile /local/workdir/RNA-Seq_workshop/eag/ncbi-genomes-2023-03-30/GCA_000004515.5_Glycine_max_v4.0_genomic.gtf --sjdbOverhang 63 



### Now we map our RNAseq reads to the indexed reference genome. Use shell script, star_align.sh

#!/bin/bash -l 

#SBATCH –nodes=1

#SBATCH –mem=16000

#SBATCH –time=60:00

#SBATCH –partition=regular

#SBATCH –chdir=/local/workdir/RNA-Seq_workshop/eag/exercise1

#SBATCH –job-name=EAG_STARmap

#SBATCH –output=STARtest.out

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_C1.fastq --outFileNamePrefix SY_C1_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_C2.fastq --outFileNamePrefix SY_C2_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_C3.fastq --outFileNamePrefix SY_C3_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_C4.fastq --outFileNamePrefix SY_C4_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_N1.fastq --outFileNamePrefix SY_N1_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_N2.fastq --outFileNamePrefix SY_N2_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_N3.fastq --outFileNamePrefix SY_N3_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_N4.fastq --outFileNamePrefix SY_N4_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_P1.fastq --outFileNamePrefix SY_P1_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_P2.fastq --outFileNamePrefix SY_P2_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_P3.fastq --outFileNamePrefix SY_P3_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_P4.fastq --outFileNamePrefix SY_P4_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_PN1.fastq --outFileNamePrefix SY_PN1_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_PN2.fastq --outFileNamePrefix SY_PN2_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_PN3.fastq --outFileNamePrefix SY_PN3_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 

/programs/STAR-2.7.10b/STAR --quantMode GeneCounts --genomeDir indexgenome --runThreadN 24  --readFilesIn /local/workdir/RNA-Seq_workshop/eag/trim_crop_SY_PN4.fastq --outFileNamePrefix SY_PN4_ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate 


### Look at the log.final.out files to see if mapped reads #'s make sense! 
less SY_C1_Log.final.out    

less SY_C1_ReadsPerGene.out.tab  
  
### SY_C1_Log.final.out has basicstatistics of the alignment & SY_C1_ReadsPerGene.out.tab has read counts per gene. column 1 = gene ID,  column 2= counts for unstranded RNA-seq, column 3= counts for reads aligned with plus strand of RNA, column 4= counts for reads aligned with minus strand of RNA 

### We will now create a file that we will use in our differencial abundance analysis.
### We need the gene name in column 1, and then we will use the unstranded RNA-seq results that are in column 2 of the ReadsPerGene.out.tab files. 
### We have 16 files. I can paste them together, then cut the 2nd column information for each sample: Therefore, I count by 4's from 2 to 62, remove the stats columns and output a new file. 
#### paste: Merge the 16 files by columns; cut -f1,2,6,10,14,18,22,etc: Extract columns 1,2,6,10,14,18,22,etc. from the merged data (column 1 is the gene name, columns 2-62 are second column from each individual file); tail -n +5: Discard the first 4 lines of statistics summary and start from line 5; >gene_count.txt: Write the result into a file gene_count.txt

paste *_ReadsPerGene.out.tab| cut -f1,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62 | tail -n +5 > gene_counts.txt 


### This file is used for the next step in R
