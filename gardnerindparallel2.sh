#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load gnu_parallel/201612222
module load trimmomatic/0.35
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12
module load samtools/1.2

## Recommended submission:
        # 12-18 cores
        # 2 hours
        # 16gb

#makes a sub-directory to run untrimmed reads to avoid confusion with
#previous data:
mkdir /scratch/aubcls03/group/ind

### Change directory to the sub-directory
# cd /scratch/YOUR_ID/group
cd /scratch/aubcls03/group/ind

#####   copy all .fastq.gz to  scratch directory in parallel using GNU parallel
#  ls = make a (memory) list of the .fastq.gz files in that directory
#  | in parallel using as many cores as possible (one job on each core) but no more jobs, cp the file in the list
# options: -- eta give the estimate time --dry-run to see if being parsed correctly
# You can do this seperately for the files you need to run using this code from before changing the specific file names

#This copies the contatenated, untrimmed lizard data, as well as the already indexed Anolis genome to the ind sub-directory within the group directory in scratch, in parallel
ls /scratch/aubcls03/group/*All_R1.fastq.gz | time parallel -j+0 --eta 'cp {} .'
ls /scratch/aubcls03/group/Anolis_carolinensis* | time parallel -j+0 --eta 'cp {} .'

#------------------------------------------------------------------------------------------------------------------------------
#Now for mapping (The genome reference has already been indexed (for how to index, see gardnergroupparallel2.sh))




#Setting up for loop to run the mapping: 

ls *All_R1.fastq.gz |cut -d "_" -f 1,2 | sort | uniq > list6
numfiles=`less list6 | wc -l`
files=($(ls *All_R1.fastq.gz |cut -d "_" -f 1,2 | sort | uniq))

#you can change the i=0 to whatever file this stops on so it will pick back up running where it leaves off!
for ((i=10; i<=$numfiles-1; i++))
do 
echo Now processing : ${files[i]} ...

bwa mem -t 20  Anolis_carolinensis_genome ${files[i]}_All_R1.fastq.gz > ${files[i]}.sam
bwa mem -t 20  Anolis_carolinensis_genome ${files[i]}_CTTGTA_All_R1.fastq.gz > ${files[i]}_in.sam


done 

# Output =  2_TGACCA.sam
# Output =  10-In_CTTGTA.sam
# Output =  1_In_in.sam
#will give several blank files as a result of not being able to find the ones matching the _in files
#again, remove any empty files

find . -size 0 -delete


###############  Using samtools to process the bam Input: HS06_GATCAG_All.sam
#  Make the list of prefixes for all the .sam files we want to process with Samtools
ls | grep "trimmed.fastq.gz" |cut -d "_" -f 1,2 | sort | uniq  > list7

# Use a loop to process through the names in the list using samtools
while read i;
do
   	## convert .sam to .bam and sort the alignments
        # -@ is the number of threads
samtools view -@ 12 -bS ${i}.sam  | samtools sort -@ 12 -  ${i}_sorted  
 # Example Input: 2_TGACCA.sam; Output: 2_TGACCA_sorted.bam
        ## index the sorted .bam
samtools index  ${i}_sorted.bam
        ## Tally counts of reads mapped to each transcript; and calcuate the stats.
samtools idxstats   ${i}_sorted.bam     >	${i}_Counts.txt
samtools flagstat	${i}_sorted.bam         >	${i}_Stats.txt
done<list7

#again, remove any empty files

find . -size 0 -delete


#have to re-run for the 1_in files

while read i;
do
        ## convert .sam to .bam and sort the alignments
        # -@ is the number of threads
samtools view -@ 12 -bS ${i}_in.sam  | samtools sort -@ 12 -  ${i}_in_sorted
 # Example Input: 2_TGACCA.sam; Output: 2_TGACCA_sorted.bam
        ## index the sorted .bam
samtools index  ${i}_in_sorted.bam
        ## Tally counts of reads mapped to each transcript; and calcuate the stats.
samtools idxstats   ${i}_in_sorted.bam     >       ${i}_in_Counts.txt
samtools flagstat       ${i}_in_sorted.bam         >       ${i}_in_Stats.txt
done<list7

#again, remove any empty files

find . -size 0 -delete

#have to setup another loop to delete count and stat files that have the _in we dont need:
#This makes a list of all the files we want to delete so we are left with the ones that actually matter to us:

ls -s *_in_Counts.txt | sort | head -n 11 | cut -d " " -f 4  > list9
ls -s *_in_Stats.txt | sort | head -n 11 | cut -d " " -f 3 > list10
cat list9 list10 > list11

#The loop that deletes the unnecessary files:
while read z;
do
rm ${z}

done<list11

#removes the ones that we dont need:
rm 1_In_Counts.txt
rm 1_In_Stats.txt

##### Make a directory for my results in my home folder
# mkdir /home/YOUR_ID/class_shared/YOUR_NAME/fastqc
#  *** EDIT ***
mkdir /home/aubcls03/class_shared/Steven/group/ind
mkdir /home/aubcls03/class_shared/Steven/group/ind/BWA_Counts
mkdir /home/aubcls03/class_shared/Steven/group/ind/BWA_Stats

### Move the results output to my directory for safe keeping
## /home/YOUR_ID/class_shared/YOUR_NAME/group/BWA_Counts
#  *** EDIT ***
cp *Counts.txt /home/aubcls03/class_shared/Steven/group/ind/BWA_Counts
cp *Stats.txt /home/aubcls03/class_shared/Steven/group/ind/BWA_Stats

# Move to your home directory and tarball the folders so they are ready to transfer back to your computer.
#  *** EDIT ***
cd /home/aubcls03/class_shared/Steven/group/ind
tar -cvzf BWA_Counts.tar.gz /home/aubcls03/class_shared/Steven/group/ind/BWA_Counts
tar -cvzf BWA_Stats.tar.gz /home/aubcls03/class_shared/Steven/group/ind/BWA_Stats

#----------------------------------------------------------------------------------------
#runs final clean up for files we dont need from the samtools portion of the code:

#move back to your scratch directory
#cd /scratch/YOURID/group/
cd /scratch/aubcls03/group/ind

ls -s *_sorted* | sort | head -n 24 | cut -d " " -f 8  > list12



#The loop that deletes the unnecessary files:
while read z;
do
rm ${z}

done<list12

#now we are left with only the files from the trimming,mapping,samtools steps we need!

#-------------------------------------------------------------------------
#now for the stats part...

