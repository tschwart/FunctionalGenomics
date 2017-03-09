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

#makes directory for our stuff(change your code to your scratch):
mkdir /scratch/aubcls03/group

### Change directory to the scratch directory
# cd /scratch/YOUR_ID/group
cd /scratch/aubcls03/group

#####   copy all .fastq.gz to  scratch directory in parallel using GNU parallel
#  ls = make a (memory) list of the .fastq.gz files in that directory
#  | in parallel using as many cores as possible (one job on each core) but no more jobs, cp the file in the list
# options: -- eta give the estimate time --dry-run to see if being parsed correctly
# You can do this seperately for the files you need to run using this code from before changing the specific file names

#This copies our lizard data to the group directory in our scratch, in parallel
        #  cp /home/aubcls03/class_shared/Exp2_Lizard_EpigeneticDataset/*.fastq.gz .
ls /home/aubcls03/class_shared/Exp2_Lizard_EpigeneticDataset/*.fastq.gz | time parallel -j+0 --eta 'cp {} .'

#### Create list of names:
# ls (list) contents of directory with fastq files, cut the names of the files at
        #underscore characters and keep the first three chunks (i.e. fields; -f 1,2,3),
        #sort names and keep only the unique ones (there will be duplicates of all
        #file base names because of PE reads), then send the last 6 lines to a file
        #called list with tail
# HS03_TTAGGC_L005_R1_001.fastq.gz
                # 1 = HS03
                # 2 = TTAGGC
                # 3 = L005
                # 4 = R1
                # 5 = 001

ls | grep ".fastq.gz" |cut -d "_" -f 1,2 | sort | uniq > list

#### Make the list then use that list to cat the files together. We have singld reads, but I ran it twice because 1 set of
#the files doesnt work using the first two fields with the _L00* code
ls | grep ".fastq.gz" |cut -d "_" -f 1,2 | sort | uniq | parallel cat {1}_L00*_R1_*.fastq.gz '>' {1}_All_R1.fastq.gz ::: ${i}
ls | grep ".fastq.gz" |cut -d "_" -f 1,2 | sort | uniq | parallel cat {1}_CTTGTA_L00*_R1_*.fastq.gz '>' {1}_CTTGTA_All_R1.fastq.gz ::: ${i}

#This command removes all the resulting blank files from the lists above, (it tries to match the files that dont work with
# the *_L00* code and vice-versa) and prevents from running fastqc on a blank file, savingtime

find . -size 0 -delete


##  Run fastqc on the All files in parallel (This works on all the files)
ls *_All_R1.fastq.gz | time parallel -j+0 --eta 'fastqc {}'
 

##### Make a directory for my results in my home folder
# mkdir /home/YOUR_ID/class_shared/YOUR_NAME/group
mkdir /home/aubcls03/class_shared/Steven/group/
mkdir /home/aubcls03/class_shared/Steven/group/fastqc1



### Move the results output to my directory for safe keeping
# cp *.fastqc.zip /home/YOUR_ID/class_shared/YOUR_NAME/group/fastqc1
cp *fastqc.zip /home/aubcls03/class_shared/Steven/testing/fastqc1
cp *.txt /home/aubcls03/class_shared/Steven/testing/fastqc1

#moves back to your home directory and makes a tarball:
cd /home/aubcls03/class_shared/Steven/group/
tar -czvf fastqc1.tgz /home/aubcls03/class_shared/Steven/group/fastqc1

# moves back to scratch directory:
cd /scratch/aubcls03/group/

#------------------------------------------------------------------------------------------------------------------------------
#now for trimming:

##################  Now for the Cleaning ################################
# copy over the fasta file with the adapter file to use for screening
#I assume we use the same one
cp /home/aubcls03/class_shared/code/AdaptersToTrim.fa .

#### Create list of names (see above for description) and put into file called list.
ls *_All_R1.fastq.gz | cut -d "_" -f 1,2 | sort | uniq > list2

#Note that if the script times out you can edit this list using "tail -n _ list " to start back at where it leaves off
#(This is what I had to do when testing it) (The _ represents whatever line you want to start from)


### while loop to process through the names in the list
while read i
do
############ Trimmomatic #############
############  Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
#MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
#SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
#requiredQuality: specifies the average quality required.
# -threads  is the option to define the number of threads (cores) to use. For this to be effective you need to request those cores at subm$
#  ON HOPPER: trimmomatic-0.36

#ours are single end reads. Again, the code has to be ran to accomodate the files that dont match the _All_R1* code, and
# vice-versa
java -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 6 -phred33 "$i"_All_R1.fastq.gz "$i"_All_R1_trimmed.fastq.gz ILLUMINACLIP:AdaptersToTrim.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36 
java -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 6 -phred33 "$i"_CTTGTA_All_R1.fastq.gz "$i"_CTTGTA_All_R1_trimmed.fastq.gz ILLUMINACLIP:AdaptersToTrim.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

done<list2

#again, remove any empty files before running fastqc

find . -size 0 -delete


############### Now assess Quality again
#fastqc on the cleaned paired fastq files in parallel
ls *_All_R1_trimmed.fastq.gz | parallel -j+0  --eta 'fastqc {}'



##### Make a directory for my results in my home folder
# mkdir /home/YOUR_ID/class_shared/YOUR_NAME/group/trimmed
mkdir /home/aubcls03/class_shared/Steven/group/trimmed

### Move the results output to my directory for safe keeping
mv *trimmed_fastqc.zip /home/aubcls03/class_shared/Steven/group/trimmed

#moves to group directory and makes tarball:
cd /home/aubcls03/class_shared/Steven/group/
tar -czvf trimmed.tgz /home/aubcls03/class_shared/Steven/group/trimmed

#moves back to scratch:
cd /scratch/aubcls03/group/

#----------------------------------------------------------------------------------------------------------------------------------
#Now for mapping
# Copy the reference genome to your scratch directory
ls /apps/bio/unzipped/ensembl/pub/release-87/anolis_carolinensis/dna/*.dna.toplevel.fa | time parallel -j+0 --eta 'cp {} .'

echo "The copying is complete..."

#Indexing reference library for BWA mapping: (I had to index it using the bwtsw command below as it kept giving 
#error messages when running the (is) way)
        # -p is the prefix
        # -a is the algorithm (is) then the input file
bwa index -p  Anolis_carolinensis_genome -a bwtsw  Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa



#Setting up for loop to run the mapping: (I couldn't get it to work with her code, so I had to look
#up how to do it and then I could only get it to work by using the bwa code below:

ls | grep "trimmed.fastq.gz" |cut -d "_" -f 1,2 | sort | uniq > list6
numfiles=`less list6 | wc -l`
files=($(ls | grep "trimmed.fastq.gz" |cut -d "_" -f 1,2 | sort | uniq))

#you can change the i=0 to whatever file this stops on so it will pick back up running where it leaves off!
for ((i=0; i<=$numfiles-1; i++))
do 
echo Now processing : ${files[i]} ...

bwa mem -t 20  Anolis_carolinensis_genome ${files[i]}_All_R1_trimmed.fastq.gz > ${files[i]}.sam
bwa mem -t 20  Anolis_carolinensis_genome ${files[i]}_CTTGTA_All_R1_trimmed.fastq.gz > ${files[i]}_in.sam


done 

# Output =  2_TGACCA.sam
# Output =  10-In_CTTGTA.sam
# Output =  1_In_in.sam
#will give several blank files as a result of not being able to find the ones matching the _in files
#again, remove any empty files

find . -size 0 -delete


###############  Using samtools to process the bam Input: HS06_GATCAG_All.sam
#  Make the list of prefixes for all the .sam files we want to process with Samtools
ls | grep "trimmed.fastq.gz" |cut -d "_" -f 1,2 | sort | uniq  > list3

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
done<list3

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
done<list3

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
#mkdir /home/aubcls03/class_shared/Steven/group/BWA_Counts
mkdir /home/aubcls03/class_shared/Steven/group/BWA_Counts
mkdir /home/aubcls03/class_shared/Steven/group/BWA_Stats

### Move the results output to my directory for safe keeping
## /home/YOUR_ID/class_shared/YOUR_NAME/group/BWA_Counts
#  *** EDIT ***
cp *Counts.txt /home/aubcls03/class_shared/Steven/group/BWA_Counts
cp *Stats.txt /home/aubcls03/class_shared/Steven/group/BWA_Stats

# Move to your home directory and tarball the folders so they are ready to transfer back to your computer.
#  *** EDIT ***
cd /home/aubcls03/class_shared/Steven/group/
tar -cvzf BWA_Counts.tar.gz /home/aubcls03/class_shared/Steven/group/BWA_Counts
tar -cvzf BWA_Stats.tar.gz /home/aubcls03/class_shared/Steven/group/BWA_Stats

#----------------------------------------------------------------------------------------
#runs final clean up for files we dont need from the samtools portion of the code:

#move back to your scratch directory
#cd /scratch/YOURID/group/
cd /scratch/aubcls03/group/

ls -s *_sorted* | sort | head -n 24 | cut -d " " -f 8  > list12



#The loop that deletes the unnecessary files:
while read z;
do
rm ${z}

done<list12

#now we are left with only the files from the trimming,mapping,samtools steps we need!

#-------------------------------------------------------------------------
#now for the stats part...

