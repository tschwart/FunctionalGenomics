#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load trimmomatic/0.35
module load gnu_parallel/201612222

## Recommended submission:
        # 6 cores
        # 2 hours
        # 16gb
       
       
### Change directory to the scratch directory
# cd /scratch/YOUR_ID/fastqc
cd /scratch/aubtss/fastqc_GNU_parallel3

##################  Now for the Cleaning ################################
# copy over the fasta file with the adapter file to use for screening
cp /home/aubtss/class_shared/code/AdaptersToTrim.fa .

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

ls | grep ".fastq" |cut -d "_" -f 1,2 | sort | uniq > list

### while loop to process through the names in the list
while read i
do

############ Trimmomatic #############
############  Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
#MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
#SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
#requiredQuality: specifies the average quality required.
# -threads  is the option to define the number of threads (cores) to use. For this to be effective you need to request those cores at submission
#  ON HOPPER: trimmomatic-0.36

java -jar /opt/asn/apps/trimmomatic_0.35/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 6 -phred33 "$i"_All_R1.fastq "$i"_All_R2.fastq "$i"_All_R1_paired_threads.fastq "$i"_All_R1_unpaired_threads.fastq "$i"_All_R2_paired_threads.fastq "$i"_All_R2_unpaired_threads.fastq ILLUMINACLIP:AdaptersToTrim.fa:2:30:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36 

done<list

############### Now assess Quality again
#fastqc on the cleaned paired fastq files in parallel
ls *_R1_paired_threads.fastq | parallel -j+0  --eta 'fastqc {}'
ls *_R2_paired_threads.fastq | parallel -j+0  --eta 'fastqc {}'


##### Make a directory for my results in my home folder
# mkdir /home/YOUR_ID/class_shared/YOUR_NAME/fastqc
mkdir /home/aubtss/class_shared/Tonia/GNU_parallel3

### Move the results output to my directory for safe keeping
# cp *.fastqc.zip /home/YOUR_ID/class_shared/YOUR_NAME/fastqc
cp *fastqc.zip /home/aubtss/class_shared/Tonia/GNU_parallel3
cp *.txt /home/aubtss/class_shared/Tonia/GNU_parallel3


