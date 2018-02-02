#/bin/sh

#-- Auburn University High Performance and Parallel Computing
#-- Hopper Cluster Sample Job Submission Script

#-- This script provides the basic scheduler directives you
#-- can use to submit a job to the Hopper scheduler.
#-- Other than the last two lines it can be used as-is to
#-- send a single node job to the cluster. Normally, you
#-- will want to modify the #PBS directives below to reflect
#-- your workflow...

####-- For convenience, give your job a name

#PBS -N Trinotate_22Sept

#-- Provide an estimated wall time in which to run your job
#-- The format is DD:HH:MM:SS.  

#PBS -l walltime=03:00:00:00 

#-- Indicate if\when you want to receive email about your job
#-- The directive below sends email if the job is (a) aborted, 
#-- when it (b) begins, and when it (e) ends

#PBS -m abe tss0019@auburn.edu

#-- Inidicate the working directory path to be used for the job.
#-- If the -d option is not specified, the default working directory 
#-- is the home directory. Here, we set the working directory
#-- current directory

#PBS -d /home/tss0019/GarterSnakes/Annotation/Reference_set_26782_2013_11

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional

#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.

#PBS -l nodes=1:ppn=20

#-- Now issue the commands that you want to run on the compute nodes.

#-- With the -V option, you can load any software modules
#-- either before submitting, or in the job submission script.

#-- You should modify the lines below to reflect your own
#-- workflow...

#module load <myprogram_modulefile>
#module load fastqc/11.5
#module load trimmomatic/0.36
#module load samtools/1.3.1
#module load gcc/5.1.0
#module load bowtie2/2.2.9
#module load tophat/2.1.1
#module bcftools/1.3.1
#module load hisat/2.0.5
module load stringtie/1.3.2d
module load trinotate/3.0.2
module load sqlite/3190300
module load ncbi-blast/2.6.0
module load hmmer/3.0
module load pfamscan/9592
module load transdecoder/3.0.1
module load signalp-4.1/signalp
module load signalp/4.1/signalp
module load signalp/4.1
module load tmhmm/2.0c
module load rnammer/1.2
module load hmmer/3.1b2
module load perl/5.26.0

#./myprogram <parameters>

#--  After saving this script, you can submit your job to the queue
#--  with...

#--  qsub sample_job.sh
##########################################


#PBS -j oe
#PBS -q debug

# Define DATADIR to be where the input files are
DATADIR=/home/tss0019/GarterSnakes/Annotation/Reference_set_26782_2013_11

#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x
#
