#!/bin/bash

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## FASTQC output is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##		Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##			The last line of this script will make a tarball of the output directory to bring back to your computer
## For running the script on the Alabama Super Computer.
                ##For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
        ## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
        ## then run the script by using "run_script [script name]"
        ## suggested paramenters are below to submit this script.
                ## queue: medium 
                ## core: 6
                ## time limit (HH:MM:SS): 02:00:00  (may need to increase, if so run on medium queue)
                ## Memory: 12gb
                ## run on asax
###############################################

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

# Modules
	#  load the module
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

## STOP. You need to replace the [number] with YOUR paths to 
##       make variables for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsf0051                        ## Example: MyID=aubrmg001

# Variables: raw data directory (DD), working directory(WD), Quality after cleaning (PCQ), name of file containing the adpaters.
WD=/scratch/${MyID}/FinalProj                          ## Example: WD=/scratch/$MyID/PracticeRNAseq
DD=/scratch/${MyID}/FinalProj/RawData                          ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
CD=/scratch/${MyID}/FinalProj/CleanData                       ## Example: CD=/scratch/$MyID/PracticeRNAseq/CleanData
PCQ=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
				## In the future, for your data, you will likely need to edit this for other projects based on how your libraries 
				## were made to search for the correct adapters for your project
				
## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
## Note, if the parent directory exists do not use -p
mkdir ${CD}
mkdir ${WD}/${PCQ}

### making directory to store non cat SRR files after they run in the loop
SRRD=/scratch/${MyID}/FinalProj/SRRdirectory
mkdir ${SRRD}

#Move to data directory
cd ${DD}

#####Concatenating Files############
Male1=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934759" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Male1: ${#Male1[@]}"
Male2=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934760" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Male2: ${#Male2[@]}"
Male3=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934761" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Male3: ${#Male3[@]}"

Female1=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934616" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Female1: ${#Female1[@]}"
Female2=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934617" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Female2: ${#Female2[@]}"
Female3=($(awk -F "," '{gsub(/\r/,"")} $5 == "SAMN31934618" {print $1}' "${WD}/BeagleSamples.csv"))
	echo "Female3: ${#Female3[@]}"

for SRR in "${Male1[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934759_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934759_R2.fastq
done

for SRR in "${Male2[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934760_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934760_R2.fastq
done

for SRR in "${Male3[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934761_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934761_R2.fastq
done

for SRR in "${Female1[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934616_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934616_R2.fastq
done

for SRR in "${Female2[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934617_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934617_R2.fastq
done

for SRR in "${Female3[@]}"; do
        cat ${SRR}*_1*.fastq >> SAMN31934618_R1.fastq
        cat ${SRR}*_2*.fastq >> SAMN31934618_R2.fastq
done

####List out the SAMN files created#######

ls -lh SAMN*_R1.fastq SAMN*_R2.fastq

###############################################################

## Move SRR files out of RawData into SRRdirectory
mv ${DD}/SRR*.fastq ${SRRD}/

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list
	echo "list size:" $(wc -l < list)

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
cp /home/${MyID}/graze_class/AdaptersToTrim_All.fa .

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
                 #/apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar 
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
	PE -threads 6 -phred33 \
        "$i"_R1.fastq "$i"_R2.fastq  \
        ${CD}/"$i"_R1_paired.fastq ${CD}/"$i"_R1_unpaired.fastq  ${CD}/"$i"_R2_paired.fastq ${CD}/"$i"_R2_unpaired.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.

	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

#####Condifence check for trimmomatic
	 echo "Trim done: $i"
    ls ${CD}/${i}_R1_paired.fastq ${CD}/${i}_R2_paired.fastq

fastqc ${CD}/"$i"_R1_paired.fastq --outdir=${WD}/${PCQ}
fastqc ${CD}/"$i"_R2_paired.fastq --outdir=${WD}/${PCQ}

#######confidence check for fastqc
echo "FastQC done: $i"

done<list			# This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd ${WD}

########Confidence check looking at QC files made
echo "QC files:" $(ls ${WD}/${PCQ} | wc -l)

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${PCQ}.tar.gz ${PCQ}

## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate the quality of the data.

cp ${DD}/SAMN* /home/${MyID}/FinalProject
cp ${PCQ}.tar.gz /home/${MyID}/FinalProject
