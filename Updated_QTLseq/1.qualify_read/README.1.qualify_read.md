File creadted by Lester Young, 12Feb2019, modified 14Mar2019

(The `1.trim_and_filter_reads` directory will replace the `1.qualify_reads` directory)-for now I'm keeping it the same

In this part of the pathway, the readfiles will be trimmed using trimmomatic and files for each individual will be 
created. A foward and reverese readfile will be made for paired and unpaired reads. We're interested in the forward 
paired and reverse paired readfiles.

The parental readfiles will be used to generate a secondary reference, if desired by the user. To implement this option, change 
`Key0_Make_secondary_reference` in `config.txt` to equal "yes". Otherwise change it to "no".

The input files for the secondary reference should be in ../parental_readfiles, as gzipped fastq files. The name of the parental 
readfiles should be the same as `${Key1_My_cultivar_sample_name}` in `config.txt` and have the following naming format:
	`${Key1_My_cultivar_sample_name}_R1.fastq.gz` for the forward reads
	`${Key1_My_cultivar_sample_name}_R2.fastq.gz` for the reverse reads

The output will be placed in 1.qualify_read/secondary_readfiles (which will be generated if necessary) and will have the following 
names
	`${Key1_My_cultivar_sample_name}-paired_[0-9]*_R1.fastq.gz` for the paired foward reads
        `${Key1_My_cultivar_sample_name}_[0-9]*_R2.fastq.gz` for the paired reverse reads
        `${Key1_My_cultivar_sample_name}_[0-9]*_R1.fastq.gz` for the unpaired foward reads
        `${Key1_My_cultivar_sample_name}_[0-9]*_R2.fastq.gz` for the unpaired reverse reads


The pipeline will expect to find the readfiles for the bulks in the ../mybulk_A and ../mybulk_B directories. The pipeline will generate 
a 1.qualify_read/mybulk_A and 1.qualify_read/mybulk_B directory if it doesn't exist. The common directory name "mybulk" and the "A" and 
"B" can be changed by the user by altering the `Key1_Bulked_sample_name`, `Key1_Bulked_sample_Type_A` and `Key1_Bulked_sample_Type_B` 
in `config.txt`. The ../mybulk_A and ../mybulk_B readfiles need to be gzipped fastq files and have the following naming format:

	`readfilename_[0-9]*_R1.fastq.gz` for the forward reads
	`readfilename_[0-9]*_R2.fastq.gz` for the reverse reads

The output of the trimmomatic will be put into 1.qualify_read/mybulk_A or 1.qualify_read/mybulk_B.
	`readfilename-paired_[0-9]*_R1.fastq.gz` for the paired foward reads
	`readfilename-paired_[0-9]*_R2.fastq.gz` for the paired reverse reads
	`readfilename-unpaired_[0-9]*_R1.fastq.gz` for the unpaired foward reads
	`readfilename-unpaired_[0-9]*_R2.fastq.gz` for the unpaired reverse reads

The A and B bulks need to have the same number of individuals in each bulk

The default filtering values for trimmomatic are set in `config.txt` and can be altered by the user. The values are:
	`Key1_Phred_quality_score_for_my_cultivar=30`	= Phred cutoff for secondary reference
	`Key1_Min_length_my_cultivar_reads=90`		= minimum of 90 bp read length
	
	`Key1_Phred_quality_score_for_bulked=30`	= Phred cutoff for bulked reads
	`Key1_Min_length_bulked_reads=90`		= minimum length of bulk reads

Trimmomatic automatically determines the phred scoring system used (assuming illumina reads are being used) as well as the number of 
threads it can use.
