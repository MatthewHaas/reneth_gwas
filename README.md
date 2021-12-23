# README for reneth_gwas
Code relating to Reneth's GWAS project

## Directory setup
In my initial analysis that I called [**nov_2021_gbs**](https://github.com/MatthewHaas/nov_2021_gbs), I created a plain text file called `nov_2021_gbs_directory_setup.txt`. It is easy enough for me to understand, but I thought the steps might make more sense to everyone else if I put it ito a format that was easier to understand. I think it is particularly important because this was done interactively and involves switching back and forth between the `bash` command line and the `R` statistical environment. Since markdown files (`.md`) enable code blocks, I thought it would help with the interactive steps outlined here.

The data are available here in the following directory:<br>
```bash
/home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/
```

If you want to count the number of `FASTQ` files there are (as a type of pre-check), you can use this command:
```bash
find /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/*fastq.gz | wc -l 
```

The analysis was carried out in the following working direcotry:
```bash
/scratch.global/haasx092/reneth_gwas
```

The first step is to make a ```CSV``` file containing the full paths to the data.
```bash
# Write the names of gzipped FASTQ files to a CSV file
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/*fastq.gz > 211117_gbs_filenames.csv
```
The resulting file (`211117_gbs_filenames.csv`) contains _all_ of the `FASTQ` files from the `Kimball_Project_008` data release. This includes Reneth's GWAS populations, the 309 & 310 populations, and Claudia's disease-resistant and -susceptible samples. For purposes of this anaysis, all non-GWAS samples were manually removed and saved as `211222_reneth_gwas_gbs_filenames.csv`.

The next step is to move into the `R` statistical environment to go from file names to a workable `CSV` file that will be used in the next step of the directory structure setup
```R
# Read in data using the data.table package
library(data.table)
fread("211222_reneth_gwas_gbs_filenames.csv", header=F) -> x

# Change the column name from V1 to something more informative (filename)
setnames(x, "filename")
# Add a new column called sample_number. It will initially contain the entire filename, but we will work to retain only the sample number
x[, sample_number := filename]
# Strip off first part of filename until sample number begins (S) but do not include it.
x[, sample_number := sub("^.+[S]", "", sample_number)]
# Strip off end of the filename (after the sample number) ... begins with "_R1" or "_R2"
x[, sample_number := sub("_[R1].+$", "", sample_number)]

# Convert sample numbers to numerical and add leading zeros to all samples (to help with sorting).
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]

# Reorder rows in ascending order
x[order(sample_number)] -> x

# Set column order (my personal preference for sample_number to come first)
setcolorder(x, c("sample_number", "filename")) -> x

# Write output to CSV
write.csv(x, file="211222_reneth_gwas_sample_names_and_numbers.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)

# Save table as an R object
save(x, file="211222_reneth_gwas_sample_names_and_numbers.Rdata")
```
After that is done, use the `CSV` file using `bash` to create the directory structure.<br>
**Note:** The `echo $i` part is not really necessary. I just included it to watch the progress.
```bash
cat 211222_reneth_gwas_sample_names_and_numbers.csv | cut -f 1 -d , \
	| while read i; do
	d=Sample_$i
	echo $i
	mkdir -p $d
	done
```
Once that is done, you will probably notice that there is a directory called `Sample_sample_number` which is an artefact of the code. I probably could change the code so that the header isn't interpreted as a sample name, but it's also super easy to just remove it after the code finishes. You can easily remove it with a one-liner:
```bash
rm -rf Sample_sample_number
```
Next, you should make a file with the list of directories. This `txt` file will come in handy for future steps of the GBS analysis.
```bash
ls Sample*/ -d | tr -d / > 211222_reneth_gwas_sample_directory_list.txt
```
This next step is necessary because we are working with paired-end reads. We are doing it because the file `211222_reneth_gwas_sample_names_and_numbers.csv` contains 2 lines per sample (one for the forward read and one for the reverse read).
```bash
awk 'FNR%2' 211222_reneth_gwas_sample_names_and_numbers.csv > 211222_reneth_gwas_file_list_every_other.csv
```
_Make sure you open the resulting file using_ `vi` _to manually remove the header line_. Once that is done, we can make symbolic links (symlinks) to point to the data rather than take up disk space by needlessly duplicating the original files. **Note** that when I analyzed the original dataset, I set `n` to start at 73 and then increased that value by 1 with each iteration of the `while` loop. Since this iteration of the analysis only contains the GWAS samples, there are gaps in the sequence of sample numbers, necessitating a different approach. The approach I used involves extracting the sample number (`Snumber`) from the file name and using that rather than relying on counting iterations through the loop.
```bash
# Make symlinks to GBS data
cat 211222_reneth_gwas_file_list_every_other.csv | cut -f 9 -d / \
	| while read i; do
	STEM=$(echo $i | rev | cut -f 3,4,5,6,7 -d "_" | rev)
	Snumber=$(echo $i | rev | cut -f 3 -d "_"| rev | sed 's/^S//g')
	n=$(printf "%04d\n" $Snumber)
	echo $STEM
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R1_001.fastq.gz Sample_$n/Sample_${n}_R1.fq.gz
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R2_001.fastq.gz Sample_$n/Sample_${n}_R2.fq.gz
	done
```
In the next step, we will move back to the `R` statistical environment to create a sample key.
```R
# Move back to R
library(data.table)

# Read in data
x <- fread("211222_reneth_gwas_sample_names_and_numbers.csv")

# Change column names
setnames(x, c("sample_number", "sample_name"))

# Add leading zeros
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]
# Add "Sample_" to each sample number
x[, sample_number := paste0("Sample_", sample_number)]

# Remove beginning the beginning part of the filename to remove the part of the path that is no longer necessary to keep
x[, sample_name := sub("^.+Project_008/", "", sample_name)]

# Remove trailing part of filenames (sample names)---ultimately, we only need one line per sample, not two (a consequence of having 2 files per sample for paired-end reads)
x[, sample_name := sub("_[R1].+$", "", sample_name)]
x[, sample_name := sub("_[R2].+$", "", sample_name)]

# Retain unique values only
x <- unique(x)

# Save to CSV
write.csv(x, file="211222_reneth_gwas_sample_key.csv", row.names = FALSE, sep=",", quote=FALSE)
```

## Adapter trimming
The next step in the process is to trim the adapters. Since this is my second time processing this dataset, there is no reason to run the FastQC quality reports.

The script to submit for the adapter trimming is [run_cutadapt.sh](adapter_trimming/run_cutadapt.sh) which depends on/calls the script [cutadapt_wrapper_script.sh](adapter_trimming/cutadapt_wrapper_script.sh). That means they need to be in the same directory in order to work properly.

## Read alignment
After you have trimmed the adapters from the reads, the next step is to align the reads to the genome. We use the Burrows-Wheeler Aligner Maximal Exact Match (BWA-MEM). I decided to speed up the step by running multiple processes in parallel; however, rather than use [GNU Parallel](https://www.gnu.org/software/parallel/), I chose to break the alignment step into 5 separate scripts, segregated by GWAS population membership (Barron, FY-C20, Itasca-C12, Itasca-C20, and K20. This is why there are 5 separate BWA scripts. **Note:** This is only appropriate for the alignment step. Don't try to do the same thing with the SNP-calling step.
1) [run_bwa_Barron.sh](alignment/run_bwa_Barron.sh) which requires [Barron_samples.txt](helper_files/Barron_samples.txt)
2) [run_bwa_FYC20.sh](alignment/run_bwa_FYC20.sh) which requires [FYC20_samples.txt](helper_files/FYC20_samples.txt)
3) [run_bwa_ItascaC12.sh](alignment/run_bwa_ItascaC12.sh) which requires [ItascaC12_samples.txt](helper_files/ItascaC12_samples.txt)
4) [run_bwa_ItascaC20.sh](alignment/run_bwa_ItascaC20.sh) which requires [ItascaC20_samples.txt](helper_files/ItascaC20_samples.txt)
5) [run_bwa_K2.sh](alignment/run_bwa_K2.sh) which requries [K2_samples.txt](helper_files/K2_samples.txt)
