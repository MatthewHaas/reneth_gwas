# README for reneth_gwas
Code relating to Reneth's GWAS project

## Directory setup
In my initial analysis that I called **nov_2021_gbs**, I created a plain text file called `nov_2021_gbs_directory_setup.txt`. It is easy enough for me to understand, but I thought the steps might make more sense to everyone else if I put it ito a format that was easier to understand. I think it is particularly important because it involves switching back and forth between the `bash` command line and the `R` statistical environment. Since markdown files (`.md`) enable code blocks, I thought it would help with the interactive steps outlined here.

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
Once that is done, you will probably notice that there is a directory called `Sample_sample_number` which is obviously a mistake. You can easily remove it with a one-liner:
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
_Make sure you open the resulting file using_ `vi` _to manually remove the header line_. Once that is done, we can make symbolic links (symlinks) to point to the data rather than take up disk space by needlessly duplicating the original files.
```bash
# Make symlinks to GBS data
# It should start at Sample_0073 because that's the first sample number in our dataset
# !!! IMPORTANT !!! I also opened the CSV file with vim (```vi``` command) and removed the header line because its presence screwed up the numbering system. Dead/broken symlinks were placed into the Sample_0073 folder and the links for Sample_0073 were placed in the Sample_0074 folder. Each directory/sample ID was therefore off by 1.
n=$(printf "%04d\n" "$((73))")
cat 211222_reneth_gwas_file_list_every_other.csv | cut -f 9 -d / \
	| while read i; do
	STEM=$(echo $i | rev | cut -f 3,4,5,6,7 -d "_" | rev)
	echo $STEM
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R1_001.fastq.gz Sample_$n/Sample_${n}_R1.fq.gz
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R2_001.fastq.gz Sample_$n/Sample_${n}_R2.fq.gz
	exitStatus=$?
	if [ $exitStatus != 0 ]; then
		n=$(printf "%04d\n" "$((10#$n+1))")
	fi
	n=$(printf "%04d\n" "$((10#$n+1))")
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
