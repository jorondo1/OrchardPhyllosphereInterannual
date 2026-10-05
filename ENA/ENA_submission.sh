# See https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html

APPLE=/net/nfs-bio/jbod2/def-ilafores/analysis/2026_AppleMicrobiome
cd $APPLE
module purge
# cd /home/def-ilafores/programs # Install aspera latest :
# sh ibm-aspera-connect_4.2.12.780_linux_x86_64.sh 
# export PATH=/home/ronj2303/.aspera/connect/bin:$PATH # or add to bashrc
ml java

####################################
### Submitting raw (clean) samples #
####################################

mkdir -p ENA_submission/run_manifests

tail -n +2 ENA_submission/samples-2026-09-28T17_44_05.csv | while IFS= read -r line; do
	
	id=$(echo "$line" | awk -F',' '{print $1}')
    sample=$(echo "$line" | awk -F',' '{print $2}')
    
	file="ENA_submission/run_manifests/${id}_manifest.txt"
	fq1=$(find $APPLE/data/*/2_cutadapt -maxdepth 2 -type f -name "*${sample}_R1.fastq.gz")
	fq2=$(find $APPLE/data/*/2_cutadapt -maxdepth 2 -type f -name "*${sample}_R2.fastq.gz")
	md5_1=$(find $APPLE/data/*/2_cutadapt -maxdepth 2 -type f -name "*${sample}_R1.fastq.md5")
	md5_2=$(find $APPLE/data/*/2_cutadapt -maxdepth 2 -type f -name "*${sample}_R2.fastq.md5")
	
	echo -e "SAMPLE\t${id}" > "$file"
	echo -e "STUDY\tPRJEB127428" >> "$file"
	echo -e "NAME\tOrchard Phyllosphere ${sample}" >> "$file"    # EDIT THIS
	echo -e "INSTRUMENT\tIllumina MiSeq" >> "$file"              # EDIT THIS
	echo -e "INSERT_SIZE\t300" >> "$file"                        # EDIT THIS
	echo -e "LIBRARY_NAME\t2025" >> "$file"                      # EDIT THIS
	echo -e "LIBRARY_SOURCE\tMETAGENOMIC" >> "$file"
	echo -e "LIBRARY_SELECTION\tPCR" >> "$file"                  # EDIT THIS
	echo -e "LIBRARY_STRATEGY\tAMPLICON" >> "$file"              # EDIT THIS
	echo -e "FASTQ\t$fq1" >> "$file"
	echo -e "FASTQ\t$fq2" >> "$file"

done 


for man in $(find ENA_submission/run_manifests -type f -name '*manifest.txt'); do
java -jar /net/nfs-bio/jbod2/def-ilafores/programs/webin-cli-9.0.3.jar \
	-context reads -userName Webin-67053 -password LRCGsnth==1 \
	-submit -ascp -manifest "$man"
done
