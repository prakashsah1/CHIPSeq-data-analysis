for sample in *.fastq.gz; do
    base=$(basename "$sample" .fastq.gz)
    bwa mem -t 4 ../reference/hg19.fa "$sample" | samtools view -Sb - > "${base}.bam"
done               
