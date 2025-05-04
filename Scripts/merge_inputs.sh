samtools merge mixed_input.bam SRR7949377.bam SRR7949383.bam
# replace inpus files names with your input file names

samtools sort -@ 4 -o  mixed_input.sorted.bam mixed_input.bam
samtools index mixed_input.sorted.ba
