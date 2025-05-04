for bam in *.bam; do

    base="${bam%.bam}"

    samtools sort -@ 4 -o "${base}.sorted.bam" "$bam"

    samtools index "${base}.sorted.bam"
done

