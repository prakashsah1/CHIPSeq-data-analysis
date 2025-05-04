for bam in *sorted.bam
do
  base="${bam%.sorted.bam}"
  bamCoverage --bam "$bam" --normalizeUsing RPKM --extendReads 200 -o "../bigwigs/${base}.bw"
done

