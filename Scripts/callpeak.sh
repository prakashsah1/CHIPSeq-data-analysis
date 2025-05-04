for bam in *sorted.bam
do
  if [[ $bam != mixed_input* ]]; then
     macs3 callpeak -t $bam -c mixed_input.sorted.bam -f BAM -n ${bam%.sorted.bam} -g hs --outdir peak
  fi
done
