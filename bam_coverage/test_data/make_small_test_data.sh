egrep "contig_7|contig_5" aln.sam > contig_7.sam
samtools view -h -o contig_7.bam contig_7.sam 
