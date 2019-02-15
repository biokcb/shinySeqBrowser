# Modify the input example data to be smaller and easier to work with

# just keep genes from original annotation file downloaded from wormbase
awk '$3 == "gene"' c_elegans.WS230.annotations.gff3 > cel_ws230_example.gff3

# sort an old wild type mRNA seq sample, remove unaligned, and subsample to 30% reads
samtools sort n2_hits.bam | samtools view -b -F 4 | samtools view -bs 0.3 > cel_example.bam

