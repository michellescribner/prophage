#!/bin/bash

# Calculate average read depth across 10bp windows for evolved populations and clones

cd /Users/mrs/Documents/prophage/coverage

### Create file with 10bp windows
samtools faidx NC_008463.fasta
awk -v OFS='\t' {'print $1,$2'} NC_008463.fasta.fai > NC_008463.txt
bedtools makewindows -g /Users/mrs/Documents/prophage/coverage/NC_008463.txt -w 10 > /Users/mrs/Documents/prophage/coverage/output/NC_008463.10.windows.bed

cd output

# Ancestor
samtools depth -a /Users/mrs/Documents/prophage/sequences/anc/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > anccov_10a.txt
paste NC_008463.10.windows.bed anccov_10a.txt > anccov_10.txt
# lasR pf5r
samtools depth -a /Users/mrs/Documents/prophage/sequences/127/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > 127cov_10a.txt
paste NC_008463.10.windows.bed 127cov_10a.txt > 127cov_10.txt
# morA pf5r
samtools depth -a /Users/mrs/Documents/prophage/sequences/128/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > 128cov_10a.txt
paste NC_008463.10.windows.bed 128cov_10a.txt > 128cov_10.txt
# lasR
samtools depth -a /Users/mrs/Documents/prophage/sequences/120/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > 120cov_10a.txt
paste NC_008463.10.windows.bed 120cov_10a.txt > 120cov_10.txt
# morA
samtools depth -a /Users/mrs/Documents/prophage/sequences/126/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > 126cov_10a.txt
paste NC_008463.10.windows.bed 126cov_10a.txt > 126cov_10.txt
#pf5r.1
samtools depth -a /Users/mrs/Documents/prophage/sequences/G1_breseq/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > G1cov_10a.txt
paste NC_008463.10.windows.bed G1cov_10a.txt > G1cov_10.txt
#pf5r.2
samtools depth -a /Users/mrs/Documents/prophage/sequences/B1_breseq/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > pf5r.2cov_10a.txt
paste NC_008463.10.windows.bed pf5r.2cov_10a.txt > pf5r.2cov_10.txt
# Planktonic Population 5
samtools depth -a /Users/mrs/Documents/prophage/sequences/P5/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > P5cov_10a.txt
paste NC_008463.10.windows.bed P5cov_10a.txt > P5cov_10.txt
# Biofilm population 1
samtools depth -a /Users/mrs/Documents/prophage/sequences/33/data/reference.bam | awk '{sum+=$3} (NR%10)==0{print sum/10; sum=0;}' > B1cov_10a.txt
paste NC_008463.10.windows.bed B1cov_10a.txt > B1cov_10.txt


# calculate ratios
samtools depth -a /Users/mrs/Documents/prophage/sequences/anc/data/reference.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' #146.644
samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/anc/data/reference.bam #146.644
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/anc/data/reference.bam #212.494

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/P5/data/reference.bam #118.348
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/P5/data/reference.bam #7295.75

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/B1/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/127/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/127/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/127/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/128/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/128/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/120/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/120/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/126/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/126/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/G1_breseq/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/G1_breseq/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/B1_breseq/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/B1_breseq/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/1306_breseq/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/1306_breseq/data/reference.bam 

samtools coverage -r NC_008463.1:1-6,537,648 /Users/mrs/Documents/prophage/sequences/1307_breseq/data/reference.bam 
samtools coverage -r NC_008463.1:4,345,126-4,355,790 /Users/mrs/Documents/prophage/sequences/1307_breseq/data/reference.bam 