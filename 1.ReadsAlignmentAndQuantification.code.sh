#generate the reference for STAR
#reference was downloaed from NCBI
#https://www.ncbi.nlm.nih.gov/genome/?term=txid10036[orgn]
#https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mesocricetus_auratus/all_assembly_versions/GCF_017639785.1_BCM_Maur_2.0/

#manually convert GCF_017639785.1_BCM_Maur_2.0_genomic.fna to GCF_017639785.1_BCM_Maur_2.0_genomic.fa
#sjdboverhang was set to 149 as our read = 2*150
cd /Softwares/deng/STAR/STARref/GoldenHamster
/Softwares/deng/STAR/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir BCM_Maur_2.0_index \
--genomeFastaFiles GCF_017639785.1_BCM_Maur_2.0_genomic.fa \
--sjdbGTFfile GCF_017639785.1_BCM_Maur_2.0_genomic.gtf \
--sjdbOverhang 149

#mapping the raw fastq to BCM_Maur_2.0####
/Softwares/deng/STAR/STAR-2.7.9a/bin/Linux_x86_64/STAR --genomeLoad LoadAndExit --genomeDir /Softwares/deng/STAR/STARref/GoldenHamster/BCM_Maur_2.0_index
for i in $(ls *_1.fq.gz | sed 's#_1.fq.gz##' |  sort -u); 
do 
/Softwares/deng/STAR/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 20 \
--genomeDir /Softwares/deng/STAR/STARref/GoldenHamster/BCM_Maur_2.0_index/ \
--twopassMode Basic \
--outSAMstrandField intronMotif \
--sjdbOverhang 149 \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand gunzip -c \
--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
--outFileNamePrefix /data2/deng/COVID19/STARMapping/$i
done
/Softwares/deng/STAR/STAR-2.7.9a/bin/Linux_x86_64/STAR --genomeLoad Remove --genomeDir /Softwares/deng/STAR/STARref/MouseSjdb149/


##validate the stran specific library
cd /data2/deng/COVID19/STARMapping
for i in ls *.bam;
do
infer_experiment.py -r /Softwares/deng/STAR/STARref/GoldenHamster/GCF_017639785.1_BCM_Maur_2.0_genomic.gtf.bed -i $i
done
#all samples were sequenced by strand-speciific


##obtain the count matrix
cd /data2/deng/COVID19/STARMapping
/Softwares/deng/subread-2.0.3-source/bin/featureCounts -p -t exon -g gene_id -s 2 -a /Softwares/deng/STAR/STARref/GoldenHamster/GCF_017639785.1_BCM_Maur_2.0_genomic.gtf -o ../FeatureCount/Hamster_counts.txt \
D30_BA2_2_1Aligned.sortedByCoord.out.bam \
D30_BA2_2_2Aligned.sortedByCoord.out.bam \
D30_BA2_2_3Aligned.sortedByCoord.out.bam \
D30_BA2_2_4Aligned.sortedByCoord.out.bam \
D30_BA2_2_5Aligned.sortedByCoord.out.bam \
D30_Cont_4_1Aligned.sortedByCoord.out.bam \
D30_Cont_4_2Aligned.sortedByCoord.out.bam \
D30_Cont_4_3Aligned.sortedByCoord.out.bam \
D30_Cont_4_4Aligned.sortedByCoord.out.bam \
D30_Cont_4_5Aligned.sortedByCoord.out.bam \
D30_H1N1_3_1Aligned.sortedByCoord.out.bam \
D30_H1N1_3_2Aligned.sortedByCoord.out.bam \
D30_H1N1_3_3Aligned.sortedByCoord.out.bam \
D30_H1N1_3_4Aligned.sortedByCoord.out.bam \
D30_H1N1_3_5Aligned.sortedByCoord.out.bam \
D30_WT_1_1Aligned.sortedByCoord.out.bam \
D30_WT_1_2Aligned.sortedByCoord.out.bam \
D30_WT_1_4Aligned.sortedByCoord.out.bam \
D30_WT_1_5Aligned.sortedByCoord.out.bam \
D7_BA2_2_1Aligned.sortedByCoord.out.bam \
D7_BA2_2_2Aligned.sortedByCoord.out.bam \
D7_BA2_2_3Aligned.sortedByCoord.out.bam \
D7_BA2_2_4Aligned.sortedByCoord.out.bam \
D7_BA2_2_5Aligned.sortedByCoord.out.bam \
D7Control1_1Aligned.sortedByCoord.out.bam \
D7Control1_2Aligned.sortedByCoord.out.bam \
D7Control1_3Aligned.sortedByCoord.out.bam \
D7Control1_4Aligned.sortedByCoord.out.bam \
D7Control1_5Aligned.sortedByCoord.out.bam \
D7_H1N1_3_1Aligned.sortedByCoord.out.bam \
D7_H1N1_3_2Aligned.sortedByCoord.out.bam \
D7_H1N1_3_3Aligned.sortedByCoord.out.bam \
D7_H1N1_3_4Aligned.sortedByCoord.out.bam \
D7_H1N1_3_5Aligned.sortedByCoord.out.bam \
D7SARSWT_1_1Aligned.sortedByCoord.out.bam \
D7SARSWT_1_2Aligned.sortedByCoord.out.bam \
D7SARSWT_1_3Aligned.sortedByCoord.out.bam \
D7SARSWT_1_4Aligned.sortedByCoord.out.bam \
D7SARSWT_1_5Aligned.sortedByCoord.out.bam



