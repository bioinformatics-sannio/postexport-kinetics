mkdir sequencing
cd sequencing
for i in {502..525}; do
prefetch SRR3710${i}
done

for i in {502..525}; do
fasterq-dump SRR3710${i}
done
cd ..

mkdir fastp_drosophila
for i in {502..525}; do
fastp -i sequencing/SRR3710${i}.fastq -o fastp_drosophila/SRR3710${i}.fastq --cut_tail --cut_front --trim_poly_g -w 10 -h fastp_drosophila/report_SRR3710${i}.html -j fastp_drosophila/report_SRR3710${i}.json
done

mkdir genome
cd genome
wget https://ftp.ensembl.org/pub/current/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gtf.gz
wget https://ftp.ensembl.org/pub/current/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz
gzip -d *.gz
cd ..

mkdir star_drosophila
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir star_drosophila --genomeFastaFiles genome/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa --sjdbGTFfile genome/Drosophila_melanogaster.BDGP6.54.115.gtf --sjdbOverhang 50 --genomeSAindexNbases 12

mkdir star_saccharomyces
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir star_saccharomyces --genomeFastaFiles genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --sjdbGTFfile genome/Saccharomyces_cerevisiae.R64-1-1.115.gtf --sjdbOverhang 50 --genomeSAindexNbases 10

mkdir aligned_drosophila
for i in {502..525}; do
STAR --outFilterMismatchNmax 3 \
     --twopassMode Basic \
     --alignEndsType EndToEnd \
     --runThreadN 3 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --alignSJDBoverhangMin 1 \
     --alignIntronMax 299999 \
     --genomeDir star_drosophila \
     --sjdbGTFfile genome/Drosophila_melanogaster.BDGP6.54.115.gtf \
     --outFileNamePrefix aligned_drosophila//SRR3710${i} \
     --readFilesIn fastp_drosophila/SRR3710${i}.fastq
done
wait

mkdir aligned_saccharomyces
for i in {502..525}; do
STAR --genomeDir star_saccharomyces \
     --readFilesIn fastp_drosophila/SRR3710${i}.fastq \
     --runThreadN 10 \
     --outFileNamePrefix aligned_saccharomyces/SRR3710${i}_ \
     --outSAMtype BAM Unsorted \
     --alignIntronMax 5000
done
wait

grep "Uniquely mapped reads number" aligned_saccharomyces/*_Log.final.out | awk -F'[:\t]+' '{split($1,a,"/"); n=a[length(a)]; sub("_Log.final.out","",n); print n "\t" $NF}' > yeast_normalization.tsv

mkdir rmats_output
cat <<EOF >"rmats_output/nuc_0_0h.txt"
aligned_drosophila/SRR3710502Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710514Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_0_5h.txt"
aligned_drosophila/SRR3710503Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710515Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_1_5h.txt"
aligned_drosophila/SRR3710504Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710516Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_3_0h.txt"
aligned_drosophila/SRR3710505Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710517Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_5_0h.txt"
aligned_drosophila/SRR3710506Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710518Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_7_5h.txt"
aligned_drosophila/SRR3710507Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710519Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_0_0h.txt"
aligned_drosophila/SRR3710508Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710520Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_0_5h.txt"
aligned_drosophila/SRR3710509Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710521Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_1_5h.txt"
aligned_drosophila/SRR3710510Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710522Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_3_0h.txt"
aligned_drosophila/SRR3710511Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710523Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_5_0h.txt"
aligned_drosophila/SRR3710512Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710524Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_7_5h.txt"
aligned_drosophila/SRR3710513Aligned.sortedByCoord.out.bam,aligned_drosophila/SRR3710525Aligned.sortedByCoord.out.bam
EOF

for i in "0_0" "0_5" "1_5" "3_0" "5_0" "7_5"; do
rmats.py --b1 rmats_output/nuc_${i}h.txt \
         --b2 rmats_output/cyt_${i}h.txt \
         --readLength 51 \
         --variable-read-length \
         --gtf genome/Drosophila_melanogaster.BDGP6.54.115.gtf \
         -t single \
         --libType fr-firststrand \
         --nthread 4 \
         --statoff \
         --od rmats_output/${i}h_out \
         --tmp rmats_output/${i}h_tmp &
done
wait

mkdir drosophila
cp rmats_output/0_0h_out/RI.MATS.JC.txt drosophila/0_0h.txt
cp rmats_output/0_5h_out/RI.MATS.JC.txt drosophila/0_5h.txt
cp rmats_output/1_5h_out/RI.MATS.JC.txt drosophila/1_5h.txt
cp rmats_output/3_0h_out/RI.MATS.JC.txt drosophila/3_0h.txt
cp rmats_output/5_0h_out/RI.MATS.JC.txt drosophila/5_0h.txt
cp rmats_output/7_5h_out/RI.MATS.JC.txt drosophila/7_5h.txt

cat <<EOF >"r_script.R"

mat = lapply(list.files("drosophila"), function(x){
  tab = read.table(paste0("drosophila/", x), header = T)[, c(2, 4, 5, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18)]
  tab$upstreamES = tab$upstreamES + 1
  tab$downstreamES = tab$downstreamES + 1
  rownames(tab) = apply(tab, 1, function(y){paste(y[1:7], collapse = "_")})
  inu = do.call(rbind, strsplit(tab$IJC_SAMPLE_1, ","))
  colnames(inu) = paste0("inclusionNU_", 1:2)
  snu = do.call(rbind, strsplit(tab$SJC_SAMPLE_1, ","))
  colnames(snu) = paste0("skippingNU_", 1:2)
  icy = do.call(rbind, strsplit(tab$IJC_SAMPLE_2, ","))
  colnames(icy) = paste0("inclusionCY_", 1:2)
  scy = do.call(rbind, strsplit(tab$SJC_SAMPLE_2, ","))
  colnames(scy) = paste0("skippingCY_", 1:2)
  tab = cbind(tab[, c(1:7, 12, 13)], inu, snu, icy, scy)
  return(tab)
})

names(mat) = sapply(strsplit(list.files("drosophila"), split = "\\."), function(x){x[1]})


nam = Reduce(intersect, lapply(mat, rownames))

mat = cbind(mat[[1]][nam, ], mat[[2]][nam, 10:17], mat[[3]][nam, 10:17], mat[[4]][nam, 10:17], mat[[5]][nam, 10:17], mat[[6]][nam, 10:17])
mat[, 4:ncol(mat)] = lapply(mat[, 4:ncol(mat)], as.numeric)
rm(nam)

ret = apply(mat[, 10:ncol(mat)], 1, function(x){
  mea = matrix(as.numeric(x), ncol = 4, byrow = TRUE)
  mea = mea[, 1:2] + mea[, 3:4]
  mea = rowMeans(mea)
  sum(mea >= 10) >= 1
})

mat = mat[ret, ]
rm(ret)

num = sort(c(seq(10, ncol(mat), 4), seq(10, ncol(mat), 4) + 1))
mat[, num] = sweep(mat[, num], 1, mat$IncFormLen, "/")

num = sort(c(seq(12, ncol(mat), 4), seq(12, ncol(mat), 4) + 1))
mat[, num] = sweep(mat[, num], 1, mat$SkipFormLen, "/")

rm(num)

nor = read_tsv("yeast_normalization.tsv")
nor = nor$counts[c(1, 13, 1, 13, 7, 19, 7, 19,
                   2, 14, 2, 14, 8, 20, 8, 20,
                   3, 15, 3, 15, 9, 21, 9, 21,
                   4, 16, 4, 16, 10, 22, 10, 22,
                   5, 17, 5, 17, 11, 23, 11, 23,
                   6, 18, 6, 18, 12, 24, 12, 24)]

for(i in 1:length(nor)){mat[, i + 9] = (mat[, i + 9] / nor[i]) * 1e6}
rm(nor, i)

mat = mat[, -c(8,9)]

mat = data.frame(
  gene = rep(mat$GeneID, each = 48),
  chr = rep(mat$chr, each = 48),
  strand = rep(mat$strand, each = 48),
  upstreamES = rep(mat$upstreamES, each = 48),
  upstreamEE = rep(mat$upstreamEE, each = 48),
  downstreamES = rep(mat$downstreamES, each = 48),
  downstreamEE = rep(mat$downstreamEE, each = 48),
  timepoint = rep(c("t0_0", "t0_5", "t1_5", "t3_0", "t5_0", "t7_5"), each = 8),
  compartment = rep(c("nuc", "cyt"), each = 4),
  version = rep(c("inclusion", "skipping"), each = 2),
  replicate = c("r1", "r2"),
  count = as.numeric(as.vector(t(mat[, 8:ncol(mat)])))
)

write.csv(mat, "gse83620_drosophila_normalized.csv", row.names = F)

EOF

Rscript r_script.R
