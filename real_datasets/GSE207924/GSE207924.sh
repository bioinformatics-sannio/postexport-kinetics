mkdir sequencing_human
cd sequencing_human
for i in {270..279} {293..302}; do
prefetch SRR20080${i}
fasterq-dump SRR20080${i}
done
cd ..

mkdir sequencing_mouse
cd sequencing_mouse
for i in {344..353} {324..333}; do
prefetch SRR20080${i}
fasterq-dump SRR20080${i}
done
cd ..

mkdir cutadapt_human
for i in {270..279} {293..302}; do
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -e .2 \
         -m 20 \
         -u 3 \
         -j 8 \
         -o cutadapt_human/SRR20080${i}_noadapt_1.fastq \
         -p cutadapt_human/SRR20080${i}_noadapt_2.fastq \
         sequencing_human/SRR20080${i}_1.fastq sequencing_human/SRR20080${i}_2.fastq

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -q 20,0 \
         --nextseq-trim=20 \
         -m 20 \
         -j 8 \
         -o cutadapt_human/SRR20080${i}_trimmed_1.fastq \
         -p cutadapt_human/SRR20080${i}_trimmed_2.fastq \
         cutadapt_human/SRR20080${i}_noadapt_1.fastq cutadapt_human/SRR20080${i}_noadapt_2.fastq

rm cutadapt_human/SRR20080${i}_noadapt_1.fastq
rm cutadapt_human/SRR20080${i}_noadapt_2.fastq

cutadapt -j 8 \
         -m 20 \
         -U -3 \
         --max-n 0 \
         -o cutadapt_human/SRR20080${i}_1.fastq \
         -p cutadapt_human/SRR20080${i}_2.fastq \
         cutadapt_human/SRR20080${i}_trimmed_1.fastq cutadapt_human/SRR20080${i}_trimmed_2.fastq

rm cutadapt_human/SRR20080${i}_trimmed_1.fastq
rm cutadapt_human/SRR20080${i}_trimmed_2.fastq
done
wait

mkdir cutadapt_mouse
for i in {344..353} {324..333}; do
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -e .2 \
         -m 20 \
         -u 3 \
         -j 8 \
         -o cutadapt_mouse/SRR20080${i}_noadapt_1.fastq \
         -p cutadapt_mouse/SRR20080${i}_noadapt_2.fastq \
         sequencing_mouse/SRR20080${i}_1.fastq sequencing_mouse/SRR20080${i}_2.fastq

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -q 20,0 \
         --nextseq-trim=20 \
         -m 20 \
         -j 8 \
         -o cutadapt_mouse/SRR20080${i}_trimmed_1.fastq \
         -p cutadapt_mouse/SRR20080${i}_trimmed_2.fastq \
         cutadapt_mouse/SRR20080${i}_noadapt_1.fastq cutadapt_mouse/SRR20080${i}_noadapt_2.fastq

rm cutadapt_mouse/SRR20080${i}_noadapt_1.fastq
rm cutadapt_mouse/SRR20080${i}_noadapt_2.fastq

cutadapt -j 8 \
         -m 20 \
         -U -3 \
         --max-n 0 \
         -o cutadapt_mouse/SRR20080${i}_1.fastq \
         -p cutadapt_mouse/SRR20080${i}_2.fastq \
         cutadapt_mouse/SRR20080${i}_trimmed_1.fastq cutadapt_mouse/SRR20080${i}_trimmed_2.fastq

rm cutadapt_mouse/SRR20080${i}_trimmed_1.fastq
rm cutadapt_mouse/SRR20080${i}_trimmed_2.fastq
done
wait

mkdir genome_human
cd genome_human
wget https://zenodo.org/records/11356463/files/K562_ensGRCh38_dm6_ercc_cat.fasta.gz
wget https://zenodo.org/records/11356463/files/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf.gz
cd ..

mkdir genome_mouse
cd genome_mouse
wget https://zenodo.org/records/11356463/files/NIH3T3_mm10_dm6_ercc_cat.fasta.gz
wget https://zenodo.org/records/11356463/files/NIH3T3_mm10_MTmod_dm6_ercc_cat.gtf.gz
cd ..

mkdir index_human
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir index_human --genomeFastaFiles genome_human/K562_ensGRCh38_dm6_ercc_cat.fasta --sjdbGTFfile genome_human/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf --sjdbOverhang 147

mkdir index_mouse
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir index_mouse --genomeFastaFiles genome_mouse/NIH3T3_mm10_dm6_ercc_cat.fasta --sjdbGTFfile genome_mouse/NIH3T3_mm10_MTmod_dm6_ercc_cat.gtf --sjdbOverhang 147

grep ">h_" genome_human/K562_ensGRCh38_dm6_ercc_cat.fasta | grep -v "h_MT" | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_human/K562_h_noMT_chrNameLength.bed
grep ">h_MT" genome_human/K562_ensGRCh38_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_human/K562_h_MTonly_chrNameLength.bed
grep ">d_" genome_human/K562_ensGRCh38_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_human/K562_d_chrNameLength.bed
grep ">e_" genome_human/K562_ensGRCh38_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_human/K562_e_chrNameLength.bed

grep ">m_" genome_mouse/NIH3T3_mm10_dm6_ercc_cat.fasta | grep -v "m_MT" | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_mouse/NIH3T3_m_noMT_chrNameLength.bed
grep ">m_MT" genome_mouse/NIH3T3_mm10_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_mouse/NIH3T3_m_MTonly_chrNameLength.bed
grep ">d_" genome_mouse/NIH3T3_mm10_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_mouse/NIH3T3_d_chrNameLength.bed
grep ">e_" genome_mouse/NIH3T3_mm10_dm6_ercc_cat.fasta | sed 's/>//' | awk '{split($2,a,"-"); print $1 "\t0\t" a[2]}' > genome_mouse/NIH3T3_e_chrNameLength.bed

mkdir aligned_human
for i in {270..279} {293..302}; do
STAR --runMode alignReads --runThreadN 8 --genomeDir index_human --outSAMtype BAM Unsorted --outFilterMismatchNmax 15 --outFilterMismatchNoverReadLmax 0.09 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --alignEndsType Local --outSAMattributes NM MD NH --outSAMstrandField intronMotif --outFileNamePrefix aligned_human/SRR20080${i}_ --readFilesIn cutadapt_human/SRR20080${i}_1.fastq cutadapt_human/SRR20080${i}_2.fastq
done
wait

mkdir aligned_mouse
for i in {344..353} {324..333}; do
STAR --runMode alignReads --runThreadN 8 --genomeDir index_mouse --outSAMtype BAM Unsorted --outFilterMismatchNmax 15 --outFilterMismatchNoverReadLmax 0.09 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --alignEndsType Local --outSAMattributes NM MD NH --outSAMstrandField intronMotif --outFileNamePrefix aligned_mouse/SRR20080${i}_ --readFilesIn cutadapt_mouse/SRR20080${i}_1.fastq cutadapt_mouse/SRR20080${i}_2.fastq
done
wait

mkdir filtered_human
for i in {270..279} {293..302}; do
samtools view -@ 8 -F 0x100 -o filtered_human/SRR20080${i}_tmp1.bam aligned_human/SRR20080${i}_Aligned.out.bam
samtools view -@ 8 -F 0x400 -o filtered_human/SRR20080${i}_tmp2.bam filtered_human/SRR20080${i}_tmp1.bam
samtools view -@ 8 -F 0x800 -o filtered_human/SRR20080${i}_tmp3.bam filtered_human/SRR20080${i}_tmp2.bam
samtools view -@ 8 -F 0x8 -o filtered_human/SRR20080${i}.bam filtered_human/SRR20080${i}_tmp3.bam
rm filtered_human/SRR20080${i}_tmp1.bam
rm filtered_human/SRR20080${i}_tmp2.bam
rm filtered_human/SRR20080${i}_tmp3.bam
done
wait

mkdir filtered_mouse
for i in {344..353} {324..333}; do
samtools view -@ 8 -F 0x100 -o filtered_mouse/SRR20080${i}_tmp1.bam aligned_mouse/SRR20080${i}_Aligned.out.bam
samtools view -@ 8 -F 0x400 -o filtered_mouse/SRR20080${i}_tmp2.bam filtered_mouse/SRR20080${i}_tmp1.bam
samtools view -@ 8 -F 0x800 -o filtered_mouse/SRR20080${i}_tmp3.bam filtered_mouse/SRR20080${i}_tmp2.bam
samtools view -@ 8 -F 0x8 -o filtered_mouse/SRR20080${i}.bam filtered_mouse/SRR20080${i}_tmp3.bam
rm filtered_mouse/SRR20080${i}_tmp1.bam
rm filtered_mouse/SRR20080${i}_tmp2.bam
rm filtered_mouse/SRR20080${i}_tmp3.bam
done
wait

mkdir divided_human
for i in {270..279} {293..302}; do
samtools view -@ 8 -b -L genome_human/K562_h_noMT_chrNameLength.bed filtered_human/SRR20080${i}.bam > divided_human/SRR20080${i}_h.bam
samtools view -@ 8 -b -L genome_human/K562_h_MTonly_chrNameLength.bed filtered_human/SRR20080${i}.bam  > divided_human/SRR20080${i}_mt.bam
samtools view -@ 8 -b -L genome_human/K562_d_chrNameLength.bed filtered_human/SRR20080${i}.bam  > divided_human/SRR20080${i}_d.bam
samtools view -@ 8 -b -L genome_human/K562_e_chrNameLength.bed filtered_human/SRR20080${i}.bam  > divided_human/SRR20080${i}_e.bam
done
wait

mkdir divided_mouse
for i in {344..353} {324..333}; do
samtools view -@ 8 -b -L genome_mouse/NIH3T3_m_noMT_chrNameLength.bed filtered_mouse/SRR20080${i}.bam > divided_mouse/SRR20080${i}_m.bam
samtools view -@ 8 -b -L genome_mouse/NIH3T3_m_MTonly_chrNameLength.bed filtered_mouse/SRR20080${i}.bam  > divided_mouse/SRR20080${i}_mt.bam
samtools view -@ 8 -b -L genome_mouse/NIH3T3_d_chrNameLength.bed filtered_mouse/SRR20080${i}.bam  > divided_mouse/SRR20080${i}_d.bam
samtools view -@ 8 -b -L genome_mouse/NIH3T3_e_chrNameLength.bed filtered_mouse/SRR20080${i}.bam  > divided_mouse/SRR20080${i}_e.bam
done
wait

mkdir sorted_human
for i in {270..279} {293..302}; do
samtools sort -@ 8 -o sorted_human/SRR20080${i}_h.bam divided_human/SRR20080${i}_h.bam
samtools index -@ 8 sorted_human/SRR20080${i}_h.bam
samtools sort -@ 8 -o sorted_human/SRR20080${i}_mt.bam divided_human/SRR20080${i}_mt.bam
samtools index -@ 8 sorted_human/SRR20080${i}_mt.bam
samtools sort -@ 8 -o sorted_human/SRR20080${i}_d.bam divided_human/SRR20080${i}_d.bam
samtools index -@ 8 sorted_human/SRR20080${i}_d.bam
samtools sort -@ 8 -o sorted_human/SRR20080${i}_e.bam divided_human/SRR20080${i}_e.bam
samtools index -@ 8 sorted_human/SRR20080${i}_e.bam
done
wait

mkdir sorted_mouse
for i in {344..353} {324..333}; do
samtools sort -@ 8 -o sorted_mouse/SRR20080${i}_m.bam divided_mouse/SRR20080${i}_m.bam
samtools index -@ 8 sorted_mouse/SRR20080${i}_m.bam
samtools sort -@ 8 -o sorted_mouse/SRR20080${i}_mt.bam divided_mouse/SRR20080${i}_mt.bam
samtools index -@ 8 sorted_mouse/SRR20080${i}_mt.bam
samtools sort -@ 8 -o sorted_mouse/SRR20080${i}_d.bam divided_mouse/SRR20080${i}_d.bam
samtools index -@ 8 sorted_mouse/SRR20080${i}_d.bam
samtools sort -@ 8 -o sorted_mouse/SRR20080${i}_e.bam divided_mouse/SRR20080${i}_e.bam
samtools index -@ 8 sorted_mouse/SRR20080${i}_e.bam
done
wait

echo -e "name\tcount" > human_ercc.tsv
for i in divided_human/*_e.bam; do
    nome=$(basename "$i" _e.bam)
    count=$(samtools view -c "$i")
    echo -e "$nome\t$count" >> human_ercc.tsv
done

echo -e "name\tcount" > mouse_ercc.tsv
for i in divided_mouse/*_e.bam; do
    nome=$(basename "$i" _e.bam)
    count=$(samtools view -c "$i")
    echo -e "$nome\t$count" >> mouse_ercc.tsv
done

mkdir rmats_human
cat <<EOF >"rmats_human/nuc_0.txt"
sorted_human/SRR20080279_h.bam,sorted_human/SRR20080302_h.bam
EOF
cat <<EOF >"rmats_human/nuc_15.txt"
sorted_human/SRR20080278_h.bam,sorted_human/SRR20080301_h.bam
EOF
cat <<EOF >"rmats_human/nuc_30.txt"
sorted_human/SRR20080277_h.bam,sorted_human/SRR20080300_h.bam
EOF
cat <<EOF >"rmats_human/nuc_60.txt"
sorted_human/SRR20080276_h.bam,sorted_human/SRR20080299_h.bam
EOF
cat <<EOF >"rmats_human/nuc_120.txt"
sorted_human/SRR20080275_h.bam,sorted_human/SRR20080298_h.bam
EOF
cat <<EOF >"rmats_human/cyt_0.txt"
sorted_human/SRR20080274_h.bam,sorted_human/SRR20080297_h.bam
EOF
cat <<EOF >"rmats_human/cyt_15.txt"
sorted_human/SRR20080273_h.bam,sorted_human/SRR20080296_h.bam
EOF
cat <<EOF >"rmats_human/cyt_30.txt"
sorted_human/SRR20080272_h.bam,sorted_human/SRR20080295_h.bam
EOF
cat <<EOF >"rmats_human/cyt_60.txt"
sorted_human/SRR20080271_h.bam,sorted_human/SRR20080294_h.bam
EOF
cat <<EOF >"rmats_human/cyt_120.txt"
sorted_human/SRR20080270_h.bam,sorted_human/SRR20080293_h.bam
EOF

for i in 0 15 30 60 120; do
rmats.py --b1 rmats_human/nuc_${i}.txt \
         --b2 rmats_human/cyt_${i}.txt \
         --readLength 148 \
         --variable-read-length \
         --gtf genome_human/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf \
         -t paired \
         --libType fr-secondstrand \
         --nthread 10 \
         --statoff \
         --od rmats_human/${i}_out \
         --tmp rmats_human/${i}_tmp &
done
wait

mkdir rmats_mouse
cat <<EOF >"rmats_mouse/nuc_0.txt"
sorted_mouse/SRR20080353_m.bam,sorted_mouse/SRR20080333_m.bam
EOF
cat <<EOF >"rmats_mouse/nuc_15.txt"
sorted_mouse/SRR20080352_m.bam,sorted_mouse/SRR20080332_m.bam
EOF
cat <<EOF >"rmats_mouse/nuc_30.txt"
sorted_mouse/SRR20080351_m.bam,sorted_mouse/SRR20080331_m.bam
EOF
cat <<EOF >"rmats_mouse/nuc_60.txt"
sorted_mouse/SRR20080350_m.bam,sorted_mouse/SRR20080330_m.bam
EOF
cat <<EOF >"rmats_mouse/nuc_120.txt"
sorted_mouse/SRR20080349_m.bam,sorted_mouse/SRR20080329_m.bam
EOF
cat <<EOF >"rmats_mouse/cyt_0.txt"
sorted_mouse/SRR20080348_m.bam,sorted_mouse/SRR20080328_m.bam
EOF
cat <<EOF >"rmats_mouse/cyt_15.txt"
sorted_mouse/SRR20080347_m.bam,sorted_mouse/SRR20080327_m.bam
EOF
cat <<EOF >"rmats_mouse/cyt_30.txt"
sorted_mouse/SRR20080346_m.bam,sorted_mouse/SRR20080326_m.bam
EOF
cat <<EOF >"rmats_mouse/cyt_60.txt"
sorted_mouse/SRR20080345_m.bam,sorted_mouse/SRR20080325_m.bam
EOF
cat <<EOF >"rmats_mouse/cyt_120.txt"
sorted_mouse/SRR20080344_m.bam,sorted_mouse/SRR20080324_m.bam
EOF

for i in 0 15 30 60 120; do
rmats.py --b1 rmats_mouse/nuc_${i}.txt \
         --b2 rmats_mouse/cyt_${i}.txt \
         --readLength 148 \
         --variable-read-length \
         --gtf genome_mouse/NIH3T3_mm10_MTmod_dm6_ercc_cat.gtf \
         -t paired \
         --libType fr-secondstrand \
         --nthread 10 \
         --statoff \
         --od rmats_mouse/${i}_out \
         --tmp rmats_mouse/${i}_tmp &
done
wait

mkdir human
cp rmats_human/0_out/RI.MATS.JC.txt human/t0.txt
cp rmats_human/15_out/RI.MATS.JC.txt human/t15.txt
cp rmats_human/30_out/RI.MATS.JC.txt human/t30.txt
cp rmats_human/60_out/RI.MATS.JC.txt human/t60.txt
cp rmats_human/120_out/RI.MATS.JC.txt human/t120.txt

mkdir mouse
cp rmats_mouse/0_out/RI.MATS.JC.txt mouse/t0.txt
cp rmats_mouse/15_out/RI.MATS.JC.txt mouse/t15.txt
cp rmats_mouse/30_out/RI.MATS.JC.txt mouse/t30.txt
cp rmats_mouse/60_out/RI.MATS.JC.txt mouse/t60.txt
cp rmats_mouse/120_out/RI.MATS.JC.txt mouse/t120.txt

cat <<EOF >"human_script.R"

mat = lapply(list.files("human"), function(x){
  tab = read.table(paste0("human/", x), header = T)[, c(2, 4, 5, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18)]
  tab = tab[grep("ENSG", tab$GeneID), ]
  tab$chr = sub("^(...)..(.*)", "\\1\\2", tab$chr)
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

names(mat) = sapply(strsplit(list.files("human"), split = "\\."), function(x){x[1]})

nam = Reduce(intersect, lapply(mat, rownames))

mat = cbind(mat[[1]][nam, ], mat[[3]][nam, 10:17], mat[[4]][nam, 10:17], mat[[5]][nam, 10:17], mat[[2]][nam, 10:17])
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

# Normalizzo per profondità di sequenziamento
nor = read_tsv("human_ercc.tsv")
nor = nor$counts[c(10, 20, 10, 20, 5, 15, 5, 15,
                   9, 19, 9, 19, 4, 14, 4, 14,
                   8, 18, 8, 18, 3, 13, 3, 13,
                   7, 17, 7, 17, 2, 12, 2, 12,
                   6, 16, 6, 16, 1, 11, 1, 11)]

for(i in 1:length(nor)){mat[, i + 9] = (mat[, i + 9] / nor[i]) * 1e6}
rm(nor, i)

mat = mat[, -c(8,9)]

mat = data.frame(
  gene = rep(mat$GeneID, each = 40),
  chr = rep(mat$chr, each = 40),
  strand = rep(mat$strand, each = 40),
  upstreamES = rep(mat$upstreamES, each = 40),
  upstreamEE = rep(mat$upstreamEE, each = 40),
  downstreamES = rep(mat$downstreamES, each = 40),
  downstreamEE = rep(mat$downstreamEE, each = 40),
  timepoint = rep(c("t0", "t15", "t30", "t60", "t120"), each = 8),
  compartment = rep(c("nuc", "cyt"), each = 4),
  version = rep(c("inclusion", "skipping"), each = 2),
  replicate = c("r1", "r2"),
  count = as.numeric(as.vector(t(mat[, 8:ncol(mat)])))
)

write.csv(mat, "gse207924_human_normalized.csv", row.names = F)

EOF

cat <<EOF >"mouse_script.R"

mat = lapply(list.files("mouse"), function(x){
  tab = read.table(paste0("mouse/", x), header = T)[, c(2, 4, 5, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18)]
  tab = tab[grep("ENSM", tab$GeneID), ]
  tab$chr = sub("^(...)..(.*)", "\\1\\2", tab$chr)
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

names(mat) = sapply(strsplit(list.files("mouse"), split = "\\."), function(x){x[1]})

nam = Reduce(intersect, lapply(mat, rownames))

mat = cbind(mat[[1]][nam, ], mat[[3]][nam, 10:17], mat[[4]][nam, 10:17], mat[[5]][nam, 10:17], mat[[2]][nam, 10:17])
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

nor = read_tsv("mouse_ercc.tsv")
nor = nor$counts[c(20, 10, 20, 10, 15, 5, 15, 5,
                   19, 9, 19, 9, 14, 4, 14, 4,
                   18, 8, 18, 8, 13, 3, 13, 3,
                   17, 7, 17, 7, 12, 2, 12, 2,
                   16, 6, 16, 6, 11, 1, 11, 1)]

for(i in 1:length(nor)){mat[, i + 9] = (mat[, i + 9] / nor[i]) * 1e6}
rm(nor, i)

mat = mat[, -c(8,9)]

mat = data.frame(
  gene = rep(mat$GeneID, each = 40),
  chr = rep(mat$chr, each = 40),
  strand = rep(mat$strand, each = 40),
  upstreamES = rep(mat$upstreamES, each = 40),
  upstreamEE = rep(mat$upstreamEE, each = 40),
  downstreamES = rep(mat$downstreamES, each = 40),
  downstreamEE = rep(mat$downstreamEE, each = 40),
  timepoint = rep(c("t0", "t15", "t30", "t60", "t120"), each = 8),
  compartment = rep(c("nuc", "cyt"), each = 4),
  version = rep(c("inclusion", "skipping"), each = 2),
  replicate = c("r1", "r2"),
  count = as.numeric(as.vector(t(mat[, 8:ncol(mat)])))
)

write.csv(mat, "gse207924_mouse_normalized.csv", row.names = F)

EOF

Rscript human_script.R
Rscript mouse_script.R
