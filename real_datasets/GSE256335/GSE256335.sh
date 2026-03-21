mkdir sequencing
cd sequencing
for i in {44..47} {67..74} {77..85} {91..99}; do
prefetch SRR280468${i}
fasterq-dump SRR280468${i}
done
cd ..

mkdir genome_mouse
cd genome_mouse
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
gzip -d GRCm39.primary_assembly.genome.fa.gz
gzip -d gencode.vM38.primary_assembly.annotation.gtf.gz
cd ..

mkdir star_index
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome_mouse/GRCm39.primary_assembly.genome.fa --sjdbGTFfile genome_mouse/gencode.vM38.primary_assembly.annotation.gtf --sjdbOverhang 75

mkdir aligned_mouse
for i in {44..47} {67..74} {77..85} {91..99}; do
STAR --outFilterMismatchNmax 3 \
     --twopassMode Basic \
     --alignEndsType EndToEnd \
     --runThreadN 40 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --alignSJDBoverhangMin 1 \
     --alignIntronMax 299999 \
     --genomeDir star_index \
     --sjdbGTFfile genome_mouse/gencode.vM38.primary_assembly.annotation.gtf \
     --outFileNamePrefix aligned_mouse/SRR280468${i} \
     --readFilesIn sequencing/SRR280468${i}_1.fastq sequencing/SRR280468${i}_2.fastq
done

mkdir rsem_reference
rsem-prepare-reference -p 40 --gtf genome_mouse/gencode.vM38.primary_assembly.annotation.gtf --star genome_mouse/GRCm39.primary_assembly.genome.fa --star-sjdboverhang 75 rsem_reference/GRCm39

mkdir rsem_counts
for i in {44..47} {67..74} {77..85} {91..99}; do
rsem-calculate-expression -p 40 \
                          --paired-end \
                          --strandedness reverse \
                          --star \
                          sequencing/SRR280468${i}_1.fastq \
                          sequencing/SRR280468${i}_2.fastq \
                          rsem_reference/GRCm39 \
                          rsem_counts/SRR280468${i}
done

mkdir rmats_output
cat <<EOF >"rmats_output/nuc_0.txt"
aligned_mouse/SRR28046883Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046884Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046885Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_30.txt"
aligned_mouse/SRR28046847Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046867Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046868Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_60.txt"
aligned_mouse/SRR28046844Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046845Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046846Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_120.txt"
aligned_mouse/SRR28046872Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046873Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046874Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/nuc_240.txt"
aligned_mouse/SRR28046869Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046870Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046871Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_0.txt"
aligned_mouse/SRR28046897Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046898Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046899Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_30.txt"
aligned_mouse/SRR28046880Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046881Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046882Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_60.txt"
aligned_mouse/SRR28046877Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046878Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046879Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_120.txt"
aligned_mouse/SRR28046894Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046895Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046896Aligned.sortedByCoord.out.bam
EOF
cat <<EOF >"rmats_output/cyt_240.txt"
aligned_mouse/SRR28046891Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046892Aligned.sortedByCoord.out.bam,aligned_mouse/SRR28046893Aligned.sortedByCoord.out.bam
EOF

for i in "0" "30" "60" "120" "240"; do
rmats.py --b1 rmats_output/nuc_${i}.txt \
         --b2 rmats_output/cyt_${i}.txt \
         --readLength 75 \
         --variable-read-length \
         --gtf genome_mouse/gencode.vM38.primary_assembly.annotation.gtf \
         -t paired \
         --libType fr-firststrand \
         --nthread 40 \
         --statoff \
         --od rmats_output/${i}_out \
         --tmp rmats_output/${i}_tmp
done

mkdir rmats_counts
cp rmats_output/0_out/RI.MATS.JC.txt rmats_counts/time_0.txt
cp rmats_output/30_out/RI.MATS.JC.txt rmats_counts/time_30.txt
cp rmats_output/60_out/RI.MATS.JC.txt rmats_counts/time_60.txt
cp rmats_output/120_out/RI.MATS.JC.txt rmats_counts/time_120.txt
cp rmats_output/240_out/RI.MATS.JC.txt rmats_counts/time_240.txt

cat <<EOF >"rmats_output/r_script.R"

mat = lapply(list.files("rmats_counts"), function(x){
  tab = read.table(paste0("rmats_counts/", x), header = T)[, c(2, 4, 5, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18)]
  tab$upstreamES = tab$upstreamES + 1
  tab$downstreamES = tab$downstreamES + 1
  rownames(tab) = apply(tab, 1, function(y){paste(y[1:7], collapse = "_")})
  inu = do.call(rbind, strsplit(tab$IJC_SAMPLE_1, ","))
  colnames(inu) = paste0("inclusionNU_", 1:3)
  snu = do.call(rbind, strsplit(tab$SJC_SAMPLE_1, ","))
  colnames(snu) = paste0("skippingNU_", 1:3)
  icy = do.call(rbind, strsplit(tab$IJC_SAMPLE_2, ","))
  colnames(icy) = paste0("inclusionCY_", 1:3)
  scy = do.call(rbind, strsplit(tab$SJC_SAMPLE_2, ","))
  colnames(scy) = paste0("skippingCY_", 1:3)
  tab = cbind(tab[, c(1:7, 12, 13)], inu, snu, icy, scy)
  return(tab)
})

names(mat) = sapply(strsplit(list.files("rmats_counts"), split = "\\."), function(x){x[1]})

nam = Reduce(intersect, lapply(mat, rownames))

mat = cbind(mat[[1]][nam, ], mat[[4]][nam, 10:21], mat[[5]][nam, 10:21], mat[[2]][nam, 10:21], mat[[3]][nam, 10:21])
mat[, 4:ncol(mat)] = lapply(mat[, 4:ncol(mat)], as.numeric)
rm(nam)

ret = apply(mat[, 10:ncol(mat)], 1, function(x){
  mea = matrix(as.numeric(x), ncol = 6, byrow = TRUE)
  mea = mea[, 1:3] + mea[, 4:6]
  mea = rowMeans(mea)
  sum(mea >= 10) >= 1
})

mat = mat[ret, ]
rm(ret)

num = sort(c(seq(10, ncol(mat), 6), seq(10, ncol(mat), 6) + 1, seq(10, ncol(mat), 6) + 2))
mat[, num] = sweep(mat[, num], 1, mat$IncFormLen, "/")

num = sort(c(seq(13, ncol(mat), 6), seq(13, ncol(mat), 6) + 1, seq(13, ncol(mat), 6) + 2))
mat[, num] = sweep(mat[, num], 1, mat$SkipFormLen, "/")

rm(num)

nor = matrix(c("SRR28046883", "SRR28046884", "SRR28046885", "SRR28046897", "SRR28046898", "SRR28046899",
               "SRR28046847", "SRR28046867", "SRR28046868", "SRR28046880", "SRR28046881", "SRR28046882",
               "SRR28046844", "SRR28046845", "SRR28046846", "SRR28046877", "SRR28046878", "SRR28046879",
               "SRR28046872", "SRR28046873", "SRR28046874", "SRR28046894", "SRR28046895", "SRR28046896",
               "SRR28046869", "SRR28046870", "SRR28046871", "SRR28046891", "SRR28046892", "SRR28046893"), 10, 3, byrow = T)

nor = cbind(nor, nor)
nor = as.vector(t(nor))

nor = sapply(nor, function(x){
  gen = read.table(paste0("rsem_counts/", x, ".isoforms.results"), header = T)
  gen = gen[gen$effective_length != 0, ]
  gen = sum(gen$expected_count / gen$effective_length)
  return(gen)
})

coe = rep(rep(c(0.39, 0.15), each = 6), 5)

for(i in 1:length(nor)){mat[, i + 9] = (mat[, i + 9] / nor[i]) * coe[i] * 1e6}
rm(nor, coe, i)

mat = mat[, -c(8,9)]

mat = data.frame(
  gene = rep(mat$GeneID, each = 60),
  chr = rep(mat$chr, each = 60),
  strand = rep(mat$strand, each = 60),
  upstreamES = rep(mat$upstreamES, each = 60),
  upstreamEE = rep(mat$upstreamEE, each = 60),
  downstreamES = rep(mat$downstreamES, each = 60),
  downstreamEE = rep(mat$downstreamEE, each = 60),
  timepoint = rep(c("t0", "t30", "t60", "t120", "t240"), each = 12),
  compartment = rep(c("nuc", "cyt"), each = 6),
  version = rep(c("inclusion", "skipping"), each = 3),
  replicate = c("r1", "r2", "r3"),
  count = as.numeric(as.vector(t(mat[, 8:ncol(mat)])))
)

write.csv(mat, "gse256335_mouse_normalized.csv", row.names = F)

EOF

Rscript r_script.R
