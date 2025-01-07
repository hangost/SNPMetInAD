#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o ./tmp/exp.o.%j
#SBATCH --err ./tmp/exp.e.%j

ml star
ml rsem/1.3.0

samtools sort -n $i -o $sort_bam
samtools bam2fq -1 ${out_dir}/1.fq/${fq1} -2 ${out_dir}/1.fq/${fq2} ${sort_bam}
#samtools bam2fq -@ 12 -1 ${out_dir}/1.fq/${fq1} -2 ${out_dir}/1.fq/${fq2} ${sort_bam}

rm ${sort_bam}


STAR --genomeDir ${ref_dir} --sjdbFileChrStartEnd ${ref_dir}/sjdbList.fromGTF.out.tab --runThreadN 12 --outFileNamePrefix ${Foutdir} --outSAMtype BAM SortedByCoordinate --readFilesIn ${out_dir}/1.fq/${fq1} ${out_dir}/1.fq/${fq2}
samtools index ${Foutdir}/Aligned.sortedByCoord.out.bam

rsem-calculate-expression --no-bam-output --paired-end -p 12 --star ${out_dir}/1.fq/${fq1} ${out_dir}/1.fq/${fq2} /uufs/chpc.utah.edu/common/home/lee-group1/REF_hangost/4.rsem/WithOutpatch/rsem.hg19 ${Foutdir_rsem}/${samp}.result

rm ${out_dir}/1.fq/${fq1}
rm ${out_dir}/1.fq/${fq2}
