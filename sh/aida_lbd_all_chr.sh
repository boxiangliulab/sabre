#!/usr/bin/env bash 
rm output.txt
for ((i=1;i<=22;i++));
do
echo Now Proccesing Chromosome $i
python main.py --bam_path /data/cy/AIDA_LBD/cellranger_outputs/LBD_SG_HEL_H331_75yo/outs/possorted_genome_bam.bam --vcf_path /data/cy/AIDA_LBD/cellranger_outputs/LBD_SG_HEL_H331_75yo/outs/variants.vcf.gz --sample_name LBD_SG_HEL_H331_75yo --restrict_chr chr$i --as_quality 0.05 --shortest_path --edge_threshold 5 --mapq_threshold 60 --remove_node auto --raw_vcf --vcf_qual 10 >> output.txt
done

python ./result_statistic.py
