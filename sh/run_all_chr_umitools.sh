#!/usr/bin/env bash 
rm output.txt
for ((i=1;i<=22;i++));
do
echo Now Proccesing Chromosome $i
python main.py --bam_path ./data/SG_HEL_B001_L002_UMITOOLS.bam --vcf_path ./data/SG_HEL_B001_L002_UMITOOLS.phased.vcf.gz --barcode_path ./data/barcodes.tsv --sample_name SG_HEL_B001_L002_UMITOOLS --restrict_chr chr$i --as_quality 0.05 --shortest_path --input_type umitools --edge_threshold 1 --mapq_threshold 60 --remove_node auto >> output.txt
done

python ./result_statistic.py
