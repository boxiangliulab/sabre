#!/usr/bin/env bash 
rm output.txt
for ((i=1;i<=22;i++));
do
echo Now Proccesing Chromosome $i
#python ./main.py --bam_path ./data/NA12878_WES_v2_phased_possorted_bam.bam --vcf_path ./data/NA12878_WES_v2_phased_variants.vcf.gz --sample_name NA12878_WES_v2 --restrict_chr chr$i --shortest_path >> output.txt
python main.py --bam_path ./data/Aligned.sortedByCoord.out.bam --vcf_path ./data/NA12878.vcf.gz --barcode_path ./data/barcodes.tsv --sample_name NA12878 --restrict_chr chr$i  --verbose --as_quality 0.05 --shortest_path --input_type umitools --edge_threshold 10 --mapq_threshold 60 >> output.txt
done

echo Now Proccesing Chromosome X
#python ./main.py --bam_path ./data/NA12878_WES_v2_phased_possorted_bam.bam --vcf_path ./data/NA12878_WES_v2_phased_variants.vcf.gz --sample_name NA12878_WES_v2 --restrict_chr chrX --shortest_path >> output.txt
python main.py --bam_path ./data/Aligned.sortedByCoord.out.bam --vcf_path ./data/NA12878.vcf.gz --barcode_path ./data/barcodes.tsv --sample_name NA12878 --restrict_chr chrX  --verbose --as_quality 0.05 --shortest_path --input_type umitools --edge_threshold 10 --mapq_threshold 60 >> output.txt
python ./result_statistic.py
