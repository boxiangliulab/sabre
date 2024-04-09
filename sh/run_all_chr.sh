#!/usr/bin/env bash 
rm output.txt
for ((i=1;i<=22;i++));
do
echo Now Proccesing Chromosome $i
#python ./main.py --bam ./data/NA12878_WES_v2_phased_possorted_bam.bam --vcf ./data/NA12878_WES_v2_phased_variants.vcf.gz --sample NA12878_WES_v2 --chr chr$i --shortest_path >> output.txt
python main.py --bam ./data/possorted_genome_bam.bam --vcf ./data/NA12878.vcf.gz --barcode_path ./data/barcodes.tsv --sample NA12878 --chr chr$i  --verbose --as_quality 0.01 --shortest_path >> output.txt
done

echo Now Proccesing Chromosome X
#python ./main.py --bam ./data/NA12878_WES_v2_phased_possorted_bam.bam --vcf ./data/NA12878_WES_v2_phased_variants.vcf.gz --sample NA12878_WES_v2 --chr chrX --shortest_path >> output.txt
python main.py --bam ./data/possorted_genome_bam.bam --vcf ./data/NA12878.vcf.gz --barcode_path ./data/barcodes.tsv --sample NA12878 --chr chrX  --verbose --as_quality 0.01 --shortest_path >> output.txt

python ./result_statistic.py
