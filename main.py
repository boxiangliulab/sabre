#!/usr/bin/env python
from utils import file_utils, algo_utils, graph_utils, output_utils
import argparse
from rich import print
import random
import numpy as np
import resource

def prettify_print_header(step, content, end=''):
    print('[bold green]----- Step {}:[/bold green] {}'.format(step, content), end=end)
    

def main(opt):

    random.seed(opt.seed)
    np.random.seed(opt.seed)
    
    print('''[purple]
        +---------------------------------------+
        |                                       |
        |      [bold green]:dna:FASER for scRNA-Seq:dna:[/bold green] [italic purple]v2.0[/italic purple]     |
        |                                       |
        +---------------------------------------+
        [/purple]''')
    
    print('[italic bold green]Running with Options:[/italic bold green]', vars(opt), '\n')

    prettify_print_header(1, 'Loading and Preprocessing VCF File...', end='\r')
    if opt.var_format != 'npy':
        processed_vcf_file = file_utils.load_vcf(opt)
    else:
        processed_vcf_file = None
    variants, vid_var_map = file_utils.generate_variants(opt, processed_vcf_file)
    bed_file = file_utils.generate_bed_file(opt, variants)
    prettify_print_header(1, 'Loading and Preprocessing VCF File [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(2, 'Loading and Preprocessing BAM File...', end='\r')
    output_sam_path = file_utils.load_bam(opt, bed_file)
    reads = file_utils.generate_reads(opt, output_sam_path)
    prettify_print_header(2, 'BAM loading and preprocessing [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    # As long as we limit the max length of alternative base pairs into 1.
    # We only need to calculate whether the $end$ of an variant lies between a read.
    prettify_print_header(3, 'Mapping Variants to Reads...', end='\r')
    read_variants_map = algo_utils.read_var_map(reads, variants, vid_var_map)
    del(reads)
    del(variants)
    prettify_print_header(3, 'Mapping Variants to Reads [pink1 bold]COMPLETED![/pink1 bold]!', '\n\n')

    prettify_print_header(4, 'Mapping Alleles to Reads...', end='\r')
    allele_linkage_map, edge_barcode_map = algo_utils.extract_allele_linkage(opt, read_variants_map)
    del(read_variants_map)
    prettify_print_header(4, 'Mapping Alleles to Reads [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(5, 'Creating the allele linkage graph...', end='\r')
    allele_linkage_graph, min_mean, min_var, min_n = graph_utils.create_graph(opt, allele_linkage_map, edge_barcode_map)
    prettify_print_header(5, 'Creating the allele linkage graph [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')
    
    prettify_print_header(6, 'Finding connected components and save them...', end='\r')
    allele_subgraphs = graph_utils.find_connected_components(allele_linkage_graph)
    prettify_print_header(6, 'Finding connected components and save them [pink1 bold]COMPLETED![/pink1 bold]', end='\n\n')

    prettify_print_header(7, 'Finding conflicted subgraphs...', end='\r')
    conflicted_graphs, nonconflicted_graphs = graph_utils.find_conflict_graphs(opt, allele_subgraphs, vid_var_map)
    prettify_print_header(7, 'Finding conflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(8, 'Reporting nonconflicted subgraphs...', end='\r')
    nonconflicted_nodes, phased_vars = graph_utils.extract_nonconflicted_nodes(nonconflicted_graphs)
    prettify_print_header(8, 'Reporting nonconflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(9, 'Resolving conflicted subgraphs...', end='\r')
    resolved_conflicted_nodes, removed_edges = graph_utils.resolve_conflict_graphs(opt, conflicted_graphs, phased_vars)
    prettify_print_header(9, 'Resolving conflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(10, 'Reporting phasing result...', end='\r')
    total_hap, correct_hap, total_predict, correct_predict, total_nodes, final_graph, predict_pairs, correct_pairs = output_utils.report_phasing_result(opt, allele_linkage_graph, nonconflicted_nodes, resolved_conflicted_nodes, vid_var_map)
    #total_hap, correct_hap, total_predict, correct_predict, total_nodes, final_graph = output_utils.report_phasing_result(opt, allele_linkage_graph, nonconflicted_nodes, [], vid_var_map)
    if opt.singular:
        output_utils.report_singular_cells(opt, removed_edges, final_graph, allele_linkage_graph, vid_var_map, mean=min_mean, var=min_var, n=min_n)
    prettify_print_header(10, 'Reporting phasing result [pink1 bold]COMPLETED![/pink1 bold]', '\n')
    print("Phasing on chromosome {} [pink1 bold]COMPLETED![/pink1 bold]".format(opt.restrict_chr))
    print("[green bold]Phased Vars:\t[/green bold] {} variants in total.".format(total_nodes))    
    print("[green bold]Overall:\t[/green bold] {} haplotypes in total, {} haplotypes are in coordinate with ground truth. Haplotype accuracy is {:.4f}%".format(total_hap, correct_hap, correct_hap/total_hap * 100))
    print("[green bold]Pairwise Metric:\t[/green bold] {} pairs in total, {} pairs are in coordinate with ground truth. Pairwise accuracy is {:.4f}%".format(predict_pairs, correct_pairs, correct_pairs/predict_pairs* 100))
    # print("[green bold]Conflicted:\t[/green bold] {}  haplotypes in total, {}  haplotypes are in coordinate with ground truth. The accuracy is {:.4f}%".format(total_predict, correct_predict, 0 if total_predict == 0 else correct_predict/total_predict * 100))
    # print("[green bold]Nonconflicted:\t[/green bold] {} haplotypes in total, {} haplotypes are in coordinate with ground truth. The accuracy is {:.4f}%".format(total_hap-total_predict, correct_hap-correct_predict, (correct_hap-correct_predict)/(total_hap-total_predict) * 100))
    print("Global maximum memory usage: [pink1 bold]{}[/pink1 bold]MB\nTime consumed: [pink1 bold]{}[/pink1 bold]s".format(round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024., 4), round((resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime), 4)))
    print('''[purple]
        +---------------------------------------+
        |                                       |
        |         [bold green]:trophy:Mission Completed!:trophy:[/bold green]        |
        |                                       |
        +---------------------------------------+
        [/purple]''')

if __name__ == '__main__':
    # # Processing args
    parser = argparse.ArgumentParser()

    parser.add_argument("--id", help="A unique run ID string (e.g. sample345)", required = True, default='scFaser_test_output')
    parser.add_argument("--bam_path", help="Indexed BAMs (comma separated) containing aligned reads", required = True, default='')
    parser.add_argument("--vcf_path", help="VCF for the sample, must be gzipped and tabix indexed.", default='')
    parser.add_argument("--sample_name", help="Sample name in VCF", required = False, default='')
    parser.add_argument("--npy_path", help="If var_format is set to npy, then this argument determines the path of npy file")
    parser.add_argument("--restrict_chr", help="To restrict phasing in a given chr on BAM & VCF",default=None, type=str)
    parser.add_argument("--raw_vcf", help="If the vcf is not filtered", action='store_true')
    parser.add_argument("--var_format", help="How variants are determined", default='vcf', choices=['vcf','npy'])
    parser.add_argument("--vcf_qual", help="The quality threshold on QUAL during processing vcf files.", default=30, type=int)
    parser.add_argument("--restrict_chr_vcf", help="To restrict phasing in a given chr on VCF, if chromosome is not named equally between BAM and VCF",default=None, type=str)
    parser.add_argument("--neglect_hla", help="Indicate whether neglect variants in HLA region",default=True)
    parser.add_argument("--black_list", help="A blacklist, not implemented yet",default=None, type=str)
    parser.add_argument("--mapq_threshold", help="A filter on bam file. Reads have mapq lower than this threshold will be omitted.",default=60, type=str)
    parser.add_argument("--fiedler_threshold", help="Nodes with corresponding value in fiedler vector lower than threshold will be removed",default=1e-2, type=float)
    parser.add_argument("--remove_node", help="Remove no more than $remove_node$ in split_graph_by_common_shortest_path",default='auto', type=str)
    parser.add_argument("--shortest_path", help="Decide whether activate split_graph_by_common_shortest_path.", action='store_true')
    parser.add_argument("--as_quality", help="A filter on alignment score in BAM files", default=0.05, type=float)
    parser.add_argument("--neglect_overlap", help="Neglect overlap when deal with reads overlaps", action='store_true')
    parser.add_argument("--edge_threshold", help="A filter on low confidence edges on graph", default=10, type=int)
    parser.add_argument("--verbose", help="Determine whether output conflicted graphs", action='store_true')
    parser.add_argument("--seed", help="Random seed", type=int, default=42)
    parser.add_argument("--input_type", help="Decide the input type, e.g. cellranger, umitools", type=str, default='cellranger', choices=['cellranger', 'umitools', 'star'])
    parser.add_argument("--singular", help="Decide whether perform singular cell detection", type=bool, default=True) 

    opt = parser.parse_args()
    
    opt.restrict_chr_vcf = opt.restrict_chr if opt.restrict_chr_vcf is None else opt.restrict_chr_vcf
    
    main(opt)

