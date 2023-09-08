#!/usr/bin/env python
from utils import file_utils, algo_utils, graph_utils, output_utils
import argparse
from rich import print
import random

def prettify_print_header(step, content, end=''):
    print('[bold green]----- Step {}:[/bold green] {}'.format(step, content), end=end)
    

def main(opt):

    random.seed(114514)
    print('''[purple]
        +---------------------------------------+
        |                                       |
        |      [bold green]:dna:FASER for scRNA-Seq:dna:[/bold green] [italic purple]v0.1[/italic purple]     |
        |                                       |
        +---------------------------------------+
        [/purple]''')
    
    print('[italic bold green]Running with Options:[/italic bold green]', vars(opt), '\n')

    prettify_print_header(1, 'Loading and Preprocessing VCF File...', end='\r')
    processed_vcf_file = file_utils.load_vcf(opt)
    variants, vid_var_map = file_utils.generate_variants(processed_vcf_file)
    bed_file = file_utils.generate_bed_file(variants)
    prettify_print_header(1, 'Loading and Preprocessing VCF File [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(2, 'Loading and Preprocessing BAM File...', end='\r')
    output_sam_path = file_utils.load_bam(opt, bed_file)
    reads = file_utils.generate_reads(opt, output_sam_path)
    prettify_print_header(2, 'BAM loading and preprocessing [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    # As long as we limit the max length of alternative base pairs into 1.
    # We only need to calculate whether the $end$ of an variant lies between a read.
    prettify_print_header(3, 'Mapping Variants to Reads...', end='\r')
    read_variants_map = algo_utils.read_var_map(reads, variants, vid_var_map)
    prettify_print_header(3, 'Mapping Variants to Reads [pink1 bold]COMPLETED![/pink1 bold]!', '\n\n')

    prettify_print_header(4, 'Mapping Alleles to Reads...', end='\r')
    allele_linkage_map = algo_utils.extract_allele_linkage(read_variants_map)
    prettify_print_header(4, 'Mapping Alleles to Reads [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(5, 'Creating the allele linkage graph...', end='\r')
    allele_linkage_graph = graph_utils.create_graph(opt, allele_linkage_map)
    prettify_print_header(5, 'Creating the allele linkage graph [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(6, 'Finding connected components and save them...', end='\r')
    allele_subgraphs = graph_utils.find_connected_components(allele_linkage_graph)
    prettify_print_header(6, 'Finding connected components and save them [pink1 bold]COMPLETED![/pink1 bold]', end='\n\n')

    prettify_print_header(7, 'Finding conflicted subgraphs...', end='\r')
    conflicted_graphs, nonconflicted_graphs = graph_utils.find_conflict_graphs(opt, allele_subgraphs, vid_var_map)
    prettify_print_header(7, 'Finding conflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(8, 'Reporting nonconflicted subgraphs...', end='\r')
    nonconflicted_nodes = graph_utils.extract_nonconflicted_nodes(nonconflicted_graphs)
    prettify_print_header(8, 'Reporting nonconflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(9, 'Resolving conflicted subgraphs...', end='\r')
    resolved_conflicted_nodes = graph_utils.resolve_conflict_graphs(opt, conflicted_graphs)
    prettify_print_header(9, 'Resolving conflicted subgraphs [pink1 bold]COMPLETED![/pink1 bold]', '\n\n')

    prettify_print_header(10, 'Reporting phasing result...', end='\r')

    total_hyp, correct_hyp, total_predict, correct_predict, total_nodes = output_utils.report_phasing_result(opt, allele_linkage_graph, nonconflicted_nodes, resolved_conflicted_nodes, vid_var_map)
    prettify_print_header(10, 'Reporting phasing result [pink1 bold]COMPLETED![/pink1 bold]', '\n')
    print("Phasing on chromosome {} [pink1 bold]COMPLETED![/pink1 bold]".format(opt.restrict_chr))
    print("[green bold]Phased Vars:\t[/green bold] {} variants in total.".format(total_nodes))
    print("[green bold]Overall:\t[/green bold] {} hyplotypes in total, {} hyplotypes are in coordinate with ground truth. The accuracy is {:.4f}%".format(total_hyp, correct_hyp, correct_hyp/total_hyp * 100))
    print("[green bold]Conflicted:\t[/green bold] {}  hyplotypes in total, {}  hyplotypes are in coordinate with ground truth. The accuracy is {:.4f}%".format(total_predict, correct_predict, 0 if total_predict == 0 else correct_predict/total_predict * 100))
    print("[green bold]Nonconflicted:\t[/green bold] {} hyplotypes in total, {} hyplotypes are in coordinate with ground truth. The accuracy is {:.4f}%".format(total_hyp-total_predict, correct_hyp-correct_predict, (correct_hyp-correct_predict)/(total_hyp-total_predict) * 100))

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

    parser.add_argument("--bam_path", help="Indexed BAMs (comma separated) containing aligned reads", required = True, default='')
    parser.add_argument("--vcf_path", help="VCF for the sample, must be gzipped and tabix indexed.", required = True, default='')
    parser.add_argument("--sample_name", help="Sample name in VCF", required = False, default='')
    parser.add_argument("--restrict_chr", help="To restrict phasing in a given chr",default=None, type=str)
    parser.add_argument("--black_list", help="A blacklist, not implemented yet",default=None, type=str)
    parser.add_argument("--mapq_threshold", help="A filter on bam file. Reads have mapq lower than this threshold will be omitted.",default=60, type=str)
    parser.add_argument("--as_quality", help="A filter on alignment score in BAM files", default=0.05, type=float)
    parser.add_argument("--edge_threshold", help="A filter on low confidence edges on graph", default=5, type=int)
    parser.add_argument("--verbose", help="Determine whether output conflicted graphs", action='store_true')
    

    opt = parser.parse_args()
    
    
    main(opt)