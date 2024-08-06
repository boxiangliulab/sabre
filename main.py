#!/usr/bin/env python
from utils import file_utils, algo_utils, graph_utils, output_utils
import argparse
from rich import print as print___
import random
import numpy as np
import resource
import statistics

def prettify_print_header(step, content, end=''):
    print___('[bold green]----- Step {}:[/bold green] {}'.format(step, content), end=end)

def none_print(step='', content='', end=''):
    pass

def print_to_file(string):
    print(string, file=open('output.txt', 'a'))

chromosome_status_dict = {}

def main(opt, status_dict, return_list = None):


    if opt.total_chr != None:
        status_dict[opt.chr] = 1
        print_ = none_print
        print__ = print_to_file
    else:
        print_ = prettify_print_header
        print__ = print___
    
    print__('''[purple]
        +---------------------------------------+
        |                                       |
        |      [bold green]:dna:FASER for scRNA-Seq:dna:[/bold green] [italic purple]v3.0[/italic purple]     |
        |                                       |
        +---------------------------------------+
        [/purple]''')
    
    print_('[italic bold green]Running with Options:[/italic bold green]', vars(opt), '\n')

    processed_vcf_file = file_utils.load_vcf(opt)
    variants, vid_var_map = file_utils.generate_variants(opt, processed_vcf_file)
    bed_file = file_utils.generate_bed_file(opt, variants)

    output_sam_paths = file_utils.load_bam(opt, bed_file)
    reads = file_utils.generate_reads(opt, output_sam_paths)
    # As long as we limit the max length of alternative base pairs into 1.
    # We only need to calculate whether the $end$ of an variant lies between a read.
    allele_linkage_map, edge_barcode_map, phasable_variants, allele_linkage_read_count_map, allele_read_count, variant_allele_map = algo_utils.read_var_map(opt, reads, variants)

    allele_linkage_graph, min_mean, min_var, min_n = graph_utils.create_graph(opt, allele_linkage_map, edge_barcode_map, allele_linkage_read_count_map, allele_read_count)

    allele_subgraphs, total_possible_pairs = graph_utils.find_connected_components(allele_linkage_graph)

    conflicted_graphs, nonconflicted_graphs = graph_utils.find_conflict_graphs(opt, allele_subgraphs, vid_var_map)

    nonconflicted_nodes, phased_vars = graph_utils.extract_nonconflicted_nodes(nonconflicted_graphs)
    
    resolved_conflicted_nodes, removed_edges = graph_utils.resolve_conflict_graphs(opt, conflicted_graphs, phased_vars)

    total_hap, correct_hap, total_predict, correct_predict, total_nodes, final_graph, predict_pairs, correct_pairs, correct_variants, genome_coverage = output_utils.report_phasing_result(opt, allele_linkage_graph, nonconflicted_nodes, resolved_conflicted_nodes, vid_var_map, variant_allele_map)
    if opt.singular:
        output_utils.report_singular_cells(opt, removed_edges, final_graph, allele_linkage_graph, vid_var_map, mean=min_mean, var=min_var, n=min_n)
    if opt.allele_linkage:
        output_utils.report_allele_linkage(opt, allele_linkage_graph)
    print__("Phasing on chromosome {} COMPLETED!".format(opt.chr))
    print__('--------------------------------------------------------------')
    print__("Overall:\n#Haplotypes:\t\t {}\n#Correct Haplotypes:\t {}\n#Total variants:\t {}\n#Phased Variants:\t {}\n#Correct Variants:\t {}\nHaplotype accuracy:\t {:.4f}%\nVariants Precision:\t {:.4f}%\nVariants Recall:\t {:.4f}%\nAverage hap length:\t {:.4f}\nGenome Coverage:\t {:.4f}".format(total_hap, correct_hap,phasable_variants, total_nodes, correct_variants, correct_hap/total_hap * 100, correct_variants/total_nodes * 100, total_nodes/phasable_variants *100, total_nodes/total_hap, sum(genome_coverage)/len(genome_coverage)))
    print__('--------------------------------------------------------------')
    print__("Pairwise Metric:\n #Phased pairs:\t\t {}\nCorrect pairs:\t\t {}\nTotal pairs:\t\t {}\nPairwise accuracy:\t {:.4f}%\nPairwise recall:\t {:.4f}%".format(predict_pairs, correct_pairs, total_possible_pairs, correct_pairs/predict_pairs* 100, predict_pairs/total_possible_pairs * 100))
    print__("Global maximum memory usage: {}MB\nTime consumed: {}s".format(round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024., 4), round((resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime), 4)))
    print__('''[purple]
        +---------------------------------------+
        |                                       |
        |         [bold green]:trophy:Mission Completed!:trophy:[/bold green]        |
        |                                       |
        +---------------------------------------+
        [/purple]''')

    if opt.total_chr != None:
        status_dict[opt.chr] = 2
    
    if return_list is not None:
        return_list[opt.chr] =(total_hap, correct_hap,phasable_variants, total_nodes, correct_variants, genome_coverage, predict_pairs, correct_pairs, total_possible_pairs)


import time
import os

def watcher(chromosome_status_dict):
    os.system('clear')
    while True:
        res = ''
        sep_count = 0
        for chr_, status in chromosome_status_dict.items():
            sep_count += 1
            res += '{} Status: [ {} ]'.format(chr_ + ' '*(6-len(chr_)), '✅' if status == 2 else ('⌛️' if status == 1 else '😴'))
            res += '\t\t' if sep_count % 4 != 0 else '\n'
        print(res)
        time.sleep(0.5)
        os.system('clear')
        if 0 not in chromosome_status_dict.values() and 1 not in chromosome_status_dict.values():
            break


if __name__ == '__main__':
    # # Processing args
    parser = argparse.ArgumentParser()

    parser.add_argument("--id", help="A unique run ID string (e.g. sample345)", default='scFaser_test_output')
    parser.add_argument("--bam", help="Indexed BAMs (comma separated) containing aligned reads", default='')
    parser.add_argument("--bam_list", help="A list of input BAM files", required = False, default=None)
    parser.add_argument("--vcf", help="VCF for the sample, must be gzipped and tabix indexed.", default='')
    parser.add_argument("--sample", help="Sample name in VCF", required = False, default='')
    parser.add_argument("--npy_path", help="If var_format is set to npy, then this argument determines the path of npy file")
    parser.add_argument("--chr", help="To restrict phasing in a given chr on BAM & VCF",default='all', type=str)
    parser.add_argument("--raw_vcf", help="If the vcf is not filtered", action='store_true')
    parser.add_argument("--var_format", help="How variants are determined", default='vcf', choices=['vcf','npy'])
    parser.add_argument("--vcf_qual", help="The quality threshold on QUAL during processing vcf files.", default=30, type=int)
    parser.add_argument("--interval_threshold", help="Alleles with interval more than this threshold will be considered disconnected.", type=int, default=5000)
    parser.add_argument("--method", help="Split method, e.g. mincut, fiedler", default='fiedler')
    parser.add_argument("--sep", help="Character used to construct split variant information", type=str, default='_')
    parser.add_argument("--non_binary_variant", help="If set, means there may be multipul ALT. for a single variant", action='store_false')
    parser.add_argument("--chr_vcf", help="To restrict phasing in a given chr on VCF, if chromosome is not named equally between BAM and VCF",default=None, type=str)
    parser.add_argument("--neglect_hla", help="Indicate whether neglect variants in HLA region", action='store_true')
    parser.add_argument("--black_list", help="A blacklist, not implemented yet",default=None, type=str)
    parser.add_argument("--mapq_threshold", '--mapq', help="A filter on bam file. Reads have mapq lower than this threshold will be omitted.",default=60, type=str)
    parser.add_argument("--fiedler_threshold", help="Nodes with corresponding value in fiedler vector lower than threshold will be removed",default=1e-2, type=float)
    parser.add_argument("--remove_node", help="Remove no more than $remove_node$ in split_graph_by_common_shortest_path",default='auto', type=str)
    parser.add_argument("--shortest_path", help="Decide whether activate split_graph_by_common_shortest_path.", action='store_true')
    parser.add_argument("--as_quality", help="A filter on alignment score in BAM files", default=0.05, type=float)
    parser.add_argument("--edge_threshold", help="A filter on low confidence edges on graph", default=10, type=int)
    parser.add_argument("--verbose", help="Determine whether output conflicted graphs", action='store_true')
    parser.add_argument("--seed", help="Random seed", type=int, default=42)
    parser.add_argument("--input_type", help="Decide the input type, e.g. cellranger, umitools", type=str, default='cellranger', choices=['cellranger', 'umitools', 'star'])
    parser.add_argument("--output_conflict", help="Decide whether to output conflict graphs", action='store_true')
    parser.add_argument("--singular", help="Decide whether perform singular cell detection", action='store_true') 
    parser.add_argument("--allele_linkage", help="Decide whether output allele linkage count", action='store_true') 
    parser.add_argument("--thread", help="Number of multithread number", type=int, default=8)
    parser.add_argument("--total_chr", help="Total chromosome count for whole genome phasing", type=int, default=None)
    parser.add_argument("--chr_prefix", help="Chromosome prefix, default chr", type=str, default='chr')
    parser.add_argument("--output_vcf", help="Decide whether output vcf or not (severe performance decrease)", action='store_true')
    parser.add_argument("--tmp_dir", help="Directory of tempfile", type=str, default='./')

    opt = parser.parse_args()
    
    opt.chr_vcf = opt.chr if opt.chr_vcf is None else opt.chr_vcf

    random.seed(opt.seed)
    np.random.seed(opt.seed)
    
    if opt.chr == 'all':
        from multiprocessing import Pool, Process, Manager
        import copy

        if os.path.exists('./output.txt'):
            os.remove('./output.txt')
        from time import gmtime, strftime
        with Pool(opt.thread) as pool, Manager() as manager, open('SABRE.metrics.output', 'a') as f:
            args = []
            chromosome_status_dict = manager.dict()
            return_list = manager.dict()
            for idx in list(range(1, opt.total_chr+1))+ ['X']:
            # for idx in [19] + list(range(1, 19)) + list(range(20,23)):
                chr_ = opt.chr_prefix+str(idx)
                chromosome_status_dict[chr_] = 0
                temp = copy.deepcopy(opt)
                temp.chr = chr_
                temp.chr_vcf = chr_
                args.append((temp, chromosome_status_dict, return_list))
            p = Process(target=watcher, args=(chromosome_status_dict,))
            p.daemon = True
            p.start()
            pool.starmap(main, args)
            total_hap, correct_hap,phasable_variants, total_nodes, correct_variants, genome_coverage, predict_pairs, correct_pairs, total_possible_pairs = 0, 0, 0, 0, 0, [], 0, 0, 0
            for chr_, list_ in return_list.items():
                total_hap += list_[0]
                correct_hap += list_[1]
                phasable_variants += list_[2]
                total_nodes += list_[3]
                correct_variants += list_[4]
                genome_coverage += list_[5]
                predict_pairs += list_[6]
                correct_pairs += list_[7]
                total_possible_pairs += list_[8]
            p.kill()
            print('--------------------------------------------------------------')
            print("Overall:\n#Haplotypes:\t\t {}\n#Correct Haplotypes:\t {}\n#Total variants:\t {}\n#Phased Variants:\t {}\n#Correct Variants:\t {}\nHaplotype accuracy:\t {:.4f}%\nVariants Precision:\t {:.4f}%\nVariants Recall:\t {:.4f}%\nAverage hap length:\t {:.4f}\nGenome Coverage Median:\t {:.4f}".format(total_hap, correct_hap,phasable_variants, total_nodes, correct_variants, correct_hap/total_hap * 100, correct_variants/total_nodes * 100, total_nodes/phasable_variants *100, total_nodes/total_hap, statistics.median(genome_coverage))) 
            print('--------------------------------------------------------------')
            print("Pairwise Metric:\n#Phased pairs:\t\t {}\nCorrect pairs:\t\t {}\nTotal pairs:\t\t {}\nPairwise accuracy:\t {:.4f}%\nPairwise recall:\t {:.4f}%".format(predict_pairs, correct_pairs, total_possible_pairs, correct_pairs/predict_pairs* 100, correct_pairs/total_possible_pairs * 100)) 
            
            print('--------------------------------------------------------------', file=f)
            print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()), file=f)
            print(opt, file=f)
            print('--------------------------------------------------------------', file=f)
            print("Overall:\n#Haplotypes:\t\t {}\n#Correct Haplotypes:\t {}\n#Total variants:\t {}\n#Phased Variants:\t {}\n#Correct Variants:\t {}\nHaplotype accuracy:\t {:.4f}%\nVariants Precision:\t {:.4f}%\nVariants Recall:\t {:.4f}%\nAverage hap length:\t {:.4f}\nGenome Coverage Median:\t {:.4f}".format(total_hap, correct_hap,phasable_variants, total_nodes, correct_variants, correct_hap/total_hap * 100, correct_variants/total_nodes * 100, total_nodes/phasable_variants *100, total_nodes/total_hap, statistics.median(genome_coverage)), file=f) 
            print('--------------------------------------------------------------', file=f)
            print("Pairwise Metric:\nPhased pairs:\t\t {}\nCorrect pairs:\t\t {}\nTotal pairs:\t\t {}\nPairwise accuracy:\t {:.4f}%\nPairwise recall:\t {:.4f}%".format(predict_pairs, correct_pairs, total_possible_pairs, correct_pairs/predict_pairs* 100, correct_pairs/total_possible_pairs * 100), file=f)
    else:
        main(opt, None)


