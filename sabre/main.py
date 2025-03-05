#!/usr/bin/env python
from sabre.utils import file_utils, algo_utils, graph_utils, output_utils
import argparse
from rich import print as print___
import random
import numpy as np
import resource
import statistics

import warnings
warnings.filterwarnings("ignore")

def prettify_print_header(step, content, end=''):
    print___('[bold green]----- Step {}:[/bold green] {}'.format(step, content), end=end)

def none_print(step='', content='', end=''):
    pass

def print_to_file(string):
    print(string, file=open('output.txt', 'a'))

chromosome_status_dict = {}

def sabre(opt, status_dict, return_list = None, status_list=None):


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
        |               [bold italic green]:dna:sabre:dna:[/bold italic green]               |
        |                                       |
        +---------------------------------------+
        [/purple]''')
    
    print_('[italic bold green]Running with Options:[/italic bold green]', vars(opt), '\n')

    try:
        processed_vcf_file = file_utils.load_vcf(opt)
        variants, somatic_variants, vid_var_map = file_utils.generate_variants(opt, processed_vcf_file)
        bed_file = file_utils.generate_bed_file(opt, variants, somatic_variants)

        output_sam_paths = file_utils.load_bam(opt, bed_file)
        reads = file_utils.generate_reads(opt, output_sam_paths)
        # As long as we limit the max length of alternative base pairs into 1.
        # We only need to calculate whether the $end$ of an variant lies between a read.
        allele_linkage_map, edge_barcode_map, phasable_variants, allele_linkage_read_count_map, allele_read_count, variant_allele_map = algo_utils.read_var_map(opt, reads, variants, somatic_variants)

        allele_linkage_graph, min_mean, min_var, min_n = graph_utils.create_graph(opt, allele_linkage_map, edge_barcode_map, allele_linkage_read_count_map, allele_read_count)

        allele_subgraphs, total_possible_pairs = graph_utils.find_connected_components(allele_linkage_graph)

        conflicted_graphs, nonconflicted_graphs = graph_utils.find_conflict_graphs(opt, allele_subgraphs, vid_var_map)

        nonconflicted_nodes, phased_vars = graph_utils.extract_nonconflicted_nodes(nonconflicted_graphs)
        
        resolved_conflicted_nodes, removed_edges = graph_utils.resolve_conflict_graphs(opt, conflicted_graphs, phased_vars)

        total_hap, correct_hap, total_predict, correct_predict, total_nodes, final_graph, predict_pairs, correct_pairs, correct_variants, genome_coverage = output_utils.report_phasing_result(opt, allele_linkage_graph, nonconflicted_nodes, resolved_conflicted_nodes, vid_var_map, variant_allele_map)
        if opt.residual_edges:
            output_utils.report_singular_cells(opt, removed_edges, final_graph, allele_linkage_graph, vid_var_map, variant_allele_map, mean=min_mean, var=min_var, n=min_n)
        if opt.allele_linkage:
            output_utils.report_allele_linkage(opt, allele_linkage_graph)
        print__("Phasing on chromosome {} COMPLETED!".format(opt.chr))

        if opt.benchmark:
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

    except Exception as e:
        if opt.total_chr == None:
            raise e
        
        status_dict[opt.chr] = -1
        status_list.append('Exception on {}: {}'.format(opt.chr, e))


import time
import os

def watcher(chromosome_status_dict, status_list):
    os.system('clear')
    while True:
        res = ''
        sep_count = 0
        for chr_, status in chromosome_status_dict.items():
            sep_count += 1
            res += '{} Status: [ {} ]'.format(chr_ + ' '*(6-len(chr_)), '✅' if status == 2 else ('⌛️' if status == 1 else ( '❌' if status == -1  else '😴')))
            res += '\t\t' if sep_count % 4 != 0 else '\n'
        print(res)
        list(map(lambda x: print('[bold red blink]❌ ERROR: [/bold red blink] {}'.format(x)), status_list))
        if 0 not in chromosome_status_dict.values() and 1 not in chromosome_status_dict.values():
            break
        time.sleep(0.5)
        os.system('clear')



def main():
    # # Processing args
    parser = argparse.ArgumentParser()

    # General Options
    parser.add_argument("--id", help="A unique run ID string (e.g. sample345)", default='sabre_test', required=True)
    parser.add_argument("--bam", help="Indexed BAMs (comma separated) containing aligned reads", default='')
    parser.add_argument("--bam_list", help="A list of input BAM files", required = False, default=None)
    parser.add_argument("--vcf", help="VCF for the sample, must be gzipped and tabix indexed.", default='')
    parser.add_argument("--sample", help="Sample name in VCF", required = True, default='')
    parser.add_argument("--chr", help="To restrict phasing in a given chr on BAM & VCF",default='all', type=str)
    parser.add_argument("--total_chr", help="Total chromosome count for whole genome phasing", type=int, default=None)
    parser.add_argument("--raw_vcf", help="If the vcf is not filtered", action='store_true')
    parser.add_argument("--sep", help="Character used to construct split variant information", type=str, default='_')
    parser.add_argument("--non_binary_variant", help="If set, means there may be multipul ALT. for a single variant", action='store_true')
    parser.add_argument("--chr_prefix", help="Chromosome prefix, default chr", type=str, default='chr')
    parser.add_argument("--chr_vcf", help="To restrict phasing in a given chr on VCF, if chromosome is not named equally between BAM and VCF",default=None, type=str)
    parser.add_argument("--neglect_hla", help="Indicate whether neglect variants in HLA region", action='store_true')
    parser.add_argument("--seed", help="Random seed", type=int, default=42)
    parser.add_argument("--tmp_dir", help="Directory of tempfile", type=str, default='./')
    parser.add_argument("--output_dir", help="The path to output directory. Default ./output", default='./output', type=str)
    parser.add_argument("--input_type", help="How umi-barcode is provided, e.g. cellranger-style, umitools-style or 're' for custom regular expression.", type=str, default='cellranger', choices=['cellranger', 'umitools', 'star', 're', 'bulk', 'smartseq'])
    parser.add_argument("--bc_re", help="The regular expression for extracting Cell Barcode in the BAM file.", type=str, default=None)
    parser.add_argument("--umi_re", help="The regular expression for extracting UMI in the BAM file.", type=str, default=None)

    # Accuracy & Sensitivity Related Options
    parser.add_argument("--mapq_threshold", '--mapq', help="A filter on bam file. Reads have mapq lower than this threshold will be omitted.",default=60, type=int)
    parser.add_argument("--vcf_qual", help="The quality threshold on QUAL during processing vcf files.", default=30, type=int)
    parser.add_argument("--interval_threshold", help="Alleles with interval more than this threshold will be considered disconnected.", type=int, default=5000)
    parser.add_argument("--base_conflict_threshold", help="Base pairs that diffs less than this threshold will be ignored.", type=float, default=0.05)
    parser.add_argument("--method", help="Split method, e.g. mincut, fiedler", default='fiedler', choices=['fiedler', 'mincut'])
    parser.add_argument("--fiedler_threshold", help="Nodes with corresponding value in fiedler vector lower than threshold will be removed",default=1e-2, type=float)
    parser.add_argument("--shortest_path", help="Decide whether activate split_graph_by_common_shortest_path.", action='store_true')
    parser.add_argument("--remove_node", help="Remove no more than $remove_node$ in split_graph_by_common_shortest_path",default='auto', type=str)
    parser.add_argument("--as_quality", help="A filter on alignment score in BAM files", default=0.05, type=float)
    parser.add_argument("--edge_threshold", help="A filter on low confidence edges on graph", default=10, type=int)
    parser.add_argument("--layers", help="Number of GNN Layer", type=int, default=1)

    # Output Options
    parser.add_argument("--verbose", help="Determine whether output conflicted graphs", action='store_true')
    parser.add_argument("--benchmark", help="If set true, will output benchmark metrics", action='store_true')
    parser.add_argument("--output_conflict", help="Decide whether to output conflict graphs", action='store_true')
    parser.add_argument("--residual_edges", help="Decide whether report residual edges.", action='store_true') 
    parser.add_argument("--allele_linkage", help="Decide whether output allele linkage count", action='store_true') 
    parser.add_argument("--output_vcf", help="Decide whether output vcf or not (severe performance decrease)", action='store_true')
    parser.add_argument("--no_vcf_id", help="If set, the 'ID' row in vcf file will not be retained", action='store_true')

    # Performance Related Options
    parser.add_argument("--thread", help="Number of multithread number", type=int, default=8)

    # Monopogen somatic related arguments
    parser.add_argument("--mono", help="The result putativeSNVs.csv of Monopogen somatic output, used for somatic variation analysis.", default=None)
    parser.add_argument("--mono_svm", help="Threshold on SVM_pos_score in putativeSNVs.csv.", default=0.5, type=float)
    parser.add_argument("--mono_ld", help="Threshold on LDrefine_merged_score in putativeSNVs.csv.", default=0.5, type=float)
    parser.add_argument("--mono_ignore_ld", help="Threshold on LDrefine_merged_score in putativeSNVs.csv.", default=0.5, type=float)
    parser.add_argument("--mono_baf_u", help="The upper bound on BAF_alt in putativeSNVs.csv.", default=0.5, type=float)
    parser.add_argument("--mono_baf_l", help="The lower bound on BAF_alt in putativeSNVs.csv.", default=0.1, type=float)
    parser.add_argument("--mono_ref", help="Threshold on Dep_ref in putativeSNVs.csv.", default=5, type=int)
    parser.add_argument("--mono_alt", help="Threshold on Dep_alt in putativeSNVs.csv.", default=5, type=int)

    opt = parser.parse_args()
    
    ## argument checking
    opt.chr_vcf = opt.chr if opt.chr_vcf is None else opt.chr_vcf
    if opt.input_type == 're':
        assert opt.bc_re is not None and opt.umi_re is not None
        try:
            import re
            re.compile(opt.bc_re)
            re.compile(opt.umi_re)
        except Exception as e:
            print___('[bold red blink]❌ ERROR: [/bold red blink]Regular Expression compiling error. Please check the vadility of input re.')
            raise e

    random.seed(opt.seed)
    np.random.seed(opt.seed)

    if not os.path.exists(opt.output_dir):
        os.mkdir(opt.output_dir)
    
    if not os.path.exists('{}/{}'.format(opt.output_dir, opt.id)):
        os.mkdir('{}/{}'.format(opt.output_dir, opt.id))
    
    if opt.vcf_qual < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--vcf_qual must be set >= 0; Current value: {}'.format(opt.vcf_qual))
        raise ValueError("--vcf_qual less than 0")

    if opt.interval_threshold < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--interval_threshold must be set >= 0; Current value: {}'.format(opt.interval_threshold))
        raise ValueError("--interval_threshold less than 0")

    if opt.mapq_threshold < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--mapq_threshold must be set >= 0; Current value: {}'.format(opt.mapq_threshold))
        raise ValueError("--mapq_threshold less than 0")

    if opt.fiedler_threshold < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--fiedler_threshold must be set >= 0; Current value: {}'.format(opt.fiedler_threshold))
        raise ValueError("--fiedler_threshold less than 0")

    if opt.as_quality < 0 or opt.as_quality >1:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--as_quality must be set >= 0 and <= 1; Current value: {}'.format(opt.as_quality))
        raise ValueError("--as_quality less than 0 or greater than 1")

    if opt.edge_threshold > 100:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--edge_threshold must be set >= 0; Current value: {}'.format(opt.edge_threshold))
        raise ValueError("--edge_threshold greater than 100")

    if opt.thread < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--thread must be set >= 0; Current value: {}'.format(opt.thread))
        raise ValueError("--thread less than 0")
    
    if opt.layers < 0:
        print___('[bold red blink]❌ ERROR: [/bold red blink]--layers must be set >= 0; Current value: {}'.format(opt.layers))
        raise ValueError("--layers less than 0")
    
    if opt.chr == 'all' and opt.total_chr is None:
        print___('[bold red blink]❌ ERROR: [/bold red blink]No chromosome is specified. Please check arguments.')
        raise ValueError("No chromosome specified")

    # Perform pre-flight check

    def perform_bam_check(bam):
        if not os.path.exists(bam):
            print___('[bold red blink]❌ ERROR: [/bold red blink]BAM file does not exist. Please check your input arguments.')
            raise RuntimeError('BAM file not exists')
        
        if not os.path.exists('{}.bai'.format(bam)):
            print___('[bold red blink]❌ ERROR: [/bold red blink]BAM file {} not indexed. Sabre requires each input BAM file indexed by samtools!'.format(bam))
            raise KeyError('BAM file not indexed')

    if not opt.bam == '':
        perform_bam_check(opt.bam)
  
    if opt.bam == '':
        if opt.bam_list == None:
            raise ValueError("No bam file specified. Please check your input arguments")
        if not os.path.exists(opt.bam_list):
            raise ValueError("Invalid bam_list file. Please check your input arguments")
        with open(opt.bam_list) as f:
            for line in f:
                bamfile = line.split(',')[0]
                try:
                    perform_bam_check(bamfile)
                except RuntimeError as e:
                    print___('[bold red blink]❌ ERROR: [/bold red blink]BAM file {} does not exist. Please check input arguments'.format(bamfile))
                    raise e
                except KeyError as e:
                    print___('[bold red blink]❌ ERROR: [/bold red blink]BAM file {} not indexed. Sabre needs each input BAM file indexed by samtools!'.format(bamfile))
                    raise e
    
    if not os.path.exists(opt.vcf):
        print___('[bold red blink]❌ ERROR: [/bold red blink]VCF file not exists. Please check input arguments.')
        raise RuntimeError('VCF file not exists!')
    
    if not os.path.exists('{}.tbi'.format(opt.vcf)):
        print___('[bold red blink]❌ ERROR: [/bold red blink]VCF file {} not indexed. Sabre requires input VCF file index by tabix'.format(opt.vcf))
        raise RuntimeError('VCF file not indexed!')
    
    if opt.mono is not None and not os.path.exists(opt.mono):
        print___('[bold red blink]❌ ERROR: [/bold red blink]Monopogen output file not exists. Please check input arguments.')
        raise RuntimeError('Monopogen file not exists!')

    if opt.chr == 'all':
        from multiprocessing import Pool, Process, Manager
        import copy

        if os.path.exists('./SABRE.metrics.output'):
            os.remove('./SABRE.metrics.output')
        from time import gmtime, strftime
        with Pool(opt.thread) as pool, Manager() as manager, open('SABRE.metrics.output', 'a') as f:
            args = []
            chromosome_status_dict = manager.dict()
            status_list = manager.list()
            return_list = manager.dict()
            
            for idx in list(range(1, opt.total_chr+1))+ ['X']:
            # for idx in [19] + list(range(1, 19)) + list(range(20,23)):
                chr_ = opt.chr_prefix+str(idx)
                chromosome_status_dict[chr_] = 0
                temp = copy.deepcopy(opt)
                temp.chr = chr_
                temp.chr_vcf = chr_
                args.append((temp, chromosome_status_dict, return_list, status_list))
                
            p = Process(target=watcher, args=(chromosome_status_dict, status_list))
            p.daemon = True
            p.start()
            pool.starmap(sabre, args)
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
        
        if opt.benchmark:
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
        sabre(opt, None)


if __name__ == '__main__':
    main()