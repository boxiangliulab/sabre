from utils.file_utils import Read, Variant
import collections
from rich import print
import numpy as np

def binary_search(pos:int, variants:list[Variant], right_bound=0):
    '''
    Typical binary search,
    Time complexity: O(logn)
    '''

    left, right = 0, len(variants) - 1
    while left <= right:
        mid = left + (right - left) // 2
        value = variants[mid].end
        if value == pos:
            return mid + right_bound
        elif value < pos:
            left = mid + 1
        else:
            right = mid - 1
    return left

def read_var_map(reads: list[Read], variants: list[Variant], vid_var_map):
    '''
    Map variants onto the reads.
    First all the variants are sorted by their positions (ascending order).
    Then utilizes binary_search to find the upper bound and lower bound of a given read to all the reads.
    Make intersections and get the rest variants are located on the read.
    '''
    # First, sort all the variants by their starting point on the genome.
    variants = sorted(variants, key=lambda x: x.end)
    read_variants_map = collections.defaultdict(list)
    total_match = 0
    match_once = 0

    # Then, map read on variants by binary search
    for read in reads:
        for span in read.read_spans:
            left_end = binary_search(span[0], variants, right_bound = 0)
            right_end = binary_search(span[1], variants, right_bound = 1)

            if right_end <= left_end:
                continue
            total_match += right_end - left_end
            read_variants_map[read] += variants[left_end:right_end]
    
    print('Received {} reads in total, {} reads are taken, {} reads are omitted. {} matches are found.'\
          .format(len(reads), len(read_variants_map),len(reads) - len(read_variants_map), total_match))

    return read_variants_map

def extract_allele_linkage(opt, read_variants_map: dict):
    '''
    As we need to construct a allele graph, we need to figure out the mapping relationship between alleles and reads.
    For scRNA data, variants on reads with the same umi-barcode are actually on the same haplotype, thus shall be considered connected.
    '''
    allele_read_matchs = 0
    false_read_matchs = 0
    vars_ = set()
    edge_barcode_map = collections.defaultdict(dict)
    # TO deal with paired reads. Map alleles to qnames instead of reads.
    qname_alleles_map = collections.defaultdict(list) 
    for read, variants in read_variants_map.items():
        for var in variants:
            # get the base pair of given read
            allele = read.get_base_pair_by_var_pos(var.end)
            # if the allele=='.', means the loci is deleted/introns or simply with low confidence.
            # we thus consider this var is not on the read.
            if allele == None:
                false_read_matchs += 1
                continue
            times = len(allele)
            allele_read_matchs += 1
            geno = var.get_geno_by_allele(allele[0])
            qname_alleles_map[read.umi_barcode].append(var.unique_id+':'+str(geno)+'*'+str(times))
    # with open('bu_var_map_count.txt', 'a') as f:
    #    for barcode, allele_list in qname_alleles_map.items():
    #        f.write('{}\n'.format(len(allele_list)))
    # exit(0)
    # there's two ways of implementation
    # first is just link the closest pair of alleles on reads
    # second is link a allele with all the alleles on a same read.
    # the first is O(N)
    # the second is O(N^2)
    # but the second will not discard any read information
    allele_linkage_map = collections.defaultdict(int)
    for qname, allele_list in qname_alleles_map.items():
        allele_list = sorted(list(set(allele_list)))
        _, barcode = qname.split('.')
        for i in range(0, len(allele_list)-1):
            for j in range(i+1, len(allele_list)):
                (allele_1, times_1), (allele_2, times_2) = allele_list[i].split('*'), allele_list[j].split('*')
                allele_linkage_map[(allele_1, allele_2)] += 1
                vars_.add(allele_list[i].split(':')[0])
                vars_.add(allele_list[j].split(':')[0])
                if barcode not in edge_barcode_map[(allele_1, allele_2)].keys():
                    edge_barcode_map[(allele_1, allele_2)][barcode] = 0
                edge_barcode_map[(allele_1, allele_2)][barcode] += 1
    print('There are {} pseudo matches, among which {} are considered matched and {} are false matches, {} variants has at least 1 neighbors.'\
          .format(allele_read_matchs+false_read_matchs, allele_read_matchs, false_read_matchs, len(vars_)))
    # with open('./metric/{}_neighbored_vars.txt'.format(opt.restrict_chr), 'w') as f:
        # list(map(lambda x: f.write('{}\n'.format(x)), vars_))
    return allele_linkage_map, edge_barcode_map
