from utils.file_utils import Read, Variant, get_base_pair_by_var_pos, get_geno_by_allele
import collections
from rich import print
import numpy as np
from bidict import bidict
from pympler import asizeof

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

def read_var_map(reads, variants: list[Variant]):
    '''
    Map variants onto the reads.
    First all the variants are sorted by their positions (ascending order).
    Then utilizes binary_search to find the upper bound and lower bound of a given read to all the reads.
    Make intersections and get the rest variants are located on the read.
    '''
    # First, sort all the variants by their starting point on the genome.
    variants = sorted(variants, key=lambda x: x.end)
    total_match = 0
    match_once = 0
    read_count = 0
    taken_read_count = 0
    
    allele_read_matchs = 0
    false_read_matchs = 0
    vars_ = set()
    edge_barcode_map = collections.defaultdict(dict)
    # TO deal with paired reads. Map alleles to qnames instead of reads.
    allele_linkage_map = collections.defaultdict(int)

    barcodes_bid_map = bidict()
    alleles_aid_map = bidict()

    incremental_aid = 0
    incremental_bid = 0

    # Then, map read on variants by binary search
    for read in reads:
        read_count += 1
        mapped_variants = []
        for span in read.read_spans:
            left_end = binary_search(span[0], variants, right_bound = 0)
            right_end = binary_search(span[1], variants, right_bound = 1)

            if right_end <= left_end:
                continue
            total_match += right_end - left_end
            taken_read_count += 1
            mapped_variants += variants[left_end:right_end]

        alleles = []
        for var in mapped_variants:
            # get the base pair of given read
            allele = get_base_pair_by_var_pos(read, var.end)
            # if the allele=='.', means the loci is deleted/introns or simply with low confidence.
            # we thus consider this var is not on the read.
            if allele == None:
                false_read_matchs += 1
                continue
            times = len(allele)
            allele_read_matchs += 1
            geno = get_geno_by_allele(var, allele[0])
            alleles.append(var.unique_id+':'+str(geno)+'*'+str(times))
        # there's two ways of implementation
        # first is just link the closest pair of alleles on reads
        # second is link a allele with all the alleles on a same read.
        # the first is O(N)
        # the second is O(N^2)
        # but the second will not discard any read information
        allele_list = sorted(list(set(alleles)))
        _, barcode = read.umi_barcode.split('.')
        for i in range(0, len(allele_list)-1):
            for j in range(i+1, len(allele_list)):
                (allele_1, times_1), (allele_2, times_2) = allele_list[i].split('*'), allele_list[j].split('*')
                if allele_1 not in alleles_aid_map.keys():
                    alleles_aid_map[allele_1] = incremental_aid
                    incremental_aid += 1
                if allele_2 not in alleles_aid_map.keys():
                    alleles_aid_map[allele_2] = incremental_aid
                    incremental_aid += 1
                aid_1, aid_2 = alleles_aid_map[allele_1], allele_linkage_map[allele_2]         
                allele_linkage_map[(aid_1, aid_2)] += 1
                vars_.add(allele_list[i].split(':')[0])
                vars_.add(allele_list[j].split(':')[0])

                if barcode not in barcodes_bid_map.keys():
                    barcodes_bid_map[barcode] = incremental_bid
                    incremental_bid += 1
                bid = barcodes_bid_map[barcode]
                if barcode not in edge_barcode_map[(aid_1, aid_2)].keys():
                    edge_barcode_map[(aid_1, aid_2)][bid] = 0
                edge_barcode_map[(aid_1, aid_2)][bid] += 1

    print('Received {} reads in total, {} reads are taken, {} reads are omitted. {} matches are found.'\
          .format(read_count, taken_read_count ,read_count - taken_read_count, total_match))
    
    print('There are {} pseudo matches, among which {} are considered matched and {} are false matches, {} variants has at least 1 neighbors.'\
          .format(allele_read_matchs+false_read_matchs, allele_read_matchs, false_read_matchs, len(vars_)))
    
    print('Size of allele_linkage_map: {}, Size of edge_barcode_map: {}'.format(asizeof(allele_linkage_map), asizeof(edge_barcode_map)))
    return allele_linkage_map, edge_barcode_map, len(vars_), barcodes_bid_map, alleles_aid_map
