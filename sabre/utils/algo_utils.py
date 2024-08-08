import sabre.utils.file_utils as file_utils
import collections
from rich import print

def binary_search(pos:int, variants, left, right, right_bound=0):
    '''
    Typical binary search,
    Time complexity: O(logn)
    '''
    while left <= right:
        mid = (left + right) >> 1
        value = variants[mid].end
        if value == pos:
            return mid + right_bound
        elif value < pos:
            left = mid + 1
        else:
            right = mid - 1
    return left

def read_var_map(opt, reads, variants):
    '''
    Map variants onto the reads.
    First all the variants are sorted by their positions (ascending order).
    Then utilizes binary_search to find the upper bound and lower bound of a given read to all the reads.
    Make intersections and get the rest variants are located on the read.
    '''
    # First, sort all the variants by their starting point on the genome.
    variants = sorted(variants, key=lambda x: x.end)
    variants_len = len(variants) - 1
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
    allele_linkage_read_count_map = collections.defaultdict(int)
    allele_read_count = collections.defaultdict(int)
    mapped_variants_cache = collections.defaultdict(int)

    variant_allele_map = collections.defaultdict(set)

    # Then, map read on variants by binary search
    for read in reads:
        read_count += 1
        mapped_variants = []
        for span in read.read_spans:
            left_end = mapped_variants_cache.get(span[0])
            right_end = mapped_variants_cache.get(span[1])
            if left_end is None:
                left_end = binary_search(span[0], variants, 0, variants_len, right_bound = 0)
                mapped_variants_cache[span[0]] = left_end
            if right_end is None:
                right_end = binary_search(span[1], variants, 0, variants_len, right_bound = 1)
                mapped_variants_cache[span[1]] = right_end

            if right_end <= left_end:
                continue
            total_match += right_end - left_end
            taken_read_count += 1
            mapped_variants += variants[left_end:right_end]

        alleles = []
        for var in mapped_variants:
            # get the base pair of given read
            allele = file_utils.get_base_pair_by_var_pos(read, var.end)
            # if the allele=='.', means the loci is deleted/introns or simply with low confidence.
            # we thus consider this var is not on the read.
            if allele == None:
                false_read_matchs += 1
                continue
            times = len(allele)
            allele_read_matchs += 1
            geno = file_utils.get_geno_by_allele(opt, var, allele[0], opt.sep)
            alleles.append(var.unique_id+':'+str(geno)+'*'+str(times))
            variant_allele_map[var.unique_id].add(geno)
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
                edge_raw_read_count = min(int(times_1), int(times_2))
                allele_linkage_map[(allele_1, allele_2)] += 1
                allele_linkage_read_count_map[(allele_1, allele_2)] += edge_raw_read_count
                
                allele_read_count[allele_1] += 1
                allele_read_count[allele_2] += 1
                
                vars_.add(allele_list[i].split(':')[0])
                vars_.add(allele_list[j].split(':')[0])

                if barcode not in edge_barcode_map[(allele_1, allele_2)].keys():
                    edge_barcode_map[(allele_1, allele_2)][barcode] = 0
                edge_barcode_map[(allele_1, allele_2)][barcode] += 1

    print('Received {} reads in total, {} reads are taken, {} reads are omitted. {} matches are found.'\
          .format(read_count, taken_read_count ,read_count - taken_read_count, total_match))
    
    print('There are {} pseudo matches, among which {} are considered matched and {} are false matches, {} variants has at least 1 neighbors.'\
          .format(allele_read_matchs+false_read_matchs, allele_read_matchs, false_read_matchs, len(vars_)))
    
    return allele_linkage_map, edge_barcode_map, len(vars_), allele_linkage_read_count_map, allele_read_count, variant_allele_map
