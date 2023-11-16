#!/usr/bin/env python
import tempfile
import gzip
import numpy as np
import subprocess
import os
from rich import print
import collections
import itertools

class Variant:
    def __init__(self, col_chr, col_pos, col_id, col_ref, col_alt, col_qual, genotype_string, is_phased) -> None:
        
        self.col_chr = col_chr
        self.col_pos = col_pos
        self.col_id = col_id
        self.col_ref = col_ref
        self.col_alt = col_alt
        self.col_qual = col_qual
        self.genotype_string = genotype_string
        self.is_phased = is_phased

        self.start = int(self.col_pos)-1
        self.end = int(self.col_pos)

        self.unique_id = self.create_unique_id()

    def create_unique_id(self):
        return '_'.join([self.col_chr, self.col_pos, self.col_id, self.col_ref, self.col_alt])
    
    def create_node(self):
        return self.unique_id+':0', self.unique_id+':1'
    
    def get_geno_by_allele(self, allele):
        if allele == self.col_ref:
            return 0 
        return 1
    
    def get_allele_by_geno(self, geno):
        if int(geno) == 1:
            return self.col_alt
        return self.col_ref

class Bamline:
    def __init__(self,col_qname, col_flag, col_rname, col_pos, col_mapq, col_cigar,\
                  col_rnext, col_pnext, col_tlen, col_seq, col_qual, alignment_score,\
                      col_rx, col_qx, col_bx, col_bc, col_qt) -> None:
        self.col_qname = col_qname
        self.col_flag = col_flag
        self.col_rname = col_rname
        self.col_pos = col_pos
        self.col_mapq = col_mapq
        self.col_cigar = col_cigar
        self.col_rnext = col_rnext
        self.col_pnext = col_pnext
        self.col_tlen = col_tlen
        self.col_seq = col_seq
        self.col_qual = col_qual
        self.alignment_score = alignment_score
        self.col_rx = col_rx
        self.col_qx = col_qx
        self.col_bx = col_bx
        self.col_bc = col_bc
        self.col_qt = col_qt


class Read:

    lengths = []
    def __init__(self, umi_barcode:str, bamline_list:list[Bamline], alignment_filter, qual_threshold=20, neglect_overlap=True) -> None:
        self.umi_barcode = umi_barcode
        self.umi, self.barcode = self.umi_barcode.split('.')
        self.bamline_list = bamline_list
        self.qual_threshold = qual_threshold
        self.alignment_filter = alignment_filter
        self.neglect_overlap = True
        Read.lengths.append(len(self.bamline_list))

        self.extract_read_seq()
        self.filter_seq_by_quality()
        if len(self.read_seqs) > 1:
            self.merge_seqs()

    def extract_read_seq(self):
        '''
        Extract seq and span from each bamline by adjust_read_seq, extracted read_seqs, qual_seqs and covering_spans are stored.
        '''
        self.read_seqs = []
        self.read_quals = []
        self.read_spans = []
        for bamline in self.bamline_list:
            if int(bamline.alignment_score) <= self.alignment_filter:
                continue
            temp_read_seq, temp_qual_seq, temp_read_span = Read.adjust_read_seq(bamline)
            self.read_seqs += temp_read_seq
            self.read_quals += temp_qual_seq
            self.read_spans += temp_read_span

    def filter_seq_by_quality(self):
        '''
        We first filter the read by quality
        Any read base pair with QUAL lower than a given quality threshold will be replaced by '.'
        '''
        filtered_seqs = []

        for seq, qual in zip(self.read_seqs, self.read_quals):
            temp_seq = ''
            for i, q in enumerate(qual):
                if ord(q) - 33 < self.qual_threshold:
                    temp_seq += '.'
                else:
                    temp_seq += seq[i]
            filtered_seqs.append(temp_seq)

        self.read_seqs = filtered_seqs

    def adjust_read_seq(line:Bamline):
        '''
        To align variants onto reads, we need to get the real coverage of the given read to the genome, by analyzing the CIGAR string.
        For every bamline, or bamlines splited by $N$, we extract each corresponding read_seq, qual_seq and covering_span from it.
        Thus, a single read may contain multiple read_seqs, qual_seqs and covering_spans, determined by the number bamlines sharing the same umi_barcode,
        and the number of $N$s in CIGAR.
        '''
        adjusted_read_seq = []
        adjusted_qual_seq = []
        spans = []
        current_cigar_number = 0
        genome_pointer = int(line.col_pos)
        seq_pointer = 0
        for char in line.col_cigar:
            if char in '0123456789':
                current_cigar_number = current_cigar_number*10 + int(char)
            else:
                if char == 'S':
                    # soft clipping
                    # consumes query
                    # just ignore the corresponding seq by incrementing the seq pointer.
                    seq_pointer += current_cigar_number
                elif char == 'M' or char == '=' or char == 'X':
                    # match
                    # consumes query and reference
                    # put the corresponding seq to the adjusted_read_seq, and adjust read_pointer
                    adjusted_read_seq.append(line.col_seq[seq_pointer:seq_pointer+current_cigar_number])
                    adjusted_qual_seq.append(line.col_qual[seq_pointer:seq_pointer+current_cigar_number])
                    # [a, b] closed intervals
                    spans.append((genome_pointer, genome_pointer+current_cigar_number-1))
                    seq_pointer += current_cigar_number
                    genome_pointer += current_cigar_number
                elif char == 'I':
                    # insertion
                    # consumes query
                    # ignore the insertion by adding seq_pointer, the adjusted_read_seq will not be affected.
                    seq_pointer += current_cigar_number
                elif char == 'D' or char == 'N':
                    # deletion or skipped region
                    # consumes reference
                    # ignore it and move the pointer of genome right.
                    genome_pointer += current_cigar_number

                # No operations on "H" nor "P", cuz neither "H" nor "P" consumes neither
                current_cigar_number = 0
        
        return adjusted_read_seq, adjusted_qual_seq, spans
    
    def merge_read_with_qual(read_pair, span_pair, qual_pair, neglect_overlap):
        '''
        merge a pair of reads, the final base pair is determined by corresponding quality

        span_pair: [a, b], [c, d]
        where a <= c <= b
        '''
        read_1, read_2 = read_pair
        (a, b), (c, d) = span_pair
        qual_1, qual_2 = qual_pair

        conflict_reads = ''
        conflict_quals = ''
        # overlapping area
        for i in range(c, min(b, d) + 1):
            base_1, q_1 = read_1[i-a], qual_1[i-a]
            base_2, q_2 = read_2[i-c], qual_2[i-c]
            # conflict_reads += base_1 if q_1 > q_2 else base_2
            conflict_reads += '.'
            conflict_quals += max(q_1, q_2)

        return read_1[:c-a] + conflict_reads+ (read_2[b-c+1:d-c+1] if d>b else read_1[d-a+1:b-a+1]), \
        qual_1[:c-a] + conflict_quals+ (qual_2[b-c+1:d-c+1] if d>b else qual_1[d-a+1:b-a+1]), (a, max(b, d))

    def merge_seqs(self):
        '''
        merge reads with overlapping spans.
        '''
        merge_result = []
        mixed_input = list(zip(self.read_seqs, self.read_quals, self.read_spans))
        mixed_input = sorted(mixed_input, key=lambda x: x[2][0], reverse=True)
        while len(mixed_input) > 1:
            (read_1, qual_1, span_1) = mixed_input.pop()
            (read_2, qual_2, span_2) = mixed_input.pop()
            if span_2[0] > span_1[1]:
                # no overlap
                merge_result.append((read_1, qual_1, span_1))
                mixed_input.append((read_2, qual_2, span_2))
            else:
                # share overlap
                mixed_input.append(Read.merge_read_with_qual((read_1, read_2), (span_1, span_2), (qual_1, qual_2), self.neglect_overlap))
        
        merge_result += mixed_input
        self.read_seqs, self.read_quals, self.read_spans = list(zip(*merge_result))


    def get_base_pair_by_var_pos(self, var_pos):
        '''
        Given a variant position on genome, return the corresponding base pair on the read.
        '''

        for read, (left, right) in zip(self.read_seqs, self.read_spans):
            if var_pos <= right and var_pos >= left:
                pos = var_pos - left
                base_pair = read[pos]

                if base_pair == '.':
                    # Means the variant actually do not appear in this read.
                    # It is covered by a deletion
                    return None
                return base_pair
        return None

def load_barcode_list(opt):
    '''
    Load whitelisted barcodes from given barcodes.csv
    '''
    barcode_file = opt.barcode_path
    whitelisted_barcodes = set()
    with open(barcode_file, 'r') as f:
        for bc in f.readlines():
            whitelisted_barcodes.add(bc.strip())
    return whitelisted_barcodes
    
def load_vcf(opt):
    '''
    To load and do simple filtering on vcf file.
    '''
    # First try to figure out which column the sample corresponds to in vcf file.
    # By literally get all the column names.
    gzip_stream = gzip.open(opt.vcf_path, 'rt')
    vcf_column_index_map = collections.OrderedDict()
    for line in gzip_stream:
        if line.strip().startswith('#CHR'):
            column_names = line.strip().split('\t')
            vcf_column_index_map = dict(map(lambda x: (x, column_names.index(x)), column_names))
            break
    gzip_stream.close()

    if opt.sample_name not in vcf_column_index_map:
        raise KeyError("sample_name {} not in vcf_column_index_map.keys(): {}".format(opt.sample_name, vcf_column_index_map.keys()))


    # To restrict read-backed phasing to a given chr.
    command_line = ''
    if opt.restrict_chr is not None:
        command_line = 'tabix -h ' + opt.vcf_path + ' ' + opt.restrict_chr_vcf + ':'
    else:
        command_line = 'gunzip -c ' + opt.vcf_path
    
    vcf_out = tempfile.NamedTemporaryFile(delete=False)
    vcf_out.close()
    processed_vcf_path = vcf_out.name


    # Mean to add black-list to prevent unnecessary calculating.
    # However I don't want to implement this
    if opt.black_list is not None:
        pass
    else:
        command_line += ' | cut -f 1-9,' + str(vcf_column_index_map[opt.sample_name]+1) + " | grep -v '0|0\|1|1' > " + processed_vcf_path
        err_code = subprocess.check_call("set -euo pipefail && "+command_line, shell=True, executable='/bin/bash')

    if err_code is None:
        raise RuntimeError("Bash returned a err {} when executing {}".format(err_code, command_line))
    
    return processed_vcf_path


def generate_variants(processed_vcf_path):
    '''
    To generate wrapped variants from processed VCF file.
    Notablly, a typical VCF File header looks like:
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_NAME
    '''
    total_record_num = 0
    filtered_record_num = 0

    variants = []
    with open(processed_vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.strip('\n').startswith('#'):
                continue

            total_record_num+=1
            columns = line.strip('\n').split('\t')
            col_chr, col_pos, col_id, col_ref, col_alt, col_qual, col_filter, col_info, col_format, col_sample = columns

            # Simple filtering by $FILTER$ == PASS
            if 'PASS' not in col_filter:
                filtered_record_num+=1
                continue

            # Try to get annotated genotype in VCF
            is_phased = False
            format = col_format.split(':')
            GT_index = None
            if 'GT' in format:
                GT_index = format.index('GT')
            else:
                filtered_record_num+=1
                continue
            
            # Try to get Genotype string in VCF, w.r.t. 1/0 or 0/1 or 1|0 or 0|1
            genotype_string = col_sample.split(':')[GT_index]
            if '.' in genotype_string:
                filtered_record_num += 1
                continue
            elif '|' in genotype_string:
                is_phased = True
            elif '/' in genotype_string:
                if genotype_string.split('/')[0] == genotype_string.split('/')[1]:
                    filtered_record_num += 1
                    continue
            
            # Restrict the length of alternative base pair into 1.
            all_alleles = col_alt.split(',') + [col_ref]
            if max([len(x) for x in all_alleles]) == 1:
                variants.append(Variant(col_chr, col_pos, col_id, col_ref, col_alt, col_qual, genotype_string, is_phased))
            else:
                filtered_record_num += 1
                continue
    print('Received {} variants in total, {} variants taken, {} variants omitted.'.\
          format(total_record_num, total_record_num-filtered_record_num, filtered_record_num))
    
    os.remove(processed_vcf_path)
    
    vid_var_map = collections.OrderedDict()
    for var in variants:
        vid_var_map[var.unique_id] = var

    return variants, vid_var_map
           

def load_bam(opt, bed_file):
    '''
    To load and do simple filtering on BAM files.
    '''
    bam_out = tempfile.NamedTemporaryFile(delete=False)
    bam_out.close()
    output_sam_path = bam_out.name

    # We don't need headers...for now.
    command = "samtools view {} '{}' -F 0x400 -q {} -L {} > {}".format(opt.bam_path, opt.restrict_chr, opt.mapq_threshold,\
                                                                             bed_file, output_sam_path)
    err_code = subprocess.check_call("set -euo pipefail && "+command, shell=True, executable='/bin/bash')

    os.remove(bed_file)

    return output_sam_path


def generate_reads(opt, output_sam_path):
    '''
    Generate reads from filtered SAM file from give BAM file.
    Notablly, each line in the BAM file is wrapped as a Bamline.
    Bamlines with the same umi-barcode consists one Read.
    '''
    umibarcode_line_map = collections.defaultdict(list)
    sam_file = open(output_sam_path, mode='r')
    alignment_scores = []
    bamline_cnt = 0
    barcode_umi_cnt = collections.defaultdict(int)
    for line in sam_file.readlines():
        columns = line.strip('\n').split('\t')
        col_qname, col_flag, col_rname, col_pos, col_mapq, col_cigar, col_rnext, col_pnext, col_tlen, col_seq, col_qual = \
            columns[:11]
        other_columns = columns[11:]
        alignment_score = -100
        col_rx, col_qx, col_umi, col_barcode, col_qt = None, None, None, None, None
        for col in other_columns:
            if col.startswith('AS'):
                alignment_score = int(col.split(':')[-1])
            elif col.startswith('RX'):
                col_rx = col.split(':')[-1]
            elif col.startswith('QX'):
                col_qx = col.split(':')[-1]
            # umi
            elif col.startswith('UB'):
                col_umi = col.split(':')[-1]
            # barcode
            elif col.startswith('CB'):
                col_barcode = col.split(':')[-1]
            elif col.startswith('QT'):
                col_qt = col.split(':')[-1]
        line = Bamline(col_qname, col_flag, col_rname, col_pos, col_mapq, col_cigar, col_rnext, col_pnext, col_tlen, col_seq,\
                     col_qual, alignment_score, col_rx, col_qx, col_umi, col_barcode, col_qt)
        if opt.input_type == 'cellranger':
            if col_umi is None or col_barcode is None :
                continue
            else:
                barcode, umi = col_barcode, col_umi
        elif opt.input_type == 'umitools':
            # SRR8551677.318476227_CTAATGGAGACTAAGT_TCCAGACCGG
            _, barcode, umi = col_qname.split('_')
        umibarcode_line_map['.'.join([umi, barcode])].append(line)
        # umibarcode_line_map['.'.join([str(bamline_cnt),str(bamline_cnt)])].append(line)
        bamline_cnt += 1
        barcode_umi_cnt[barcode] += 1
        alignment_scores.append(int(alignment_score))
    os.remove(output_sam_path)

    umi_cnt_threshold = 0
    filtered_barcode = set(filter(lambda x: barcode_umi_cnt[x]>umi_cnt_threshold, barcode_umi_cnt.keys()))

    print('Received {} barcodes in total, after filtered #UMI less than {}, {} barcodes taken, {} barcodes omitted.'.\
          format(len(barcode_umi_cnt.keys()), umi_cnt_threshold, len(filtered_barcode), len(barcode_umi_cnt.keys())-len(filtered_barcode)))

    # Filter these reads by alignment score
    alignment_score_filter = np.percentile(alignment_scores, opt.as_quality*100)
    # alignment_score_filter = -99
    reads = []
    for umi_barcode, lines in umibarcode_line_map.items():
        if umi_barcode.split('.')[-1] not in filtered_barcode:
            continue
        read = Read(umi_barcode, lines, alignment_score_filter)
        reads.append(read)
    return reads

def generate_bed_file(opt, variants:list[Variant]):
    '''
    Generate BED file, which contains var, var_start, var_end
    '''
    bed_file = tempfile.NamedTemporaryFile(delete=False, mode='wt')
    for var in variants:
        bed_file.write('{}\t{}\t{}\n'.format(opt.restrict_chr, var.start, var.end))
    bed_file.close()
    return bed_file.name