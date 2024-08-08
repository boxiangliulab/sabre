#!/usr/bin/env python
import tempfile
import gzip
import numpy as np
import subprocess
import os
from rich import print
import collections
from collections import namedtuple
import re

QUAL_THRESHOLD = 0
ALIGNMENT_FILTER = 0

Variant = namedtuple('Variant', ["unique_id", "end", "genotype_string", "is_phased"])
Bamline = namedtuple('Bamline', ['col_pos', 'col_seq', 'col_qual', 'col_cigar', 'alignment_score'])

def create_variant(sep, col_chr, col_pos, col_id, col_ref, col_alt, col_qual, genotype_string, is_phased) -> None:
    unique_id = sep.join([col_chr, col_pos, col_id, col_ref, col_alt])
    return Variant(unique_id=unique_id, end=int(col_pos), genotype_string=genotype_string, is_phased=is_phased)

def get_geno_by_allele(opt, variant, allele, sep):
    if allele == variant.unique_id.split(sep)[-2]:
        return 0 
    elif not opt.non_binary_variant:
        return 1
    elif allele == variant.unique_id.split(sep)[-1]:
       return 1
    else:
       allele = variant.unique_id.split(sep)[-1].split(',').index(allele) + 1 if allele in variant.unique_id.split(sep)[-1] else 2
       return allele

Read = namedtuple('Read', ['umi_barcode', 'read_seqs', 'read_quals', 'read_spans'])

def generate_read_from_bamline(umi_barcode:str, bamline_list:list) -> None:

    global ALIGNMENT_FILTER
    read_seqs = []
    read_quals = []
    read_spans = []
    for bamline in bamline_list:
        if int(bamline.alignment_score) <= ALIGNMENT_FILTER:
            continue
        temp_read_seq, temp_qual_seq, temp_read_span = adjust_read_seq(bamline)
        read_seqs += temp_read_seq
        read_quals += temp_qual_seq
        read_spans += temp_read_span
    return Read(umi_barcode=umi_barcode, read_seqs=tuple(read_seqs), read_quals=tuple(read_quals), read_spans=tuple(read_spans))

def adjust_read_seq(line):
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
    
    return tuple(adjusted_read_seq), tuple(adjusted_qual_seq), tuple(spans)
    

def get_base_pair_by_var_pos(read, var_pos):
    '''
    Given a variant position on genome, return the corresponding base pair on the read.
    '''
    global QUAL_THRESHOLD
    bp_qual_map = collections.defaultdict(list)
    for read, quals, (left, right) in zip(read.read_seqs, read.read_quals, read.read_spans):
        if var_pos <= right and var_pos >= left:
            pos = var_pos - left
            base_pair = read[pos]
            qual = quals[pos]
            bp_qual_map[base_pair].append(qual)

    if len(bp_qual_map) == 1:
        return [list(bp_qual_map.keys())[0] for i in range(len(list(bp_qual_map.values())[0]))] if ord(max(list(bp_qual_map.values()))[0]) - 33 >= QUAL_THRESHOLD else None
    else:
        # pros = functools.reduce(lambda a,b: a*b,list(map(lambda x: 10 ** (-(ord(x)-33)/10),list(bp_qual_map.values())[0])))
        # cons = functools.reduce(lambda a,b: a*b,list(map(lambda x: 10 ** (-(ord(x)-33)/10),list(bp_qual_map.values())[1])))
        # leading_bp = [list(bp_qual_map.keys())[0] for i in range(len(list(bp_qual_map.values())[0]))] if pros < cons else [list(bp_qual_map.keys())[1] for i in range(len(list(bp_qual_map.values())[1]))]
        # print(pros, cons, min(pros, cons)/(pros + cons), str(bp_qual_map))
        pros = sum(list(map(lambda x: -(ord(x)-33)/10, list(bp_qual_map.values())[0])))
        cons = sum(list(map(lambda x: -(ord(x)-33)/10, list(bp_qual_map.values())[1])))
        leading_bp = [list(bp_qual_map.keys())[0] for i in range(len(list(bp_qual_map.values())[0]))] if pros < cons else [list(bp_qual_map.keys())[1] for i in range(len(list(bp_qual_map.values())[1]))]
        if 10 ** -abs(pros-cons) < 0.05:
            return leading_bp
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
    gzip_stream = gzip.open(opt.vcf, 'rt')
    vcf_column_index_map = collections.OrderedDict()
    for line in gzip_stream:
        if line.strip().startswith('#CHR'):
            column_names = line.strip().split('\t')
            vcf_column_index_map = dict(map(lambda x: (x, column_names.index(x)), column_names))
            break
    gzip_stream.close()

    if opt.sample not in vcf_column_index_map:
        raise KeyError("sample {} not in vcf_column_index_map.keys(): {}".format(opt.sample, vcf_column_index_map.keys()))


    # To restrict read-backed phasing to a given chr.
    command_line = ''
    if opt.chr is not None:
        command_line = 'tabix -h ' + opt.vcf + ' ' + opt.chr_vcf + ':'
    else:
        command_line = 'gunzip -c ' + opt.vcf
    
    vcf_out = tempfile.NamedTemporaryFile(delete=False, dir=opt.tmp_dir)
    vcf_out.close()
    processed_vcf_path = vcf_out.name


    # Mean to add black-list to prevent unnecessary calculating.
    # However I don't want to implement this
    if opt.black_list is not None:
        pass
    else:
        command_line += ' | cut -f 1-9,' + str(vcf_column_index_map[opt.sample]+1) + " | grep -v '0|0\|1|1' > " + processed_vcf_path
        err_code = subprocess.check_call("set -euo pipefail && "+command_line, shell=True, executable='/bin/bash')

    if err_code is None:
        raise RuntimeError("Bash returned a err {} when executing {}".format(err_code, command_line))
    
    return processed_vcf_path

def generate_variants(opt, processed_vcf_path):
    '''
    To generate wrapped variants from processed VCF file.
    Notablly, a typical VCF File header looks like:
    #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    SAMPLE_NAME
    '''
    total_record_num = 0
    filtered_record_num = 0

    variants = []
    with open(processed_vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            total_record_num+=1
            columns = line.strip('\n').split('\t')
            col_chr, col_pos, col_id, col_ref, col_alt, col_qual, col_filter, col_info, col_format, col_sample = columns

            # Simple filtering by $FILTER$ == PASS
            if not opt.raw_vcf:
                if 'PASS' not in col_filter:
                    filtered_record_num+=1
                    continue
            else:
                if float(col_qual.strip()) < opt.vcf_qual:
                    filtered_record_num += 1
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
            if opt.neglect_hla and opt.chr=='chr6':
                if int(col_pos) >= 29602228 and int(col_pos) <= 33410226:
                    filtered_record_num += 1
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
                variants.append(create_variant(opt.sep, col_chr, col_pos, col_id, col_ref, col_alt, col_qual, genotype_string, is_phased))
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

    bam_files = []

    output_sam_paths = []
    if opt.bam_list is not None:
        with open(opt.bam_list) as f:
            bam_files = list(map(lambda x: x.strip().split(',') ,f.readlines()))
    else:
        bam_files.append([opt.bam, '1'])

    for bam_file, name in bam_files:
        bam_out = tempfile.NamedTemporaryFile(delete=False, dir=opt.tmp_dir)
        bam_out.close()
        output_sam_paths.append([bam_out.name, name])

        # We don't need headers...for now.
        command = "samtools view {} '{}' -F 0x400 -@ 6 -q {} -L {} > {}".format(bam_file, opt.chr, opt.mapq_threshold, bed_file, bam_out.name)
        err_code = subprocess.check_call("set -euo pipefail && "+command, shell=True, executable='/bin/bash')

    os.remove(bed_file)
    return output_sam_paths

def generate_reads(opt, output_sam_paths):
    '''
    Generate reads from filtered SAM file from give BAM file.
    Notablly, each line in the BAM file is wrapped as a Bamline.
    Bamlines with the same umi-barcode consists one Read.
    '''
    umibarcode_line_map = collections.defaultdict(list)
    alignment_scores = []
    bamline_cnt = 0
    global QUAL_THRESHOLD, ALIGNMENT_FILTER
    QUAL_THRESHOLD = 10

    if opt.input_type == 're':
        bc_re = re.compile(opt.bc_re)
        umi_re = re.compile(opt.umi_re)

    for output_sam_path, name in output_sam_paths:
        with open(output_sam_path) as sam_file:
            for line in sam_file:
                columns = line.strip('\n').split('\t')
                col_qname, col_flag, col_rname, col_pos, col_mapq, col_cigar, col_rnext, col_pnext, col_tlen, col_seq, col_qual = columns[:11]
                other_columns = columns[11:]
                alignment_score = -100
                col_umi, col_barcode = None, None
                for col in other_columns:
                    if col.startswith('AS'):
                        alignment_score = int(col.split(':')[-1])
                    # umi
                    elif col.startswith('UB'):
                        col_umi = col.split(':')[-1]
                    # barcode
                    elif col.startswith('CB'):
                        col_barcode = col.split(':')[-1]
                bamline = Bamline(col_pos, col_seq, col_qual, col_cigar, alignment_score)
                if opt.input_type == 'cellranger':
                    if col_umi is None or col_barcode is None:
                        continue
                    else:
                        barcode, umi = col_barcode, col_umi
                elif opt.input_type == 'umitools':
                    # SRR8551677.318476227_CTAATGGAGACTAAGT_TCCAGACCGG
                    _, barcode, umi = col_qname.split('_')
                elif opt.input_type == 'star':
                    # AGATGTACTATCAGCAACATTGGC_GTGAGGACTT_AAAAAEEEEE_SRR6750053.36514992
                    barcode, umi = col_qname.split('_')[:2]
                elif opt.input_type == 're':
                    barcode, umi = bc_re.findall(line)[0], umi_re.findall(line)[0]
                    
                # barcode, umi = str(bamline_cnt), str(bamline_cnt)
                umi_barcode = '.'.join([umi, barcode]) + '_{}'.format(name)
                umibarcode_line_map[umi_barcode].append(bamline)
                bamline_cnt += 1
                alignment_scores.append(int(alignment_score))
            os.remove(output_sam_path)

        # Filter these reads by alignment score
        ALIGNMENT_FILTER = np.percentile(alignment_scores, opt.as_quality*100)
        for umi_barcode, lines in umibarcode_line_map.items():
            yield generate_read_from_bamline(umi_barcode=umi_barcode, bamline_list=lines)

def generate_bed_file(opt, variants:list[Variant]):
    '''
    Generate BED file, which contains var, var_start, var_end
    '''
    bed_file = tempfile.NamedTemporaryFile(delete=False, mode='wt', dir=opt.tmp_dir)
    sep = opt.sep
    for var in variants:
        bed_file.write('{}\t{}\t{}\n'.format(sep.join(var.unique_id.split(sep)[:-4]), var.end-1, var.end))
    bed_file.close()
    return bed_file.name


