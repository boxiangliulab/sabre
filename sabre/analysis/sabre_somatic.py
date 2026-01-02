import networkx as nx
import os
import pickle

import pysam
import argparse

import subprocess
import bisect
import re
import functools

import collections


import time

chromosome_map = {}
gene_dict = {}

def traverse_graph(G, start_node):
    visited = set()  # 避免重复访问
    result = set()  # 存储 germline 节点
    stack = [start_node]  # 使用栈进行深度优先搜索 (DFS)
    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        for neighbor in G.neighbors(node):
            if 'somatic' in neighbor and ':1' in neighbor:
                return []  # 规则 2：如果遇到 somatic:1，则失败
            elif 'somatic' in neighbor and ':0' in neighbor:
                stack.append(neighbor)  # 规则 3：继续遍历 somatic:0
            else:
                result.add(neighbor)  # 规则 3：如果是 germline，加入结果列表
    return list(result)

def find_somatic_variants(opt):
    graph_base_dir = '{}/{}/conflict_graphs/'.format(opt.output_dir, opt.id)
    variant_infos = []

    for graph_path in os.listdir(graph_base_dir):
        if 'somatic' not in graph_path: continue
        if graph_path.split('.')[0] != opt.chr: continue
        try:
            G = pickle.load(open(os.path.join(graph_base_dir, graph_path), 'rb'))
            for node in G.nodes:
                if 'somatic' in node:
                    chr_, pos, _, ref_g, alt_g = node.split(':')[0].split('_')
                    pos = int(pos)
                    result = find_gene(chr_, pos)
                    if result is None: continue
                    gene, (chr_, start, end) = result
                    variant_infos.append(((chr_, pos, ref_g, alt_g), gene, (start, end), os.path.join(graph_base_dir, graph_path)))
        except Exception as e:
            print(e)
    
    return variant_infos

def find_gene(chromosome, pos):
    global chromosome_map, gene_dict

    if chromosome not in chromosome_map:
        return None  # 该染色体没有基因
    
    gene_list = chromosome_map[chromosome]  # 获取该染色体的基因区间
    starts = [start for start, _, _ in gene_list]  # 仅取 start 进行二分查找
    
    # 二分查找第一个大于pos的start索引
    idx = bisect.bisect_right(starts, pos) - 1
    
    if idx >= 0:
        start, end, gene = gene_list[idx]
        if start <= pos <= end:
            return gene, gene_dict[gene]  # pos 落在基因区域
    
    return None

def analysis_somatic_variant(opt, variant_infos):
    in_phase_germline_missense_pairs = set()
    out_phase_germline_missense_pairs = set()
    chromosome = opt.chr
    vcf = pysam.VariantFile(opt.vcf)
    for (chr_, pos, ref_g, alt_g), gene, (start, end), graph_path in variant_infos:

        sG = pickle.load(open(graph_path, 'rb'))

        somatic_variant = f'{chr_}_{pos}_somatic_{ref_g}_{alt_g}'
        if somatic_variant+':1' not in sG.nodes or somatic_variant+':0' not in sG.nodes:
            continue

        lines = list(vcf.fetch(chromosome, start, end))

        upper_half = []
        lower_half = []

        germline_variants = []

        for line in lines:
            if not line: continue
            if line.samples[0]['GT'] in [(0, 0), (1, 1)]:
                continue
            germline_variants.append('_'.join([chr_, str(line.pos), '.', line.ref, line.alts[0]])+':1')
            left, right = line.samples[0]['GT']
            upper_half.append('_'.join([chr_,str(line.pos),'.',line.ref, line.alts[0]]) + f':{left}')
            lower_half.append('_'.join([chr_,str(line.pos),'.',line.ref, line.alts[0]]) + f':{right}')

        upper_half_set = set(upper_half)
        lower_half_set = set(lower_half)

        if len(upper_half_set) <= 1: continue
        
        upper_half_dict = [dict() for i in upper_half]
        lower_half_dict = [dict() for i in lower_half]

        if len(germline_variants) == 0: continue

        for missense_variant in germline_variants:
            for index, upper_var in enumerate(upper_half):
                if missense_variant in upper_half_set:
                    upper_half_dict[index][missense_variant.split(':')[0]] = True
                    lower_half_dict[index][missense_variant.split(':')[0]] = False
                else:
                    upper_half_dict[index][missense_variant.split(':')[0]] = False
                    lower_half_dict[index][missense_variant.split(':')[0]] = True


        G = nx.Graph()

        for node, data in zip(upper_half, upper_half_dict):
            G.add_node(node)
            G.nodes[node]['missense'] = data

        for node, data in zip(lower_half, lower_half_dict):
            G.add_node(node)
            G.nodes[node]['missense'] = data

        for i in range(len(upper_half)-1):
            # for j in range(1, len(lower_half)):
            G.add_edge(upper_half[i], upper_half[i+1])
            G.add_edge(lower_half[i], lower_half[i+1])

        somatic_node = somatic_variant + ':0'
        for neighbor in sG.neighbors(somatic_node):
            G.add_edge(somatic_node, neighbor)
            G[somatic_node][neighbor]['barcodes'] = sG[somatic_node][neighbor]['barcodes']
            G[somatic_node][neighbor]['raw_read_count'] = sG[somatic_node][neighbor]['raw_read_count']
        
        somatic_node = somatic_variant + ':1'
        for neighbor in sG.neighbors(somatic_node):
            G.add_edge(somatic_node, neighbor)
            G[somatic_node][neighbor]['barcodes'] = sG[somatic_node][neighbor]['barcodes']
            G[somatic_node][neighbor]['raw_read_count'] = sG[somatic_node][neighbor]['raw_read_count']

        in_phase_normal_cells = set()
        in_phase_mutate_cells = set()

        out_phase_normal_cells = set()
        out_phase_mutate_cells = set()

        for missense_variant in germline_variants:
            missense_variant = missense_variant.split(':')[0]
            somatic_alt_connect_to_missense = []
            for neighbor in traverse_graph(G, somatic_variant+':1'):
                if 'somatic' in neighbor: continue
                if 'missense' not in G.nodes[neighbor]:
                    somatic_alt_connect_to_missense = []
                    break
                somatic_alt_connect_to_missense.append(G.nodes[neighbor]['missense'][missense_variant])

            if (any(somatic_alt_connect_to_missense) and not all(somatic_alt_connect_to_missense)) or len(somatic_alt_connect_to_missense) == 0: 
                continue

            support_cells = []
            support_weights = 0
            for neighbor in nx.neighbors(G, somatic_variant + ':1'):
                barcode = G[somatic_variant+':1'][neighbor]['barcodes']
                support_cells+=list(barcode.keys())
                support_weights+=G[somatic_variant+':1'][neighbor]['raw_read_count']
            
            if len(support_cells) < opt.cells or support_weights < opt.reads: continue

            if all(somatic_alt_connect_to_missense):
                # in phase
                for node in traverse_graph(G, somatic_variant+':0'):
                    if 'somatic' in node: continue
                    if 'missense' not in G.nodes[node]: continue
                    if G.nodes[node]['missense'][missense_variant]:
                        in_phase_normal_cells |= set(G[somatic_variant+':0'][node]['barcodes'].keys())
                for node in traverse_graph(G, somatic_variant+':1'):
                    if 'somatic' in node: continue
                    if 'missense' not in G.nodes[node]: continue
                    if G.nodes[node]['missense'][missense_variant]:
                        in_phase_mutate_cells |= set(G[somatic_variant+':1'][node]['barcodes'].keys())
            in_phase_germline_missense_pairs.add((opt.id, missense_variant, somatic_variant, gene, ','.join(support_cells), support_weights))

            if not any(somatic_alt_connect_to_missense):
                # out of phase
                for node in traverse_graph(G, somatic_variant+':0'):
                    if 'somatic' in node: continue
                    if 'missense' not in G.nodes[node]: continue
                    if not G.nodes[node]['missense'][missense_variant]:
                        out_phase_normal_cells |= set(G[somatic_variant+':0'][node]['barcodes'].keys())
                for node in traverse_graph(G, somatic_variant+':1'):
                    if 'somatic' in node: continue
                    if 'missense' not in G.nodes[node]: continue
                    if not G.nodes[node]['missense'][missense_variant]:
                        out_phase_mutate_cells |= set(G[somatic_variant+':1'][node]['barcodes'].keys())
            out_phase_germline_missense_pairs.add((opt.id, missense_variant, somatic_variant, gene, ','.join(support_cells), support_weights))

    with open(f'{opt.output_dir}/{opt.id}/{opt.chr}.in.phase.details.tsv', 'w') as f:
        f.write('sample\tgermline\tsomatic\tgene\tcells\tsupports\n')
        for (sample, missense_variant, somatic_variant, gene, cells, supports) in in_phase_germline_missense_pairs:
            f.write(f'{sample}\t{missense_variant}\t{somatic_variant}\t{gene}\t{cells}\t{supports}\n')

    with open(f'{opt.output_dir}/{opt.id}/{opt.chr}.out.of.phase.details.tsv', 'w') as f:
        f.write('sample\tgermline\tsomatic\tgene\tcells\tsupports\n')
        for (sample, missense_variant, somatic_variant, gene, cells, supports) in out_phase_germline_missense_pairs:
            f.write(f'{sample}\t{missense_variant}\t{somatic_variant}\t{gene}\t{cells}\t{supports}\n')

def run(opt):

    global chromosome_map, gene_dict

    # Check Arguments
    if not os.path.exists(opt.vcf):
        raise ValueError(f"VCF file not found: {opt.vcf}. Please check input argument --vcf")

    if not os.path.exists(os.path.join(opt.output_dir, opt.id)):
        raise ValueError("File not found: {}. Please check input argument --id and --output_dir".format(os.path.join(opt.output_dir, opt.id)))
    
    if not opt.gtf.endswith('.gtf') and not opt.gtf.endswith('.gtf.gz'):
        raise ValueError("Invalid argument for --gtf. Please give a valid .gtf file or .gtf.gz file.")
    
    if not os.path.exists(opt.gtf):
        raise ValueError("File not found: {}. Please check input argument --gtf".format(opt.gtf))

    # load chromosome_map, gene_dict
    if not os.path.exists(f'{opt.gtf}.gene.dict'):
        init(opt)
    

    gene_dict = pickle.load(open(f'{opt.gtf}.gene.dict', 'rb'))
    chromosome_map = pickle.load(open(f'{opt.gtf}.chr.dict', 'rb'))

    variant_infos = find_somatic_variants(opt)
    analysis_somatic_variant(opt, variant_infos)


def dfs_from_node(G, start_node, visited=None, path=None):
    if visited is None:
        visited = set()
    if path is None:
        path = []

    visited.add(start_node.split(':')[0])
    path = path + [start_node]
    neighbors = [n for n in G.neighbors(start_node) if n.split(':')[0] not in visited]
    if not neighbors:
        return [path]
    all_paths = []
    for neighbor in neighbors:
        sub_paths = dfs_from_node(G, neighbor, visited.copy(), path)
        all_paths.extend(sub_paths)
    return all_paths

def test_graph(G: nx.Graph):
    somatic_set = []
    for node in G.nodes:
        if 'somatic' in node: 
            somatic_set.append(node)
    putative_somatic_variant_edge_count_map = collections.defaultdict(lambda: [0, 0])
    for node in somatic_set:
        variant, geno = node.split(':')
        putative_somatic_variant_edge_count_map[variant][int(geno)] = len([n for n in G.neighbors(node)])

    if len(set(map(lambda x: x.split(':')[0], somatic_set))) > 1 or len(somatic_set) != 2:
        return set(), set(map(lambda x: x.split(':')[0],(filter(lambda x: 'somatic' in x, G.nodes)))), putative_somatic_variant_edge_count_map
    
    variant_color_map = collections.defaultdict(set)
    first_visit_paths = dfs_from_node(G, somatic_set[0])
    blue_pills = list(set(functools.reduce(lambda x, y: x+y, first_visit_paths)))
    for bp in blue_pills:
        variant, geno = bp.split(':')
        variant_color_map[variant].add(f'blue:{geno}')

    second_visit_paths = dfs_from_node(G, somatic_set[1])
    red_pills = list(set(functools.reduce(lambda x, y: x+y, second_visit_paths)))
    for bp in red_pills:
        variant, geno = bp.split(':')
        variant_color_map[variant].add(f'red:{geno}')

    for variant, pill_box in variant_color_map.items():
        count = len(pill_box)
        if count == 4:
            return set(), set(map(lambda x: x.split(':')[0],(filter(lambda x: 'somatic' in x, G.nodes)))), putative_somatic_variant_edge_count_map

    return set(map(lambda x: x.split(':')[0],(filter(lambda x: 'somatic' in x, G.nodes)))), set(), putative_somatic_variant_edge_count_map


def filter_(opt):
    graph_base_dir = '{}/{}/conflict_graphs/'.format(opt.output_dir, opt.id)
    output_result_list = []
    for graph_path in os.listdir(graph_base_dir):
        if 'somatic' not in graph_path: continue

        try:
            G = pickle.load(open(os.path.join(graph_base_dir, graph_path), 'rb'))
            sm_set, se_set, edge_count_map = test_graph(G)
            for sm in sm_set:
                output_result_list.append(f'{opt.id}\t{sm}\t{max(edge_count_map[sm])}\t{min(edge_count_map[sm])}\t{1}\n')

            for sm in se_set:
                output_result_list.append(f'{opt.id}\t{sm}\t{max(edge_count_map[sm])}\t{min(edge_count_map[sm])}\t{0}\n')
        except Exception as e:
            print(e)
    
    with open(f'{opt.output_dir}/{opt.id}/refined.somatic.mutations.tsv', 'w') as f:
        f.write('ID\tvariant\Ref\tAlt\tResult\n')
        for line in output_result_list:
            f.write(line)

def init(opt):
    gene_name_re = re.compile('gene_name "(.*?)";')
    gene_dict = {}
    with open(opt.gtf) as f:
        for line in f:
            if 'HAVANA\tgene' not in line:
                continue
            gene_name = gene_name_re.search(line).group(1)
            chr_, _, _, start, end = line.split('\t')[:5]
            gene_dict[gene_name] = (chr_, int(start), int(end))

    chromosome_map = {}
    for gene, (chr_, start, end) in gene_dict.items():
        if chr_ not in chromosome_map:
            chromosome_map[chr_] = []
        chromosome_map[chr_].append((start, end, gene))

    for chr_ in chromosome_map:
        chromosome_map[chr_].sort()
    
    pickle.dump(gene_dict, open(f'{opt.gtf}.gene.dict', 'wb'))
    pickle.dump(chromosome_map, open(f'{opt.gtf}.chr.dict', 'wb'))


def main():
    parser = argparse.ArgumentParser(prog='sabre-somatic')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # 子命令 init
    parser_init = subparsers.add_parser('init', help='Initialize sabre-somatic and generate gtf index file.')
    parser_init.add_argument('--gtf', help='Input GTF file.', required=True)
    parser_init.set_defaults(func=init)

    # 子命令 init
    parser_filter = subparsers.add_parser('filter', help='Filter the somatic mutation callset with phasing information.')
    parser_filter.add_argument("--id", help="Input ID", required=True)
    parser_filter.add_argument("--output_dir", help="Path to output directory, should be the same as sabre --output_dir. Default ./output", required=False, default='./output')
    parser_filter.set_defaults(func=filter_)

    # 子命令 run
    parser_run = subparsers.add_parser('run', help='Run sabre-somatic analysis')
    parser_run.add_argument("--id", help="Input ID", required=True)
    parser_run.add_argument("--vcf", help="The germline VCF file of `id`", required=True)
    parser_run.add_argument("--chr", help="Indicate the chromosome, on which the sabre-somatic will perform one-two hit analysis")
    parser_run.add_argument("--output_dir", help="Path to output directory, should be the same as sabre --output_dir. Default ./output", required=False, default='./output')
    parser_run.add_argument("--gtf", help="Input GTF file. We recommend you to input a .gtf file rather than a .gtf.gz file, because sabre-somatic will have to depress the .gtf.gz file everytime you input a compressed gtf file.", required=True)
    parser_run.add_argument("--cells", help="Threshold on number of supporting cells.", default=1, type=int)
    parser_run.add_argument("--reads", help="Threshold on number of supporting reads", default=5, type=int)
    parser_run.set_defaults(func=run)

    # 解析并调用对应的函数
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
    
