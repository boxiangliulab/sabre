import networkx as nx
import os
import pickle

in_phase_output = open('./in.phase.hits.txt', 'w')
out_of_phase_output = open('./out.of.phase.hits.txt', 'w')

def examine_in_phase(lbd_list):
    for lbd in lbd_list:
        print('Processing {}...'.format(lbd))
        graph_path = './output/{}/conflict_graphs/'.format(lbd)
        if not os.path.exists(graph_path): return
        graphs = os.listdir(graph_path)
        for graph in graphs:
            try:
                G = pickle.load(open(graph_path+graph, 'rb'))
            except:
                continue
            variants = set()
            list(map(lambda x: variants.add(x.split(':')[0]), list(G.nodes)))
            variants = list(variants)
            for i in range(len(variants)-1):
                for j in range(i, len(variants)):
                    var1, var2 = variants[i], variants[j]
                    if (var1+':0', var2+':0') in G.edges and (var1+':1', var2+':0') in G.edges and (var1+':1', var2+':1') in G.edges and (var1+':0', var2+':1') not in G.edges:
                        _00_cell_count = len(G.edges[var1+':0', var2+':0']['barcodes'].keys())
                        _01_cell_count = len(G.edges[var1+':1', var2+':0']['barcodes'].keys())
                        _11_cell_count = len(G.edges[var1+':1', var2+':1']['barcodes'].keys())
                        try:
                            in_phase_output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, G.edges[var1+':0', var2+':0']['prime_weight'], G.edges[var1+':1', var2+':0']['prime_weight'], G.edges[var1+':1', var2+':1']['prime_weight'],\
                                                                                         _00_cell_count, _01_cell_count, _11_cell_count, G.edges[var1+':0', var2+':0']['raw_read_count'], G.edges[var1+':1', var2+':0']['raw_read_count'], G.edges[var1+':1', var2+':1']['raw_read_count'],\
                                                                                            G.nodes[var1+':0']['raw_read_count'], G.nodes[var1+':1']['raw_read_count'], G.nodes[var2+':0']['raw_read_count'], G.nodes[var2+':1']['raw_read_count']))
                        except:
                            continue
                    elif (var1+':0', var2+':0') in G.edges and (var1+':0', var2+':1') in G.edges and (var1+':1', var2+':1') in G.edges and (var1+':1', var2+':0') not in G.edges:
                        _00_cell_count = len(G.edges[var1+':0', var2+':0']['barcodes'].keys())
                        _01_cell_count = len(G.edges[var1+':0', var2+':1']['barcodes'].keys())
                        _11_cell_count = len(G.edges[var1+':1', var2+':1']['barcodes'].keys())
                        try:
                            in_phase_output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, G.edges[var1+':0', var2+':0']['prime_weight'], G.edges[var1+':0', var2+':1']['prime_weight'], G.edges[var1+':1', var2+':1']['prime_weight'],\
                                                                                         _00_cell_count, _01_cell_count, _11_cell_count, G.edges[var1+':0', var2+':0']['raw_read_count'], G.edges[var1+':0', var2+':1']['raw_read_count'], G.edges[var1+':1', var2+':1']['raw_read_count'],\
                                                                                            G.nodes[var1+':0']['raw_read_count'], G.nodes[var1+':1']['raw_read_count'], G.nodes[var2+':0']['raw_read_count'], G.nodes[var2+':1']['raw_read_count']))
                        except:
                            continue

def examine_out_of_phase(lbd_list):
    for lbd in lbd_list:
        print('Processing {}...'.format(lbd))
        graph_path = './output/{}/conflict_graphs/'.format(lbd)
        if not os.path.exists(graph_path): return
        graphs = os.listdir(graph_path)
        for graph in graphs:
            try:
                G = pickle.load(open(graph_path+graph, 'rb'))
            except:
                continue
            variants = set()
            list(map(lambda x: variants.add(x.split(':')[0]), list(G.nodes)))
            variants = list(variants)
            for i in range(len(variants)-1):
                for j in range(i, len(variants)):
                    var1, var2 = variants[i], variants[j]
                    if (var1+':0', var2+':0') in G.edges and (var1+':1', var2+':0') in G.edges and (var1+':0', var2+':1') in G.edges and (var1+':1', var2+':1') not in G.edges:
                        _00_cell_count = len(G.edges[var1+':0', var2+':0']['barcodes'].keys())
                        _01_cell_count = len(G.edges[var1+':0', var2+':1']['barcodes'].keys())
                        _10_cell_count = len(G.edges[var1+':1', var2+':0']['barcodes'].keys())
                        try:
                            out_of_phase_output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, G.edges[var1+':0', var2+':0']['prime_weight'], G.edges[var1+':0', var2+':1']['prime_weight'], G.edges[var1+':1', var2+':0']['prime_weight'],\
                                                                                         _00_cell_count, _01_cell_count, _10_cell_count, G.edges[var1+':0', var2+':0']['raw_read_count'], G.edges[var1+':0', var2+':1']['raw_read_count'], G.edges[var1+':1', var2+':0']['raw_read_count'],\
                                                                                            G.nodes[var1+':0']['raw_read_count'], G.nodes[var1+':1']['raw_read_count'], G.nodes[var2+':0']['raw_read_count'], G.nodes[var2+':1']['raw_read_count']))
                        except:
                            continue

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--id", help="Input ID", required=True)
parser.add_argument("--gtf", help="Input GTF file", required=True)
opt = parser.parse_args()

examine_in_phase(opt.id)
examine_out_of_phase(opt.id)

import pandas as pd


def preprocess(output_file):

    potential_df = pd.read_csv('./{}.txt'.format(output_file), sep='\t', header=None, names=['Sample', 'File', 'Var1', 'Var2', '00', '01', '10', '00_cell', '01_cell', '10_cell', '00_raw_count', '01_raw_count', '10_raw_count', 'Var1:0', 'Var1:1', 'Var2:0', 'Var2:1'])

    potential_df = potential_df[['Sample', 'Var1', 'Var2', '00', '01', '10', '00_cell', '01_cell', '10_cell', '00_raw_count', '01_raw_count', '10_raw_count', 'Var1:0', 'Var1:1', 'Var2:0', 'Var2:1']]
    potential_df['Sample'] = potential_df['Sample'].apply(lambda x: x if not x.endswith('_002') else x[:x.index('_002')])
    potential_df = potential_df.groupby(['Sample', 'Var1', 'Var2'], as_index=False).agg(sum)

    print(potential_df.head())

    def calculate_inherit_ratio(group):
        _00_cell, _01_cell, _10_cell = group[6], group[7], group[8]
        if _10_cell > _01_cell:
            return _00_cell / _10_cell
        else:
            return _00_cell / _01_cell

    def calculate_somatic_ratio(group):
        _00_cell, _01_cell, _10_cell = group[6], group[7], group[8]
        if _10_cell > _01_cell:
            return round(_01_cell / _10_cell, 7)
        else:
            return round(_10_cell / _01_cell, 7)
        
    def calculate_minor_freq(group):
        _00_cell, _01_cell, _10_cell = group['00_cell'], group['01_cell'], group['10_cell']
        return min(_00_cell, _10_cell , _01_cell)
        

    def sum_reads(group):
        _00, _01, _10 = group['00'], group['01'], group['10']
        return sum([_00, _01, _10])

    def sum_counts(group):
        _00, _01, _10 = group['00_raw_count'], group['01_raw_count'], group['10_raw_count']
        return sum([_00, _01, _10])

    potential_df['inherit_ratio'] = potential_df.apply(calculate_inherit_ratio, axis=1)
    potential_df['somatic_ratio'] = potential_df.apply(calculate_somatic_ratio, axis=1)
    potential_df['minor_freq'] = potential_df.apply(calculate_minor_freq, axis=1)

    potential_df['reads_sum'] = potential_df.apply(sum_reads, axis=1)
    potential_df['raw_reads_count_sum'] = potential_df.apply(sum_counts, axis=1)

    potential_df.to_csv('{}.csv'.format(output_file), index=False)

def annotate(output_file):
    selected_pairs = pd.read_csv('{}.csv'.format(output_file), sep=',')
    # selected_pairs = pd.read_csv('rua.csv', sep=',')

    var1, var2 = selected_pairs['Var1'].tolist(), selected_pairs['Var2'].tolist()
    vars_ = var1 + var2
    vars_ = set(vars_)

    with open('sites.bed', 'w') as f:
        list(map(lambda x:f.write('{}\t{}\t{}\n'.format(x.split('_')[0], x.split('_')[1], int(x.split('_')[1])+1)), vars_))

    import subprocess
    cmd = 'bedtools intersect -a sites.bed -b {} -wa -wb | egrep "CDS\t" > annotated_sites.bed'.format(opt.gtf)
    subprocess.check_call(cmd, shell=True)

    import re
    gene_name_re = re.compile(r'gene_name "(.*?)";')
    f = open('annotated_sites.bed')

    pos_gene_map = {}

    for line in f.readlines():
        items = line.strip().split('\t')
        chr_, pos, gene = items[0], items[1], gene_name_re.findall(line)[0]
        pos_gene_map['_'.join([chr_, pos])] = gene

    genes = set(pos_gene_map.values())

    # list(map(lambda x: print(x, end=' '), genes))

    list(map(lambda x: genes.add(gene_name_re.findall(x)[0]), filter(lambda x: 'CDS\t' in x,f.readlines())))
    genes = list(filter(lambda x:'ENSG' not in x, genes))
    genes = list(filter(lambda x:'IG' not in x, genes))
    list(map(lambda x: print(x, end=' '), genes))

    def find_corresponding_gene(group):
        var1 = group['Var1'][:group['Var1'].index('_.')]
        if var1 not in pos_gene_map: return False
        return False if pos_gene_map[var1] not in genes else pos_gene_map[var1]

    selected_pairs['gene'] = selected_pairs.apply(find_corresponding_gene, axis=1)
    selected_pairs = selected_pairs[selected_pairs['gene']!=False]
    selected_pairs.to_csv('./{}.annotated.csv', sep='\t', index=False)

preprocess('in.phase.hits')
preprocess('out.of.phase.hits')
annotate('in.phase.hits')
annotate('out.of.phase.hits')