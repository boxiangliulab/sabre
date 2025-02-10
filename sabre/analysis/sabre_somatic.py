# python ~/scratch/sabre/Faser-scRNA-main/sabre/analysis/sabre_somatic.py --id liu-2024 --threads 47 --gtf ~/scratch/refdata-gex-GRCh38-2024-A/genes/genes.gtf

import networkx as nx
import multiprocessing.managers
import os
import pickle

import argparse
import pandas as pd

import time
from tqdm import tqdm

import copy

import multiprocessing

def update_progress_bars(progress_list, total_iterations, n_processes):
    bars = [tqdm(total=total_iterations+1, position=i) for i in range(n_processes)]
    # return
    while True:
        for i, bar in enumerate(bars):
            bar.n = progress_list[i]  # 更新tqdm的进度
            bar.refresh()  # 刷新显示
        # if sum(progress_list) >= (total_iterations-1) * n_processes:
        #     break
        time.sleep(0.1)

def check_in_phase_graph(graphs, graph_path, output_list, rawlink_list, lbd, progress_list, thread_id):
    total_graphs_nums = len(graphs)
    
    for graph in graphs:
        try:
            G = pickle.load(open(graph_path+graph, 'rb'))

            variants = set()
            list(map(lambda x: variants.add(x.split(':')[0]), list(G.nodes)))
            variants = list(variants)
            for i in range(len(variants)-1):
                for j in range(i+1, len(variants)):
                    var1, var2 = variants[i], variants[j]

                    node_1_0, node_1_1, node_2_0, node_2_1 = var1+':0', var1+':1', var2+':0', var2+':1'

                    if not (node_1_0 in G and node_1_1 in G and node_2_0 in G and node_2_1 in G): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_1)
                    G_prime.remove_node(node_2_1)
                    if not nx.has_path(G_prime, node_1_0, node_2_0): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_0)
                    G_prime.remove_node(node_2_0)
                    if not nx.has_path(G_prime, node_1_1, node_2_1): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_1)
                    G_prime.remove_node(node_2_0)
                    if nx.has_path(G_prime, node_1_0, node_2_1): 
                        G_prime = nx.Graph(G)
                        G_prime.remove_node(node_1_0)
                        G_prime.remove_node(node_2_1)
                        if nx.has_path(G_prime, node_1_1, node_2_0):
                            continue

                        _00_cell_count = len(G.edges[node_1_0, node_2_0]['barcodes'].keys()) if (node_1_0, node_2_0) in G.edges else -1
                        _01_cell_count = len(G.edges[node_1_0, node_2_1]['barcodes'].keys()) if (node_1_0, node_2_1) in G.edges else -1
                        _11_cell_count = len(G.edges[node_1_1, node_2_1]['barcodes'].keys()) if (node_1_1, node_2_1) in G.edges else -1

                        _00_cell_prime_weight = G.edges[node_1_0, node_2_0]['prime_weight'] if (node_1_0, node_2_0) in G.edges else -1
                        _01_cell_prime_weight = G.edges[node_1_0, node_2_1]['prime_weight'] if (node_1_0, node_2_1) in G.edges else -1
                        _11_cell_prime_weight = G.edges[node_1_1, node_2_1]['prime_weight'] if (node_1_1, node_2_1) in G.edges else -1

                        _00_cell_raw_count = G.edges[node_1_0, node_2_0]['raw_read_count'] if (node_1_0, node_2_0) in G.edges else -1
                        _01_cell_raw_count = G.edges[node_1_0, node_2_1]['raw_read_count'] if (node_1_0, node_2_1) in G.edges else -1
                        _11_cell_raw_count = G.edges[node_1_1, node_2_1]['raw_read_count'] if (node_1_1, node_2_1) in G.edges else -1

                        if rawlink_list is not None:
                            if (node_1_0, node_2_0) in G.edges:
                                for barcode in G.edges[node_1_0, node_2_0]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_0, node_2_0))
                            if (node_1_0, node_2_1) in G.edges:
                                for barcode in G.edges[node_1_0, node_2_1]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_0, node_2_1))
                            if (node_1_1, node_2_1) in G.edges:
                                for barcode in G.edges[node_1_1, node_2_1]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_1, node_2_1))


                        output_list.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, _00_cell_prime_weight, _01_cell_prime_weight, _11_cell_prime_weight,\
                                                                                        _00_cell_count, _01_cell_count, _11_cell_count, _00_cell_raw_count, _01_cell_raw_count, _11_cell_raw_count, \
                                                                                        G.nodes[var1+':0']['allele_read_count'], G.nodes[var1+':1']['allele_read_count'], G.nodes[var2+':0']['allele_read_count'], G.nodes[var2+':1']['allele_read_count']))
                        continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_0)
                    G_prime.remove_node(node_2_1)
                    if nx.has_path(G_prime, node_1_1, node_2_0): 
                        G_prime = nx.Graph(G)
                        G_prime.remove_node(node_1_1)
                        G_prime.remove_node(node_2_0)
                        if nx.has_path(G_prime, node_1_0, node_2_1):
                            continue

                        _00_cell_count = len(G.edges[node_1_0, node_2_0]['barcodes'].keys()) if (node_1_0, node_2_0) in G.edges else -1
                        _10_cell_count = len(G.edges[node_1_1, node_2_0]['barcodes'].keys()) if (node_1_1, node_2_0) in G.edges else -1
                        _11_cell_count = len(G.edges[node_1_1, node_2_1]['barcodes'].keys()) if (node_1_1, node_2_1) in G.edges else -1

                        _00_cell_prime_weight = G.edges[node_1_0, node_2_0]['prime_weight'] if (node_1_0, node_2_0) in G.edges else -1
                        _10_cell_prime_weight = G.edges[node_1_1, node_2_0]['prime_weight'] if (node_1_1, node_2_0) in G.edges else -1
                        _11_cell_prime_weight = G.edges[node_1_1, node_2_1]['prime_weight'] if (node_1_1, node_2_1) in G.edges else -1

                        _00_cell_raw_count = G.edges[node_1_0, node_2_0]['raw_read_count'] if (node_1_0, node_2_0) in G.edges else -1
                        _10_cell_raw_count = G.edges[node_1_1, node_2_0]['raw_read_count'] if (node_1_1, node_2_0) in G.edges else -1
                        _11_cell_raw_count = G.edges[node_1_1, node_2_1]['raw_read_count'] if (node_1_1, node_2_1) in G.edges else -1

                        if rawlink_list is not None:
                            if (node_1_0, node_2_0) in G.edges:
                                for barcode in G.edges[node_1_0, node_2_0]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_0, node_2_0))
                            if (node_1_1, node_2_0) in G.edges:
                                for barcode in G.edges[node_1_1, node_2_0]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_1, node_2_0))
                            if (node_1_1, node_2_1) in G.edges:
                                for barcode in G.edges[node_1_1, node_2_1]['barcodes'].keys():
                                    rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_1, node_2_1))

                        output_list.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, _00_cell_prime_weight, _10_cell_prime_weight, _11_cell_prime_weight,\
                                                                                        _00_cell_count, _10_cell_count, _11_cell_count, _00_cell_raw_count, _10_cell_raw_count, _11_cell_raw_count, \
                                                                                        G.nodes[var1+':0']['allele_read_count'], G.nodes[var1+':1']['allele_read_count'], G.nodes[var2+':0']['allele_read_count'], G.nodes[var2+':1']['allele_read_count']))

        except Exception as e:
            continue
                
        finally:
            progress_list[thread_id] += 1


def examine_in_phase(opt, lbd):
    in_phase_output = open('./{}.in.phase.hits.txt'.format(lbd), 'w')
    print('Processing {} in phase...'.format(lbd))
    graph_path = '{}/{}/conflict_graphs/'.format(opt.output_dir, lbd)
    if not os.path.exists(graph_path): return
    graphs = os.listdir(graph_path)
    total_iterations = int(len(graphs)/opt.threads) + 1

    with multiprocessing.Manager() as manager:
        output_list = manager.list()
        rawlink_list = manager.list() if opt.rawlink else None
        progress_list = manager.list(([0] * (opt.threads-1)))
        
        processes = []
        for i in range(opt.threads):
            p = multiprocessing.Process(target=check_in_phase_graph, args=(graphs[i::opt.threads-1], graph_path, output_list, rawlink_list, lbd, progress_list, i)) if i != opt.threads-1 else multiprocessing.Process(target=update_progress_bars, args=(progress_list, total_iterations, opt.threads-1))
            p.deamon = True if i == opt.threads-1 else False
            processes.append(p)
            p.start()

        for p in processes[:-1]:
            p.join()
        
        processes[-1].kill()

        list(map(in_phase_output.write, set(output_list)))
        if rawlink_list is not None:
            with open('./{}.in.phase.hits.raw.linkage.txt'.format(lbd), 'w') as f:
                list(map(f.write, rawlink_list))

    in_phase_output.close()


def check_out_of_phase_graph(graphs, graph_path, output_list, rawlink_list, lbd, progress_list, thread_id):
    total_graphs_nums = len(graphs)
    
    for graph in graphs:
        try:
            G = pickle.load(open(graph_path+graph, 'rb'))

            variants = set()
            list(map(lambda x: variants.add(x.split(':')[0]), list(G.nodes)))
            variants = list(variants)
            for i in range(len(variants)-1):
                for j in range(i+1, len(variants)):
                    var1, var2 = variants[i], variants[j]

                    node_1_0, node_1_1, node_2_0, node_2_1 = var1+':0', var1+':1', var2+':0', var2+':1'

                    if not (node_1_0 in G and node_1_1 in G and node_2_0 in G and node_2_1 in G): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_1)
                    G_prime.remove_node(node_2_0)
                    if not nx.has_path(G_prime, node_1_0, node_2_1): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_0)
                    G_prime.remove_node(node_2_1)
                    if not nx.has_path(G_prime, node_1_1, node_2_0): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_1)
                    G_prime.remove_node(node_2_1)
                    if not nx.has_path(G_prime, node_1_0, node_2_0): continue

                    G_prime = nx.Graph(G)
                    G_prime.remove_node(node_1_0)
                    G_prime.remove_node(node_2_0)
                    if nx.has_path(G_prime, node_1_1, node_2_1): continue

                    _00_cell_count = len(G.edges[node_1_0, node_2_0]['barcodes'].keys()) if (node_1_0, node_2_0) in G.edges else -1
                    _01_cell_count = len(G.edges[node_1_0, node_2_1]['barcodes'].keys()) if (node_1_0, node_2_1) in G.edges else -1
                    _10_cell_count = len(G.edges[node_1_1, node_2_0]['barcodes'].keys()) if (node_1_1, node_2_0) in G.edges else -1

                    _00_cell_prime_weight = G.edges[node_1_0, node_2_0]['prime_weight'] if (node_1_0, node_2_0) in G.edges else -1
                    _01_cell_prime_weight = G.edges[node_1_0, node_2_1]['prime_weight'] if (node_1_0, node_2_1) in G.edges else -1
                    _10_cell_prime_weight = G.edges[node_1_1, node_2_0]['prime_weight'] if (node_1_1, node_2_0) in G.edges else -1

                    _00_cell_raw_count = G.edges[node_1_0, node_2_0]['raw_read_count'] if (node_1_0, node_2_0) in G.edges else -1
                    _01_cell_raw_count = G.edges[node_1_0, node_2_1]['raw_read_count'] if (node_1_0, node_2_1) in G.edges else -1
                    _10_cell_raw_count = G.edges[node_1_1, node_2_0]['raw_read_count'] if (node_1_1, node_2_0) in G.edges else -1

                    if rawlink_list is not None:
                        if (node_1_0, node_2_0) in G.edges:
                            for barcode in G.edges[node_1_0, node_2_0]['barcodes'].keys():
                                rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_0, node_2_0))
                        if (node_1_1, node_2_0) in G.edges:
                            for barcode in G.edges[node_1_1, node_2_0]['barcodes'].keys():
                                rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_1, node_2_0))
                        if (node_1_0, node_2_1) in G.edges:
                            for barcode in G.edges[node_1_0, node_2_1]['barcodes'].keys():
                                rawlink_list.append('{}\t{}\t{}\n'.format(barcode, node_1_0, node_2_1))

                    output_list.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(lbd, graph, var1, var2, _00_cell_prime_weight, _01_cell_prime_weight, _10_cell_prime_weight,\
                                                                                    _00_cell_count, _01_cell_count, _10_cell_count, _00_cell_raw_count, _01_cell_raw_count, _10_cell_raw_count, \
                                                                                    G.nodes[var1+':0']['allele_read_count'], G.nodes[var1+':1']['allele_read_count'], G.nodes[var2+':0']['allele_read_count'], G.nodes[var2+':1']['allele_read_count']))

        except Exception as e:
            print(e)
            continue
                
        finally:
            progress_list[thread_id] += 1

def examine_out_of_phase(opt, lbd):
    out_of_phase_output = open('./{}.out.of.phase.hits.txt'.format(lbd), 'w')
    print('Processing {} out-of-phase...'.format(lbd))
    graph_path = '{}/{}/conflict_graphs/'.format(opt.output_dir, lbd)
    if not os.path.exists(graph_path): return
    graphs = os.listdir(graph_path)
    total_iterations = int(len(graphs)/(opt.threads-1)) + 1

    with multiprocessing.Manager() as manager:
        output_list = manager.list()
        rawlink_list = manager.list() if opt.rawlink else None
        progress_list = manager.list(([0] * (opt.threads-1)))
        
        processes = []
        for i in range(opt.threads):
            p = multiprocessing.Process(target=check_out_of_phase_graph, args=(graphs[i::opt.threads-1], graph_path, output_list, rawlink_list, lbd, progress_list, i)) if i != opt.threads-1 else multiprocessing.Process(target=update_progress_bars, args=(progress_list, total_iterations, opt.threads-1))
            p.deamon = True if i == opt.threads-1 else False
            processes.append(p)
            p.start()

        for p in processes[:-1]:
            p.join()
        
        print('All processes Finished')
        processes[-1].kill()

        list(map(out_of_phase_output.write, set(output_list)))
        if rawlink_list is not None:
            with open('./{}.out.of.phase.hits.raw.linkage.txt'.format(lbd), 'w') as f:
                list(map(f.write, rawlink_list))

    out_of_phase_output.close()



def preprocess(output_file):
    print('Preprocessing {}...'.format(output_file))
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

    potential_df.to_csv('{}.csv'.format(output_file), index=False)

def annotate(opt, output_file, gtf):
    print('Annotating {}...'.format(output_file))
    selected_pairs = pd.read_csv('{}.csv'.format(output_file), sep=',')
    # selected_pairs = pd.read_csv('rua.csv', sep=',')

    var1, var2 = selected_pairs['Var1'].tolist(), selected_pairs['Var2'].tolist()
    vars_ = var1 + var2
    vars_ = set(vars_)

    with open('sites.bed', 'w') as f:
        list(map(lambda x:f.write('{}\t{}\t{}\n'.format(x.split('_')[0], x.split('_')[1], int(x.split('_')[1])+1)), vars_))

    import subprocess

    if gtf.endswith('.gtf.gz'):
        cmd = f'zcat {gtf} > {gtf[gtf.index('.gz')]}'
        gtf_for_bedtools = gtf[gtf.index('.gz')]
    else:
        gtf_for_bedtools = gtf

    if opt.cds:
        cmd = 'bedtools intersect -a sites.bed -b {} -wa -wb > annotated_sites.bed | egrep "CDS"'.format(gtf_for_bedtools)
    else:
        cmd = 'bedtools intersect -a sites.bed -b {} -wa -wb > annotated_sites.bed'.format(gtf_for_bedtools)
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


    genes = list(filter(lambda x:'ENSG' not in x, genes))
    list(map(lambda x: print(x, end=' '), genes))

    with open('./{}.genes.list'.format(output_file), 'w') as f:
        list(map(lambda x: f.write('{} '.format(x)), genes))

    def find_corresponding_gene(group):
        var1 = group['Var1'][:group['Var1'].index('_.')]
        if var1 not in pos_gene_map: return False
        return False if pos_gene_map[var1] not in genes else pos_gene_map[var1]

    selected_pairs['gene'] = selected_pairs.apply(find_corresponding_gene, axis=1)
    selected_pairs = selected_pairs[selected_pairs['gene']!=False]
    selected_pairs.to_csv('./{}.annotated.csv'.format(output_file), sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", help="Input ID", required=True)
    parser.add_argument("--output_dir", help="Path to output directory, should be the same as sabre --output_dir. Default ./output", required=False, default='./output')
    parser.add_argument("--threads", help="Multithread number", type=int, default=1)
    parser.add_argument("--rawlink", help="Output raw linkage per cell per variant pair", action='store_true')
    parser.add_argument("--gtf", help="Input GTF file. We recommend you to input a .gtf file rather than a .gtf.gz file, because sabre-somatic will have to depress the .gtf.gz file everytime you input a compressed gtf file.", required=True)
    parser.add_argument("--cds", help="If specified, only mutations on CDS will be phased.", action='store_true')
    opt = parser.parse_args()

    # Check Arguments
    if not os.path.exists(os.path.join(opt.output_dir, opt.id)):
        raise ValueError("File not found: {}. Please check input argument --id and --output_dir".format(opt.output_dir, opt.id))
    
    if opt.threads < 1:
        raise ValueError("Thread number cannot be less than 1!")
    
    if not opt.gtf.endswith('.gtf'):
        raise ValueError("Invalid argument for --gtf. Please give a valid .gtf file or .gtf.gz file.")
    
    if not os.path.exists(opt.gtf):
        raise ValueError("File not found: {}. Please check input argument --gtf".format(opt.gtf))

    examine_in_phase(opt, opt.id)
    preprocess('{}.in.phase.hits'.format(opt.id))
    annotate(opt, '{}.in.phase.hits'.format(opt.id), opt.gtf)

    examine_out_of_phase(opt, opt.id)
    preprocess('{}.out.of.phase.hits'.format(opt.id))
    annotate(opt, '{}.out.of.phase.hits'.format(opt.id), opt.gtf)

if __name__ == '__main__':
    main()