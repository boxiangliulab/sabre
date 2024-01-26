import networkx as nx
import matplotlib.pyplot as plt
import collections
import scipy.stats
import math
import os

def get_opposite_allele(allele, is_opposite=True):
    '''
    Input: A:0
    Output: A:1
    '''
    if is_opposite:
        return allele.split(':')[0]+':{}'.format(1-int(allele.split(':')[1]))
    return allele

def check_haplotype(ground_truth, predict):
    '''
    #1 
    Input : 
        ground_truth: A0B0C0
        predict: A1B1C1
    Output:
        True

    #2
    Input : 
        ground_truth: A0B0C1
        predict: A1B1C1
    Output:
        False
    '''

    is_opposite = False
    is_first_allele = True
    for i, gt in enumerate(ground_truth):
        if gt == '-':
            continue
        if is_first_allele:
            is_first_allele = False
            if gt != predict[i]:
                is_opposite = True
                continue
        else:
            if (not is_opposite and gt != predict[i]) or (is_opposite and gt == predict[i]):
                return False
    return True

def de_duplicate(alleles_list:list[list[str]]):
    '''
    As many of these haplotypes are actually not disjointed from each other, like
        1---2
    0---1   
    it should be 
    0---1---2
    Thus, we solve this problem by graph.
    '''

    graph = nx.Graph()
    # alleles_list = sorted(map(sorted, alleles_list))
    for alleles in alleles_list:
        is_opposite = False
        for allele in alleles:
            opposite_allele = get_opposite_allele(allele)
            if opposite_allele in graph.nodes:
                is_opposite = True

        for i in range(len(alleles)-1):
            graph.add_edge(get_opposite_allele(alleles[i], is_opposite), get_opposite_allele(alleles[i+1], is_opposite))
    
    return graph

def draw_graph_with_weights(folder, G:nx.Graph):
    '''
    Demonstrate failed graphs
    '''
    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
    nx.draw_networkx_nodes(G, pos, node_size=10)
    nx.draw_networkx_edges(G, pos, width=1)

    nx.draw_networkx_labels(G, pos, font_size=5, font_family="sans-serif")
    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos, edge_labels)
    plt.savefig('./output/{}_graph_pdfs/{}.pdf'.format(folder, list(G.nodes)[0]))
    plt.clf()

def report_phasing_result(opt, G, nonconflicted_nodes, resolved_conflicted_nodes, vid_var_map):
    '''
    Report the phasing result of Faser.
    Nonconflicted nodes and predicted conflicted nodes are reported separately.
    '''
    
    total_hap = 0
    correct_hap = 0
    total_predict = 0
    correct_predict = 0
    already_reported_set = []
    total_nodes = 0

    final_graph = de_duplicate(nonconflicted_nodes + resolved_conflicted_nodes)
    final_haplotypes = list(nx.connected_components(final_graph))

    if not os.path.exists('./output'):
        os.mkdir('./output')

    with open('./output/chr_{}_haplotypes.tsv'.format(opt.restrict_chr), 'w') as f:
        f.write('haplotype\tpredicted_phasing\tgt_phasing\tcorrection\n')
        for nodes in final_haplotypes:
            var_phasing_list = list(map(lambda x: x.split(':'), nodes))
            total_nodes += len(var_phasing_list)
            vids, phasing = zip(*var_phasing_list)
            var_set = set(vids)
            if len(var_set) <= 1:
                continue
            if var_set in already_reported_set:
                continue
            already_reported_set.append(var_set)
            vids_string = ','.join(vids)
            predicted_phasing_str = ''.join(phasing)

            #Get ground_truth phasing result
            gt_phasing_str = ''.join(list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', vids)))
            is_correct = 0
            if check_haplotype(gt_phasing_str, predicted_phasing_str):
                is_correct = 1

            total_hap += 1
            correct_hap += is_correct

            f.write('{}\t{}\t{}\t{}\n'.format(vids_string, predicted_phasing_str, gt_phasing_str, is_correct))

    return total_hap, correct_hap, total_predict, correct_predict, total_nodes, final_graph



def deprecated_report_singular_cells(opt, removed_edges:list[nx.Graph], final_graph:list[nx.Graph], mean, var, n):
    '''
    Report removed edges during resolving confliced graphs.
    '''
    barcode_pair_map = collections.defaultdict(dict)
    for sg in removed_edges:
        barcode_node_map = collections.defaultdict(list)
        for node in list(sg.nodes):
            for barcode, count in sg.nodes[node]['cells']:
                barcode_node_map[barcode].append((node, count))
        for barcode, node_list in barcode_node_map.items():
            if len(node_list) == 1:
                continue
            node_list = sorted(node_list, key=lambda x: int(x[0].split('_')[1]))
            for i in range(0, len(node_list) - 1):
                for j in range(i+1, len(node_list)):
                    (left_node, left_count), (right_node, right_count) = node_list[i], node_list[j]
                    if(left_node, right_node) not in sg.edges:
                        continue
                    if left_node.split(':')[0] == right_node.split(':')[0]:
                        continue
                    oppo_left_node, oppo_right_node = get_opposite_allele(left_node), get_opposite_allele(right_node)
                    if left_node in final_graph.nodes and right_node in final_graph.nodes:
                        if nx.has_path(final_graph, left_node, right_node):
                            continue
                    if oppo_left_node in final_graph.nodes and oppo_right_node in final_graph.nodes:
                        if nx.has_path(final_graph, oppo_left_node, oppo_right_node):
                            continue                
                    link_value = min(left_count, right_count)
                    if (left_node, right_node) in barcode_pair_map[barcode].keys():
                        barcode_pair_map[barcode][(left_node, right_node)] += link_value
                    elif (oppo_left_node, oppo_right_node) in barcode_pair_map[barcode].keys():
                        barcode_pair_map[barcode][(oppo_left_node, oppo_right_node)] += link_value
                    else:
                        barcode_pair_map[barcode][(left_node, right_node)] = link_value

    # output
    df = n-1
    T = scipy.stats.t(df)
    all_pairs, selected_pairs = 0, 0
    with open('singular_cell_linkage.txt', 'w') as f:
        f.write('barcode\tvar\tgeno\tsupport\tp-value\n')
        for barcode, allele_pairs in barcode_pair_map.items():
            for allele_pair, link_value in allele_pairs.items():
                p_value = T.cdf((mean-link_value)/math.sqrt(var/10))
                if p_value > 0.05:
                    continue
                f.write('{}\t{}\t{}\t{}\t{}\n'.format(barcode, ','.join(map(lambda x:x.split(':')[0], allele_pair)), ''.join(map(lambda x:x.split(':')[1], allele_pair)),link_value, p_value))
    
    print()

def report_singular_cells(opt, removed_sub_graphs:list[nx.Graph], final_graph:nx.Graph, allele_linkage_graph:nx.Graph, vid_var_map, mean, var, n):
    '''
    Report removed edges during resolving confliced graphs.
    '''

    barcode_pair_map = collections.defaultdict(dict)
    barcode_pair_neg_map = collections.defaultdict(dict)
    for sg in removed_sub_graphs:
        calculated_pairs = set()
        for l, r in list(sg.edges):
            left_node, right_node = max(l, r), min(l, r)
            oppo_left_node, oppo_right_node = get_opposite_allele(left_node), get_opposite_allele(right_node)
            if (left_node, right_node) in calculated_pairs:
                continue
            calculated_pairs.add((left_node, right_node))

            if left_node in final_graph.nodes and right_node in final_graph.nodes:
                if nx.has_path(final_graph, left_node, right_node):
                    continue
            if oppo_left_node in final_graph.nodes and oppo_right_node in final_graph.nodes:
                if nx.has_path(final_graph, oppo_left_node, oppo_right_node):
                    continue

            if (oppo_left_node, right_node) in allele_linkage_graph.edges:
                barcode_weight_map:dict = allele_linkage_graph.edges[(oppo_left_node, right_node)]['barcodes']
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    if (left_node, right_node) not in barcode_pair_neg_map[barcode].keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
                    if barcode in barcode_weight_map.keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] += barcode_weight_map[barcode]
            elif (left_node, oppo_right_node) in allele_linkage_graph.edges:
                barcode_weight_map:dict = allele_linkage_graph.edges[(left_node, oppo_right_node)]['barcodes']
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    if (left_node, right_node) not in barcode_pair_neg_map[barcode].keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
                    if barcode in barcode_weight_map.keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] += barcode_weight_map[barcode]
            else:
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
            
            for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                if (left_node, right_node) not in barcode_pair_map[barcode].keys():
                    barcode_pair_map[barcode][(left_node, right_node)] = count
                else:
                    barcode_pair_map[barcode][(left_node, right_node)] += count
    
    # output
    df = n-1
    T = scipy.stats.t(df)
    all_pairs, selected_pairs = 0, 0
    if not os.path.exists('./output/singular_cells'):
        os.mkdir('./output/singular_cells')
    with open('./output/singular_cells/singular_cell_linkage_{}.txt'.format(opt.restrict_chr), 'w') as f:
        f.write('barcode\tvar\tgeno\tsupport\tp-value\toppo_support\tcorrect\n')
        for barcode, allele_pairs in barcode_pair_map.items():
            for allele_pair, link_value in allele_pairs.items():
                p_value = T.cdf((mean-link_value)/math.sqrt(var/10))
                if p_value > 0.05:
                    continue
                gt_phasing_str = ''.join(list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', map(lambda x:x.split(':')[0], allele_pair))))
                is_correct = 0
                if check_haplotype(gt_phasing_str, ''.join(map(lambda x:x.split(':')[1], allele_pair))):
                    is_correct = 1
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(barcode, ','.join(map(lambda x:x.split(':')[0], allele_pair)), ''.join(map(lambda x:x.split(':')[1], allele_pair)),link_value, p_value, barcode_pair_neg_map[barcode][allele_pair], is_correct))
