import networkx as nx
import matplotlib.pyplot as plt
import collections
import scipy.stats
import math

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
    
    return list(nx.connected_components(graph))

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

    final_haplotypes = de_duplicate(nonconflicted_nodes + resolved_conflicted_nodes)

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

    return total_hap, correct_hap, total_predict, correct_predict, total_nodes


def report_singular_cells(opt, removed_edges:list[nx.Graph], mean, var, n):
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
            node_list = sorted(node_list, key=lambda x: x[0])
            for i in range(0, len(node_list) - 1):
                for j in range(1, len(node_list)):
                    (left_node, left_count), (right_node, right_count) = node_list[i], node_list[j]
                    oppo_left_node, oppo_right_node = get_opposite_allele(left_node), get_opposite_allele(right_node)
                    link_value = min(left_count, right_count)
                    if (left_node, right_node) in barcode_pair_map[barcode].keys():
                        barcode_pair_map[barcode][(left_node, right_node)] += link_value
                    elif (oppo_left_node, oppo_right_node) in barcode_pair_map[barcode].keys():
                        barcode_pair_map[barcode][(oppo_left_node, oppo_right_node)] += link_value
                    else:
                        barcode_pair_map[barcode][(left_node, right_node)] = 0

    # output
    df = n-1
    T = scipy.stats.t(df)
    with open('singular_cell_linkage.txt', 'w') as f:
        f.write('barcode\tvar\tgeno\tp-value\n')
    for barcode, allele_pairs in barcode_pair_map.items():
        for allele_pair, link_value in allele_pairs:
            p_value = T.cdf((mean-link_value)/math.sqrt(var/10))
            f.write('{}\t{}\t{}\t{}\n'.format(barcode, ','.join(map(lambda x:x.split(':')[0], allele_pair)), ''.join(map(lambda x:x.split(':')[1], allele_pair)), p_value))