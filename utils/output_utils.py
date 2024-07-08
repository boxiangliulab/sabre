import networkx as nx
import matplotlib.pyplot as plt
import collections
import scipy.stats
import math
import os
from utils.graph_utils import find_conflict_alleles

def get_opposite_allele(allele, is_opposite=True):
    '''
    Input: A:0
    Output: A:1
    '''
    if is_opposite:
        return allele.split(':')[0] + ':{}'.format(0 if allele.split(':')[1] != '0' else 1)
    return allele

def check_haplotype(ground_truth, predict):
    '''
    #1 
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
    predict = predict.replace('2', '1')
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

def de_duplicate(G:nx.Graph, alleles_list:list[list[str]], flag):
    '''
    As many of these haplotypes are actually not disjointed from each other, like
        1---2
    0---1   
    it should be 
    0---1---2
    Thus, we solve this problem by graph.
    '''

    # alleles_list = sorted(map(sorted, alleles_list))
    for alleles in alleles_list:
        is_opposite = False
        is_opposite_list = []
        for allele in alleles:
            opposite_allele = get_opposite_allele(allele)
            if opposite_allele in G.nodes:
                is_opposite_list.append(True)
            elif allele in G.nodes:
                is_opposite_list.append(False)
        if not (all(is_opposite_list) or all(map(lambda x : not(x), is_opposite_list))):
            continue
        is_opposite = all(is_opposite_list)
        for i in range(len(alleles)-1):
            G.add_edge(get_opposite_allele(alleles[i], is_opposite), get_opposite_allele(alleles[i+1], is_opposite), flag=flag)
    
    return G

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
    total_var_set = set()

    predict_pairs = 0
    correct_pairs = 0

    genome_coverage = []
    correct_variants = 0

    final_graph = nx.Graph()
    final_graph = de_duplicate(final_graph, nonconflicted_nodes, 'non')
    final_graph = de_duplicate(final_graph, resolved_conflicted_nodes, 'con')

    final_haplotypes = []
    final_haplotypes = list(nx.connected_components(final_graph))

    if not os.path.exists('./output/{}'.format(opt.id)):
        os.mkdir('./output/{}'.format(opt.id))

    with open('./output/{}/{}.output.SABRE'.format(opt.id, opt.chr), 'w') as g:
        for nodes in final_haplotypes:
            # var_phasing_list: [(chr1_147242777_._G_C, 1), ...]
            var_phasing_list = list(map(lambda x: x.split(':'), nodes))
            var_phasing_list = sorted(var_phasing_list, key=lambda x: int(x[0].split(opt.sep)[1]))
            vids, phasing = zip(*var_phasing_list)
            var_set = set(vids)
            if len(var_set) <= 1:
                continue
            if var_set in already_reported_set:
                continue
            already_reported_set.append(var_set)

            splited_haplotypes = [[var_phasing_list[0]]]

            for pair in var_phasing_list[1:]:
                pos_pre = int(splited_haplotypes[-1][-1][0].split(opt.sep)[1])
                pos_now = int(pair[0].split(opt.sep)[1])
                if pos_now - pos_pre >= opt.interval_threshold:
                    splited_haplotypes.append([pair])
                else:
                    splited_haplotypes[-1].append(pair)

            idx = 0
            for hap in splited_haplotypes:

                if len(hap) <= 1: continue

                vids, phasing = list(zip(*hap))

                total_var_set = total_var_set | set(vids)
                
                vids_string = ','.join(vids)
                predicted_phasing_str = ''.join(phasing)

                #Get ground_truth phasing result
                gt_phasing_str = ''.join(list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', vids)))
                is_correct = 0
                if check_haplotype(gt_phasing_str, predicted_phasing_str):
                    is_correct = 1

                total_hap += 1
                correct_hap += is_correct

                predict_pairs += math.comb(len(predicted_phasing_str), 2)
                if is_correct:
                    correct_variants += len(predicted_phasing_str)
                    correct_pairs += math.comb(len(predicted_phasing_str), 2)
                else:
                    cnt = 0
                    for al1, al2 in zip(gt_phasing_str, predicted_phasing_str):
                        if al1 != al2:
                            cnt += 1
                    correct_pairs += math.comb(len(predicted_phasing_str)-cnt, 2)
                    correct_pairs += math.comb(cnt, 2)
                    correct_variants += max(len(predicted_phasing_str)-cnt, cnt)

                genome_poses = list(map(lambda x: int(x.split(opt.sep)[1]), vids))
                genome_coverage.append(max(genome_poses) - min(genome_poses))

                g.write('BLOCK: offset: {} len: {} phased: {} SPAN: {} correct: {}\n'.format(vid_var_map[vids[0]].end, len(vids), len(vids), vid_var_map[vids[-1]].end - vid_var_map[vids[0]].end, is_correct))
                for vid, p_, g_ in zip(vids, phasing, list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', vids))):
                    g.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(idx, p_, abs(1-int(p_)), opt.chr, vid_var_map[vid].end, vid_var_map[vid].unique_id.split(opt.sep)[-2 + int(p_)], vid_var_map[vid].unique_id.split(opt.sep)[-2 + abs(1-int(p_))], '{}|{}'.format(g_, abs(1-int(g_))) if g_ != '-' else '-|-'))
                    idx += 1
                g.write('****************\n')

    return total_hap, correct_hap, total_predict, correct_predict, len(total_var_set), final_graph, predict_pairs, correct_pairs, correct_variants, genome_coverage 



def report_singular_cells(opt, removed_sub_graphs:list[nx.Graph], final_graph:nx.Graph, allele_linkage_graph:nx.Graph, vid_var_map, mean, var, n):
    '''
    Report removed edges during resolving confliced graphs.
    '''

    barcode_pair_map = collections.defaultdict(dict)
    barcode_pair_global_neg_map = collections.defaultdict(int)
    barcode_pair_neg_map = collections.defaultdict(dict)
    pair_neg_count_map = collections.defaultdict(int)
    pair_singular_map = {}
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

            pair_singular_map[(left_node, right_node)] = 0
            if oppo_left_node in final_graph.nodes and right_node in final_graph.nodes:
                if nx.has_path(final_graph, oppo_left_node, right_node):
                    pair_singular_map[(left_node, right_node)] = 1
            elif left_node in final_graph.nodes and oppo_right_node in final_graph.nodes:
                if nx.has_path(final_graph, left_node, oppo_right_node):
                    pair_singular_map[(left_node, right_node)] = 1

            if (oppo_left_node, right_node) in allele_linkage_graph.edges:
                barcode_weight_map:dict = allele_linkage_graph.edges[(oppo_left_node, right_node)]['barcodes']
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    if (left_node, right_node) not in barcode_pair_neg_map[barcode].keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
                    if barcode in barcode_weight_map.keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] += barcode_weight_map[barcode]
                if (left_node, right_node) not in pair_neg_count_map.keys():
                    pair_neg_count_map[(left_node, right_node)] += allele_linkage_graph.edges[(oppo_left_node, right_node)]['prime_weight']
            if (left_node, oppo_right_node) in allele_linkage_graph.edges:
                barcode_weight_map:dict = allele_linkage_graph.edges[(left_node, oppo_right_node)]['barcodes']
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    if (left_node, right_node) not in barcode_pair_neg_map[barcode].keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
                    if barcode in barcode_weight_map.keys():
                        barcode_pair_neg_map[barcode][(left_node, right_node)] += barcode_weight_map[barcode]
                pair_neg_count_map[(left_node, right_node)] += allele_linkage_graph.edges[(left_node, oppo_right_node)]['prime_weight']
            elif (oppo_left_node, right_node) not in allele_linkage_graph.edges:
                for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                    barcode_pair_neg_map[barcode][(left_node, right_node)] = 0
            
            for barcode, count in sg.edges[(left_node, right_node)]['barcodes'].items():
                if (left_node, right_node) not in barcode_pair_map[barcode].keys():
                    barcode_pair_map[barcode][(left_node, right_node)] = count
                else:
                    barcode_pair_map[barcode][(left_node, right_node)] += count

            if (left_node, right_node) not in barcode_pair_global_neg_map.keys():
                barcodes = set()
                assert (left_node, right_node) in allele_linkage_graph.edges
                if (left_node, right_node) in allele_linkage_graph.edges:
                    barcode_weight_map = allele_linkage_graph.edges[(left_node, right_node)]['barcodes']
                    for barcode in barcode_weight_map.keys():
                        barcodes.add(barcode)
                if (oppo_left_node, right_node) in allele_linkage_graph.edges:
                    barcode_weight_map = allele_linkage_graph.edges[(oppo_left_node, right_node)]['barcodes']
                    for barcode in barcode_weight_map.keys():
                        barcodes.add(barcode)
                if (left_node, oppo_right_node) in allele_linkage_graph.edges:
                    barcode_weight_map = allele_linkage_graph.edges[(left_node, oppo_right_node)]['barcodes']
                    for barcode in barcode_weight_map.keys():
                        barcodes.add(barcode)
                if (oppo_left_node, oppo_right_node) in allele_linkage_graph.edges:
                    barcode_weight_map = allele_linkage_graph.edges[(oppo_left_node, oppo_right_node)]['barcodes']
                    for barcode in barcode_weight_map.keys():
                        barcodes.add(barcode)
                barcode_pair_global_neg_map[(left_node, right_node)] = len(barcodes)
    # output
    df = n-1
    T = scipy.stats.t(df)
    all_pairs, selected_pairs = 0, 0
    if not os.path.exists('./output/{}/singular_cells'.format(opt.id)):
        os.mkdir('./output/{}/singular_cells'.format(opt.id))
    with open('./output/{}/singular_cells/singular_cell_linkage_{}.txt'.format(opt.id, opt.chr), 'w') as f:
        f.write('barcode\tvar\tgeno\tsupport\tp-value\toppo_support\ttotal_covered_barcodes\tglobal_oppo_support\tcorrect\tis_singular\n')
        for barcode, allele_pairs in barcode_pair_map.items():
            for allele_pair, link_value in allele_pairs.items():
                p_value = T.cdf((mean-link_value)/math.sqrt(var/10))
                if p_value > 0.05:
                    pass #continue
                gt_phasing_str = ''.join(list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', map(lambda x:x.split(':')[0], allele_pair))))
                is_correct = 0
                if check_haplotype(gt_phasing_str, ''.join(map(lambda x:x.split(':')[1], allele_pair))):
                    is_correct = 1
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(barcode, ','.join(map(lambda x:x.split(':')[0], allele_pair)), ''.join(map(lambda x:x.split(':')[1], allele_pair)),link_value, p_value, barcode_pair_neg_map[barcode][allele_pair], barcode_pair_global_neg_map[allele_pair], pair_neg_count_map[allele_pair],is_correct,pair_singular_map[allele_pair]))
    with open('./output/{}/cell_allele_connections_{}.txt'.format(opt.id, opt.chr), 'w') as f:
        f.write('barcode\tvar\tgeno\tsupport\n')
        for edge in allele_linkage_graph.edges:
            barcode_weight_map = allele_linkage_graph.edges[edge]['barcodes']
            for barcode, support in barcode_weight_map.items():
                f.write('{}\t{}\t{}\t{}\n'.format(barcode, ','.join(map(lambda x:x.split(':')[0], edge)), ''.join(map(lambda x:x.split(':')[1], edge)), support))
