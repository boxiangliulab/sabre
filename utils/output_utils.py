import networkx as nx
import matplotlib.pyplot as plt
import itertools

def get_opposite_allele(allele, is_opposite=True):
    if is_opposite:
        return allele.split(':')[0]+':{}'.format(1-int(allele.split(':')[1]))
    return allele

def check_hyplotype(ground_truth, predict):

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
    As many of these hyplotypes are actually not disjointed from each other, like
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
    
    total_hyp = 0
    correct_hyp = 0
    total_predict = 0
    correct_predict = 0
    already_reported_set = []
    total_nodes = 0

    final_hyplotypes = de_duplicate(nonconflicted_nodes + resolved_conflicted_nodes)

    # already_seen_var = set()
    # for nodes in final_hyplotypes:
    #     for node in nodes:
    #         if node.split(':')[0] in already_seen_var:
    #             print(node)
    #         else:
    #             already_seen_var.add(node.split(':')[0])


    with open('./output/chr_{}_haplotypes.tsv'.format(opt.restrict_chr), 'w') as f:
        f.write('haplotype\tpredicted_phasing\tgt_phasing\tcorrection\n')
        for nodes in final_hyplotypes:
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
            if check_hyplotype(gt_phasing_str, predicted_phasing_str):
                is_correct = 1

            total_hyp += 1
            correct_hyp += is_correct

            f.write('{}\t{}\t{}\t{}\n'.format(vids_string, predicted_phasing_str, gt_phasing_str, is_correct))


    with open('./output/chr_{}_conflicted_hyplotypes.tsv'.format(opt.restrict_chr), 'w') as f:
        for nodes in resolved_conflicted_nodes:
            var_phasing_list = list(map(lambda x: x.split(':'), nodes))
            vids, phasing = zip(*var_phasing_list)
            var_set = set(vids)
            already_reported_set.append(var_set)
            vids_string = ','.join(vids)
            predicted_phasing_str = ''.join(phasing)

            #Get ground_truth phasing result
            gt_phasing_str = ''.join(list(map(lambda x: vid_var_map[x].genotype_string.split('|')[0] if vid_var_map[x].is_phased else '-', vids)))
            is_correct = 0
            if check_hyplotype(gt_phasing_str, predicted_phasing_str):
                is_correct = 1

            total_predict += 1
            correct_predict += is_correct
            
            f.write('{}\t{}\t{}\t{}\t1\n'.format(vids_string, predicted_phasing_str, gt_phasing_str, is_correct))

    return total_hyp, correct_hyp, total_predict, correct_predict, total_nodes