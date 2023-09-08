import networkx as nx
import matplotlib.pyplot as plt
import collections
from rich import print
import math
import numpy as np
from scipy.stats import binom
import pickle

def create_graph(opt, allele_linkage_map):
    '''
    Use allele linkage to create allele linkage graph.
    Notablly, in comparison with phaser, infomations of linkage by read are completely reserved.
    '''
    G = nx.Graph()
    for (alle_1, alle_2), weight in allele_linkage_map.items():
        G.add_edges_from([(alle_1, alle_2, {'weight': weight})])
    
    nx.write_gexf(G, './output/original_graph_{}.gexf'.format(opt.restrict_chr))
    G = graph_aggregation_and_update(G)
    np.save('gcn-weight',nx.get_edge_attributes(G, 'weight'), allow_pickle=True)
    edge_weights = list(nx.get_edge_attributes(G, 'weight').values())
    cut_threshold = np.percentile(edge_weights, opt.edge_threshold)
    edges_to_remove = [(a,b) for a,b,attrs in G.edges(data=True) if attrs['weight']<=cut_threshold]
    # Unfreeze the graph, otherwise will raise Frozen graph can't be modified Error.
    G.remove_edges_from(edges_to_remove)    

    return G

def visualize_graph(G:nx.Graph, save_name):
    '''
    Save the visualization of graph in PDF format.
    '''
    if len(G.nodes) == 1:
        return
    
    nx.write_gexf(G, './output/graph_pdfs/{}.gexf'.format(save_name))

def find_connected_components(G:nx.Graph):
    '''
    Find connected components in allele linkage graph.
    '''
    subgraphs = []
    for component in nx.connected_components(G):
        subgraphs.append(G.subgraph(component))
    return subgraphs

def find_conflict_graphs(opt, subgraphs:list[nx.Graph], vid_var_map):
    '''
    Determine whether a graph contains conflict variants, w.r.t. 2 alleles related to 1 variant appears in a single graph.
    Implement this by get the node list of a subgraph, and find the corelated variant of it.
    If #alleles == #variants, then no conflict is found.

    Conflicted graphs needs further phasing, non-Conflicted graphs report phasing result.
    
    v0.4 Update: Pre-split the graph, by removing edges with weight lower than the 5th percentile.
    '''
    conflicted_graphs, nonconflicted_graphs = [], []
    for sg in subgraphs:
        alleles = list(sg.nodes)
        variants = set(map(lambda x: x.split(':')[0], alleles))
        if len(alleles) == len(variants):
            nonconflicted_graphs.append(sg)
            continue

        conflicted_graphs.append(sg)
        if opt.verbose:
            vid_pos_map = {}
            node_pos_map = {}
            for idx, node in enumerate(sg.nodes):
                vid, allele = node.split(':')
                if vid not in vid_pos_map:
                    vid_pos_map[vid] = idx
                pos_x = vid_pos_map[vid]

                var = vid_var_map[vid]
                pos_y = 0 if var.genotype_string.split('|')[0] == allele else 40
                node_pos_map[node] = (pos_x, pos_y)
                sg.nodes[node]['x'] = float(pos_x)
                sg.nodes[node]['y'] = float(pos_y)
            nx.write_graphml(sg, './output/succeed_graph_pdfs/{}.graphml'.format(id(sg)))

            # nx.draw_networkx_nodes(sg, node_pos_map, node_size=5)
            # # edges
            # nx.draw_networkx_edges(sg, node_pos_map, edgelist=sg.edges, width=2, connectionstyle="arc3,rad=0.1", arrows=True)
            # edge_labels = nx.get_edge_attributes(sg, 'weight')
            # nx.draw_networkx_edges(sg, node_pos_map, edge_labels)
            # ax = plt.gca()
            # ax.margins(0.08)
            # plt.axis("off")
            # plt.tight_layout()
            # plt.show()
    
    return conflicted_graphs, nonconflicted_graphs

def find_conflict_alleles(G: nx.Graph):
    '''
    Find the conflict allele in G, by count the appearance of each variant.
    Vairants with appearance more than once are considered conflicted.
    '''
    alleles_nodes = list(G.nodes)
    variant_list = list(map(lambda x: x.split(':')[0], alleles_nodes))
    variant_count = collections.Counter(variant_list)
    return list(filter(lambda x: variant_count[x] == 2, variant_count.keys()))

def split_graph_by_min_cut(sg: nx.Graph, graph_name=''):
    '''
    Recursively cut the graph by min-cut/max-flow algorithm till there is no conflict allele inside.
    '''

    final_partitions = []

    conflict_variants = find_conflict_alleles(sg)
    # resolve alleles with higer degree in higher priority.
    variants_degree_map = {}
    for var in conflict_variants:
        degree_0 = sg.degree(var+':0')
        degree_1 = sg.degree(var+':1')
        variants_degree_map[var] = degree_0 + degree_1
    
    resolve_node = sorted(variants_degree_map.keys(), key=lambda x:variants_degree_map[x], reverse=False)[0]
    
    cut_value, partitions = nx.minimum_cut(sg, resolve_node+':0', resolve_node+':1', capacity='weight')

    if find_conflict_alleles(sg.subgraph(partitions[0])) == []:
        final_partitions.append(partitions[0])
    else:
        final_partitions+=split_graph_by_min_cut(sg.subgraph(partitions[0]))
    
    if find_conflict_alleles(sg.subgraph(partitions[1])) == []:
        final_partitions.append(partitions[1])
    else:
        final_partitions+=split_graph_by_min_cut(sg.subgraph(partitions[1]))
    
    return final_partitions


def calculate_graph_entropy(G:nx.Graph):
    '''
    Calculate the confusion level of the given graph, by using p=softmax(weight^{-1}).
    '''
    adjusted_edge_weights = []
    for edge in G.edges:
        adjusted_edge_weights.append(G.edges[edge]['weight'])
    
    adjusted_edge_weights = 1/np.array(adjusted_edge_weights)
    adjusted_edge_weights = adjusted_edge_weights - max(adjusted_edge_weights)
    p = np.exp(adjusted_edge_weights) / np.sum(np.exp(adjusted_edge_weights))
    return np.sum(p*np.log(p+1e-99))

def resolve_conflict_graphs(opt, subgraphs: list[nx.Graph]):
    '''
    By now, this function can only resolve conflict graphs with only one conflict allele.
    The underlying algorithm is a Max-flow/Min-Cut graph algorithm, namely preflow-push algorithm.
    The time complexity is greatly reduced compared with Phaser, with O(n^2 \dot \sqrt(m)) for Faser and O(2^n) in Phaser. 

    We first find conflict alleles in the graph, then consider them to be the source node and sink node for the flow.
    Then utilize networkx.minimum_cut() to split the subgraph by minimum cut.

    If the conflict graph is:
    A_0:----i----:B:----j----:A_1
    We will do a binomial test to determine whether this graph is significant.
    '''
    resolved_nodes = []
    for sg in subgraphs:
        
        graph_name = list(sg.nodes)[0]

        #Do GCN
        # sg = graph_aggregation_and_update(sg)
        # sgs = split()
        entropy = calculate_graph_entropy(sg)
        final_partitions = split_graph_by_min_cut(sg, graph_name)

        for partition in final_partitions:
            resolved_subgraph = sg.subgraph(partition)
            if opt.verbose:
                visualize_graph(resolved_subgraph, '{}_0'.format(graph_name))
            if len(partition) > 1: resolved_nodes.append([list(partition), entropy])
    
    return resolved_nodes

def extract_nonconflicted_nodes(subgraphs: list[nx.Graph]):
    '''
    Extract node sets from subgraphs.
    Those node sets can be considered as hyplotypes.
    '''
    nonconflicted_nodes = []
    for sg in subgraphs:
        nonconflicted_nodes.append(list(sg.nodes))
    
    return nonconflicted_nodes

def graph_aggregation_and_update(G:nx.Graph):
    '''
    This algorithm is much like a GCN.
    Step 1: Calculate \sum(weight_e) for each node, which is the popularity of each node P_n.
    Step 2: Normalize each edge weight by the sum of the popularity of its paired end.
    '''
    
    # We first calculate weight sum of each node
    # nx.set_node_attributes(G, 0, 'popularity')

    for node in G.nodes:
        node_popularity = 0
        for neighbor_node in G.neighbors(node):
            node_popularity += G.edges[node, neighbor_node]['weight']
        G.nodes[node]['popularity'] = node_popularity
    
    for n_a, n_b in G.edges:
        G.edges[n_a, n_b]['weight'] = round(G.edges[n_a, n_b]['weight'] / max(G.nodes[n_a]['popularity'] / G.nodes[n_b]['popularity'], G.nodes[n_b]['popularity']/ G.nodes[n_a]['popularity']), 4)
        # G.edges[n_a, n_b]['weight'] = round(G.edges[n_a, n_b]['weight'] / max((G.nodes[n_a]['popularity']*G.degree[n_b]) / (G.nodes[n_b]['popularity']*G.degree[n_a]) , (G.nodes[n_b]['popularity']*G.degree[n_a]) / (G.nodes[n_a]['popularity']*G.degree[n_b]) ), 4)
    return G