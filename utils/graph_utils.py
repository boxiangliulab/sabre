import networkx as nx
import collections
import numpy as np
import math
import tqdm
from memory_profiler import profile

def output_graph_weights(opt, G:nx.Graph, vid_var_map):

    right_edge_weights = []
    wrong_edge_weights = []
    vid_pos_map = {}
    node_pos_map = {}
    for idx, node in enumerate(G.nodes):
        vid, allele = node.split(':')
        if vid not in vid_pos_map:
            vid_pos_map[vid] = idx
        pos_x = vid_pos_map[vid]

        var = vid_var_map[vid]
        pos_y = 0 if var.genotype_string.split('|')[0] == allele else 1
        node_pos_map[node] = (pos_x, pos_y)

        G.nodes[node]['x'] = float(pos_x)
        G.nodes[node]['y'] = float(pos_y)
    '''
    for a, b, data in G.edges.data():
        if node_pos_map[a][1] != node_pos_map[b][1]:
            G.edges[a,b]['right'] = 0
            wrong_edge_weights.append(data['weight'])
            continue
        G.edges[a,b]['right'] = 1
        right_edge_weights.append(data['weight'])
    '''
    nx.write_graphml(G, './output/original_graph_{}.graphml'.format(opt.restrict_chr))
    
    np.save('non-gcn-weight', nx.get_edge_attributes(G, 'weight'))
    np.save('right_edge_weights', right_edge_weights)
    np.save('wrong_edge_weights', wrong_edge_weights)

def create_graph(opt, allele_linkage_map, var_barcode_map, vid_var_map):
    '''
    Use allele linkage to create allele linkage graph.
    Notablly, in comparison with phaser, informations of linkage by read are completely reserved.
    '''
    G = nx.Graph()
    barcode_link_weights = []
    for (alle_1, alle_2), weight in allele_linkage_map.items():
        barcode_weight_map = var_barcode_map[(alle_1, alle_2)]
        barcode_link_weights += list(barcode_weight_map.values())
        G.add_edges_from([(alle_1, alle_2, {'prime_weight': weight, 'barcodes':barcode_weight_map})])
    G = graph_aggregation_and_update(G)
    return G, np.mean(barcode_link_weights), np.var(barcode_link_weights, ddof=1), len(barcode_link_weights)

def visualize_graph(G:nx.Graph, save_name):
    '''
    Save the visualization of graph in PDF format.
    '''
    if len(G.nodes) == 1:
        return
    
    nx.write_graphml(G, './output/graph_pdfs/{}.graphml'.format(save_name))

def find_connected_components(G:nx.Graph):
    '''
    Find connected components in allele linkage graph.
    '''
    subgraphs = []
    
    temp_G = nx.Graph()
    total_possible_pairs = 0

    for left_node, right_node in G.edges:
        temp_G.add_edge(left_node.split(':')[0], right_node.split(':')[0])
    for component in nx.connected_components(temp_G):
        total_possible_pairs += math.comb(len(component), 2)

    for component in nx.connected_components(G):
        subgraphs.append(G.subgraph(component))
        total_possible_pairs += math.comb(len(component), 2)
    return subgraphs, total_possible_pairs

def find_conflict_graphs(opt, subgraphs:list[nx.Graph], vid_var_map):
    '''
    Determine whether a graph contains conflict variants, w.r.t. 2 alleles related to 1 variant appears in a single graph.
    Implement this by get the node list of a subgraph, and find the corelated variant of it.
    If #alleles == #variants, then no conflict is found.

    Conflicted graphs needs further phasing, non-Conflicted graphs report phasing result.
    '''
    conflicted_graphs, nonconflicted_graphs = [], []
    for sg in subgraphs:
        alleles = list(sg.nodes)
        if len(alleles) == 1:
            continue
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
                pos_y = 0 if var.genotype_string.split('|')[0] == allele else 1
                node_pos_map[node] = (pos_x, pos_y)
                sg.nodes[node]['x'] = float(pos_x)
                sg.nodes[node]['y'] = float(pos_y)

            nx.write_graphml(sg, './output/succeed_graph_pdfs/{}.graphml'.format(id(sg)))
    
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

def split_graph_by_common_shortest_path(sg: nx.Graph, graph_name='', max_remove_node=1):
    '''
    Find all the conflicted alleles, and calculate every shortest-path between them.
    For every node in the graph, count its appearence time in all extracted shortest-paths.
    Sort these nodes descendingly, and remove #max_remove_node nodes accordingly.
    '''
    sg = nx.Graph(sg)
    common_nodes_on_short_paths = []
    conflicted_variants = find_conflict_alleles(sg)
    for cv in conflicted_variants:
        allele_0, allele_1 = '{}:0'.format(cv), '{}:1'.format(cv)
        common_nodes_on_short_paths += nx.shortest_path(sg, allele_0, allele_1, weight='weight')
        common_nodes_on_short_paths.remove(allele_0)
        common_nodes_on_short_paths.remove(allele_1)
    
    node_counter = dict(collections.Counter(common_nodes_on_short_paths))
    candidate_nodes = sorted(node_counter.keys(), key=lambda x: node_counter[x], reverse=True)

    for i, node in enumerate(candidate_nodes):
        sg.remove_node(node)
        is_conflicted = False
        for component in nx.connected_components(sg):
            if find_conflict_alleles(sg.subgraph(component)) != []:
                is_conflicted = True
                break
        
        if not is_conflicted or i == max_remove_node-1: return nx.connected_components(sg)
    
    return []

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

def split_graph_by_fiedler_vector(sg:nx.Graph, graph_name='', threshold=1e-2):
    '''
    Split the graph by fiedler vector to gain global optimized result.
    Fiedler cutting can reach lowest cut value and #node balance at the same time.
    '''
    if len(sg.nodes) <= 1:
        return []

    if find_conflict_alleles(sg) == []:
        return nx.connected_components(sg)

    final_partitions = []

    fiedler_vector = nx.fiedler_vector(sg)
    partitions = ([],[])
    for i, node in enumerate(sg.nodes):
        if abs(fiedler_vector[i]) < threshold:
            continue
        if fiedler_vector[i] < 0:
            partitions[0].append(node)
            continue
        partitions[1].append(node)
    
    partition_results = []
    for i in nx.connected_components(sg.subgraph(partitions[0])):
        partition_results.append(i)
    for i in nx.connected_components(sg.subgraph(partitions[1])):
        partition_results.append(i)

    for i in partition_results:
        if find_conflict_alleles(sg.subgraph(i)) == []:
            final_partitions.append(i)
        else:
            final_partitions+=split_graph_by_fiedler_vector(sg.subgraph(i))
    
    return final_partitions

def split_graph_by_extracting_singular_cells(opt, sg:nx.Graph):
    '''
    If removing one cell can resolve the conflict, then why not?
    '''
    if len(sg.nodes) <= 1:
        return []
    
    conflict_variants = find_conflict_alleles(sg)
    candidate_barcodes = []
    for cv in conflict_variants:
        allele_0, allele_1 = '{}:0'.format(cv), '{}:1'.format(cv)
        sp = nx.shortest_path(sg, allele_0, allele_1)
        if len(sp) == 3:
            critical_node = sp[1]
            critical_barcodes_0 = list(set(sg.nodes[allele_0]['cells'].keys()).intersection(set(sg.nodes[critical_node]['cells'].keys())))
            critical_barcodes_1 = list(set(sg.nodes[allele_1]['cells'].keys()).intersection(set(sg.nodes[critical_node]['cells'].keys())))
            candidate_barcodes.append((critical_barcodes_0, critical_barcodes_1))
    nodes = list(map(lambda x: list(sg.nodes)[x], conflict_variants))
    
def calculate_graph_difference(G:nx.Graph, H:nx.Graph):
    '''
    Remove all edges in H from G.
    '''
    G = nx.Graph(G)
    G.remove_edges_from(list(H.edges))
    return G

def resolve_conflict_graphs(opt, subgraphs: list[nx.Graph], phased_vars:set[str]):
    '''
    First, remove edges with 5% minimum edges weights.
    For sub-graphs with more than 100 nodes, remove no more than 0.2 * #nodes by split_graph_by_common_shortest_path.
    
    We then recursively cut the graph by fiedler cut till no conflict_alleles are presented.
    '''
    resolved_nodes = []
    removed_edges = []
    for sg in subgraphs:

        sg = nx.Graph(sg)
        if len(sg.nodes) > 2:
            edge_weights = list(nx.get_edge_attributes(sg, 'weight').values())
            cut_threshold = np.percentile(edge_weights, opt.edge_threshold)
            edges_to_remove = [(a,b) for a,b,attrs in sg.edges(data=True) if attrs['weight']<=cut_threshold]
            sg.remove_edges_from(edges_to_remove)
        
        residual_graph = nx.Graph() 
        for components in nx.connected_components(sg):
            cleared_sg = sg.subgraph(components)
            if len(cleared_sg.nodes) <= 1:
                continue

            graph_name = list(cleared_sg.nodes)[0]
            final_partitions = []
            if opt.shortest_path:
                if opt.remove_node == 'auto':
                    if len(cleared_sg.nodes) < 100:
                        num_remove_node = 1
                    else:
                        num_remove_node = int(math.floor(len(cleared_sg.nodes)*0.2))
                else:
                    num_remove_node = int(opt.remove_node)
                _1_partitions = split_graph_by_common_shortest_path(cleared_sg, graph_name, max_remove_node=num_remove_node)
                for partition in _1_partitions:
                
                    if find_conflict_alleles(sg.subgraph(partition)) == []:
                        final_partitions.append(partition)
                        continue
                    final_partitions += split_graph_by_fiedler_vector(sg.subgraph(partition), graph_name, threshold=opt.fiedler_threshold)
            else:
                final_partitions += split_graph_by_fiedler_vector(cleared_sg, graph_name, threshold=opt.fiedler_threshold)
            
            residual_graph = nx.Graph(sg)

            for partition in final_partitions:
                resolved_subgraph = cleared_sg.subgraph(partition)
                if opt.verbose:
                    visualize_graph(sg, '{}_1'.format(graph_name))
                    visualize_graph(resolved_subgraph, '{}_0'.format(graph_name))
                if len(partition) > 1: resolved_nodes.append(list(partition))
                residual_graph = calculate_graph_difference(residual_graph, sg.subgraph(partition))

        for component in nx.connected_components(residual_graph):
            if len(component) > 1:
                removed_edges.append(residual_graph.subgraph(component))        
        
    return resolved_nodes, removed_edges

def extract_nonconflicted_nodes(subgraphs: list[nx.Graph]):
    '''
    Extract node sets from subgraphs.
    Those node sets can be considered as haplotypes.
    '''
    nonconflicted_nodes = []
    phased_vars = []
    for sg in subgraphs:
        nonconflicted_nodes.append(list(sg.nodes))
        phased_vars += list(map(lambda x: x.split(':')[0], list(sg.nodes)))
    
    return nonconflicted_nodes, set(phased_vars)

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
            node_popularity += G.edges[node, neighbor_node]['prime_weight']
        G.nodes[node]['popularity'] = node_popularity
    
    for n_a, n_b in G.edges:
        G.edges[n_a, n_b]['weight'] = round(G.edges[n_a, n_b]['prime_weight'] / max(G.nodes[n_a]['popularity'] / G.nodes[n_b]['popularity'], G.nodes[n_b]['popularity']/ G.nodes[n_a]['popularity']), 4)
    return G
