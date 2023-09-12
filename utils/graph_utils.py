import networkx as nx
import collections
import numpy as np

def output_graph_weights(G:nx.Graph, vid_var_map):

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
    for a, b, data in G.edges.data():
        if node_pos_map[a][1] != node_pos_map[b][1]:
            wrong_edge_weights.append(data['weight'])
            continue
        right_edge_weights.append(data['weight'])
    np.save('gcn-weight', nx.get_edge_attributes(G, 'weight'))
    np.save('right_edge_weights-gcn', right_edge_weights)
    np.save('wrong_edge_weights-gcn', wrong_edge_weights)


def create_graph(opt, allele_linkage_map, vid_var_map):
    '''
    Use allele linkage to create allele linkage graph.
    Notablly, in comparison with phaser, infomations of linkage by read are completely reserved.
    '''
    G = nx.Graph()
    for (alle_1, alle_2), weight in allele_linkage_map.items():
        G.add_edges_from([(alle_1, alle_2, {'weight': weight})])
    
    # nx.write_gexf(G, './output/original_graph_{}.gexf'.format(opt.restrict_chr))
    G = graph_strenthen(G)
    # G = graph_aggregation_and_update(G)
    nx.write_gexf(G, './output/original_graph_{}.gexf'.format(opt.restrict_chr))

    if opt.verbose:
        output_graph_weights(G, vid_var_map)

    cut_threshold = 1
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
    
    nx.write_graphml(G, './output/graph_pdfs/{}.graphml'.format(save_name))

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
    Find all the conflicted alleles, find common node in short paths between conflicted nodes and remove them.
    '''

    sg = nx.Graph(sg)
    common_nodes_on_short_paths = []
    conflicted_variants = find_conflict_alleles(sg)
    for cv in conflicted_variants:
        allele_0, allele_1 = '{}:0'.format(cv), '{}:1'.format(cv)
        common_nodes_on_short_paths += nx.shortest_path(sg, allele_0, allele_1, weight='weight')
    
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
    '''

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

def resolve_conflict_graphs(opt, subgraphs: list[nx.Graph], phased_vars:set[str]):
    '''
    By now, this function can only resolve conflict graphs with only one conflict allele.
    The underlying algorithm is a Max-flow/Min-Cut graph algorithm, namely preflow-push algorithm.
    The time complexity is greatly reduced compared with Phaser, with O(n^2 \dot \sqrt(m)) for Faser and O(2^n) in Phaser. 

    We first find conflict alleles in the graph, then consider them to be the source node and sink node for the flow.
    Then utilize networkx.minimum_cut() to split the subgraph by minimum cut.
    '''
    resolved_nodes = []
    for sg in subgraphs:

        sg = nx.Graph(sg)
        remove_nodes = []
        for node in sg.nodes:
            if node.split(':')[0] in phased_vars:
                remove_nodes.append(node)
        sg.remove_nodes_from(remove_nodes)

        for components in nx.connected_components(sg):
            cleared_sg = sg.subgraph(components)
            if len(cleared_sg.nodes) <= 1:
                continue

            graph_name = list(cleared_sg.nodes)[0]
            # entropy = calculate_graph_entropy(sg)
            final_partitions = split_graph_by_common_shortest_path(cleared_sg, graph_name, max_remove_node=2)
            # final_partitions = split_graph_by_min_cut(sg, graph_name)
            final_partitions = split_graph_by_fiedler_vector(cleared_sg, graph_name, threshold=1e-2)

            for partition in final_partitions:
                resolved_subgraph = cleared_sg.subgraph(partition)
                if opt.verbose:
                    visualize_graph(sg, '{}_1'.format(graph_name))
                    visualize_graph(resolved_subgraph, '{}_0'.format(graph_name))
                if len(partition) > 1: resolved_nodes.append(list(partition))
        
    return resolved_nodes

def extract_nonconflicted_nodes(subgraphs: list[nx.Graph]):
    '''
    Extract node sets from subgraphs.
    Those node sets can be considered as hyplotypes.
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
            node_popularity += G.edges[node, neighbor_node]['weight']
        G.nodes[node]['popularity'] = node_popularity
    
    for n_a, n_b in G.edges:
        G.edges[n_a, n_b]['weight'] = (G.edges[n_a, n_b]['weight'] / (G.nodes[n_a]['popularity'] / G.degree[n_a])) * (G.edges[n_a, n_b]['weight'] / (G.nodes[n_b]['popularity'] / G.degree[n_b]))
    return G

def graph_strenthen(G:nx.Graph):
    '''
    If:
    A:0 ----1---- B:1
    A:1 ----3---- B:0
    Then:
    A:0 ----4---- B:1
    A:1 ----4---- B:0
    '''

    var_conn_set = set()
    for a, b, data in G.edges.data():
        var_a, allele_a = a.split(':')
        var_b, allele_b = b.split(':')

        if (a,b) not in var_conn_set:
            counter_a, counter_b = '{}:{}'.format(var_a, 1-int(allele_a)), '{}:{}'.format(var_b, 1-int(allele_b))
            var_conn_set.add((a,b))
            var_conn_set.add((counter_b,counter_a))
            if (counter_a, counter_b) in G.edges:
                temp_weight = data['weight']
                G.edges[(a,b)]['weight'] = G.edges[(a,b)]['weight'] + G.edges[(counter_a, counter_b)]['weight']
                G.edges[(counter_a, counter_b)]['weight'] = G.edges[(counter_a, counter_b)]['weight'] + temp_weight

    return G