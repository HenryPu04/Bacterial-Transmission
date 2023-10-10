# Firstname Lastname
# NetID
# COMP 182 Spring 2021 - Homework 5, Problem 4
# You can import any standard library, as well as Numpy and Matplotlib.
# You can use helper functions from provided.py, and autograder.py,
# but they have to be copied over here.
# Your code here...
from collections import deque
from copy import *
import numpy
import matplotlib
from typing import Tuple
from collections import *
from graphviz import Digraph
def bfs(graph, startnode):
    """
        Perform a breadth-first search on digraph graph starting at node startnode.
        
        Arguments:
        graph -- directed graph
        startnode - node in graph to start the search from
        
        Returns:
        The distances from startnode to each node
    """
    dist = {}
    
    # Initialize distances
    for node in graph:
        dist[node] = float('inf')
    dist[startnode] = 0
    
    # Initialize search queue
    queue = deque([startnode])
    
    # Loop until all connected nodes have been explored
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def compute_rdmst(graph, root):
    """
        This function checks if:
        (1) root is a node in digraph graph, and
        (2) every node, other than root, is reachable from root
        If both conditions are satisfied, it calls compute_rdmst_helper
        on (graph, root).
        
        Since compute_rdmst_helper modifies the edge weights as it computes,
        this function reassigns the original weights to the RDMST.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node id.
        
        Returns:
        An RDMST of graph rooted at r and its weight, if one exists;
        otherwise, nothing.
    """
    
    if root not in graph:
        print ("The root node does not exist")
        return
    
    distances = bfs(graph, root)
    for node in graph:
        if distances[node] == float('inf'):
            print ("The root does not reach every other node in the graph")
            return

    rdmst = compute_rdmst_helper(graph, root)
    
    # reassign the original edge weights to the RDMST and computes the total
    # weight of the RDMST
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = graph[node][nbr]
            rdmst_weight += rdmst[node][nbr]

    return (rdmst,rdmst_weight)

def compute_rdmst_helper(graph,root):
    """
        Computes the RDMST of a weighted digraph rooted at node root.
        It is assumed that:
        (1) root is a node in graph, and
        (2) every other node in graph is reachable from root.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node in graph.
        
        Returns:
        An RDMST of graph rooted at root. The weights of the RDMST
        do not have to be the original weights.
        """
    
    # reverse the representation of graph
    rgraph = reverse_digraph_representation(graph)
    
    # Step 1 of the algorithm
    modify_edge_weights(rgraph, root)
    
    # Step 2 of the algorithm
    rdst_candidate = compute_rdst_candidate(rgraph, root)
    
    # compute a cycle in rdst_candidate
    cycle = compute_cycle(rdst_candidate)
    
    # Step 3 of the algorithm
    if not cycle:
        return reverse_digraph_representation(rdst_candidate)
    else:
        # Step 4 of the algorithm
        
        g_copy = deepcopy(rgraph)
        g_copy = reverse_digraph_representation(g_copy)
        
        # Step 4(a) of the algorithm
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        #cstar = max(contracted_g.keys())
        
        # Step 4(b) of the algorithm
        new_rdst_candidate = compute_rdmst_helper(contracted_g, root)
        
        # Step 4(c) of the algorithm
        rdmst = expand_graph(reverse_digraph_representation(rgraph), new_rdst_candidate, cycle, cstar)
        
        return rdmst


### Functions for use in Problem 4

def infer_transmap(gen_data, epi_data, patient_id):
    """
        Infers a transmission map based on genetic
        and epidemiological data rooted at patient_id
        
        Arguments:
        gen_data -- filename with genetic data for each patient
        epi_data -- filename with epidemiological data for each patient
        patient_id -- the id of the 'patient 0'
        
        Returns:
        The most likely transmission map for the given scenario as the RDMST 
        of a weighted, directed, complete digraph
        """
    
    complete_digraph = construct_complete_weighted_digraph(gen_data, epi_data)
    return compute_rdmst(complete_digraph, patient_id)


def read_patient_sequences(filename):
    """
        Turns the bacterial DNA sequences (obtained from patients) into a list containing tuples of
        (patient ID, sequence).
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A list of (patient ID, sequence) tuples.
        """
    sequences = []
    with open(filename) as f:
        line_num = 0
        for line in f:
            if len(line) > 5:
                patient_num, sequence = line.split("\t")
                sequences.append( (int(patient_num), ''.join(e for e in sequence if e.isalnum())) )
    return sequences

def read_patient_traces(filename):
    """
        Reads the epidemiological data file and computes the pairwise epidemiological distances between patients
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A dictionary of dictionaries where dict[i][j] is the
        epidemiological distance between i and j.
    """
    trace_data = []
    patient_ids = []
    first_line = True
    with open(filename) as f:
        for line in f:
            if first_line:
                patient_ids = line.split()
                patient_ids = list(map(int, patient_ids))
                first_line = False
            elif len(line) > 5:
                trace_data.append(line.rstrip('\n'))
    return compute_pairwise_epi_distances(trace_data, patient_ids)

def compute_pairwise_gen_distances(sequences, distance_function):
    """
        Computes the pairwise genetic distances between patients (patients' isolate genomes)
        
        Arguments:
        sequences -- a list of sequences that correspond with patient id's
        distance_function -- the distance function to apply to compute the weight of the 
        edges in the returned graph
        
        Returns:
        A dictionary of dictionaries where gdist[i][j] is the
        genetic distance between i and j.
        """
    gdist = {}
    cultures = {}
    
    # Count the number of differences of each sequence
    for i in range(len(sequences)):
        patient_id = sequences[i][0]
        seq = sequences[i][1]
        if patient_id in cultures:
            cultures[patient_id].append(seq)
        else:
            cultures[patient_id] = [seq]
            gdist[patient_id] = {}
    # Add the minimum sequence score to the graph
    for pat1 in range(1, max(cultures.keys()) + 1):
        for pat2 in range(pat1 + 1, max(cultures.keys()) + 1):
            min_score = float("inf")
            for seq1 in cultures[pat1]:
                for seq2 in cultures[pat2]:
                    score = distance_function(seq1, seq2)
                    if score < min_score:
                        min_score = score
            gdist[pat1][pat2] = min_score
            gdist[pat2][pat1] = min_score
    return gdist



### HELPER FUNCTIONS. ###

def find_first_positives(trace_data):
    """
        Finds the first positive test date of each patient
        in the trace data.
        Arguments:
        trace_data -- a list of data pertaining to location
        and first positive test date
        Returns:
        A dictionary with patient id's as keys and first positive
        test date as values. The date numbering starts from 0 and
        the patient numbering starts from 1.
        """
    first_pos = {}
    for pat in range(len(trace_data[0])):
        first_pos[pat + 1] = None
        for date in range(len(trace_data)):
            if trace_data[date][pat].endswith(".5"):
                first_pos[pat + 1] = date
                break
    return first_pos



def compute_epi_distance(pid1, pid2, trace_data, first_pos1, first_pos2, patient_ids):
    """
        Computes the epidemiological distance between two patients.
        
        Arguments:
        pid1 -- the assumed donor's index in trace data
        pid2 -- the assumed recipient's index in trace data
        trace_data -- data for days of overlap and first positive cultures
        first_pos1 -- the first positive test day for pid1
        first_pos2 -- the first positive test day for pid2
        patient_ids -- an ordered list of the patient IDs given in the text file
        
        Returns:
        Finds the epidemiological distance from patient 1 to
        patient 2.
        """
    first_overlap = -1
    assumed_trans_date = -1
    pid1 = patient_ids.index(pid1)
    pid2 = patient_ids.index(pid2)
    # Find the first overlap of the two patients
    for day in range(len(trace_data)):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            first_overlap = day
            break
    if (first_pos2 < first_overlap) | (first_overlap < 0):
        return len(trace_data) * 2 + 1
    # Find the assumed transmission date from patient 1 to patient 2
    for day in range(first_pos2, -1, -1):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            assumed_trans_date = day
            break
    sc_recip = first_pos2 - assumed_trans_date

    if first_pos1 < assumed_trans_date:
        sc_donor = 0
    else:
        sc_donor = first_pos1 - assumed_trans_date
    return sc_donor + sc_recip



def compute_pairwise_epi_distances(trace_data, patient_ids):
    """
        Turns the patient trace data into a dictionary of pairwise 
        epidemiological distances.
        
        Arguments:
        trace_data -- a list of strings with patient trace data
        patient_ids -- ordered list of patient IDs to expect
        
        Returns:
        A dictionary of dictionaries where edist[i][j] is the
        epidemiological distance between i and j.
        """
    edist = {}
    proc_data = []
    # Reformat the trace data
    for i in range(len(trace_data)):
        temp = trace_data[i].split()[::-1]
        proc_data.append(temp)
    # Find first positive test days and remove the indication from the data
    first_pos = find_first_positives(proc_data)
    for pid in first_pos:
        day = first_pos[pid]
        proc_data[day][pid - 1] = proc_data[day][pid - 1].replace(".5", "")
    # Find the epidemiological distance between the two patients and add it
    # to the graph
    for pid1 in patient_ids:
        edist[pid1] = {}
        for pid2 in patient_ids:
            if pid1 != pid2:
                epi_dist = compute_epi_distance(pid1, pid2, proc_data,
                                                first_pos[pid1], first_pos[pid2], patient_ids)
                edist[pid1][pid2] = epi_dist
    return edist
def reverse_digraph_representation(graph: dict) -> dict:
    """
    Function reverse_digraph_representation(graph) that takes as input a weighted
digraph graph in the standard representation.
    Returns exactly the same weighted digraph graph but in the reversed
representation.
    """
    reversed_graph = {node: {} for node in graph}
    #iterate through each node
    for node in graph:
        for neighbor, weight in graph[node].items():
            #reverse edge direction
            reversed_graph[neighbor][node] = weight
    return reversed_graph
def modify_edge_weights(rgraph: dict, root: int) -> None:
    """
    Function that takes as input a weighted digraph graph in the reversed
representation and a node root, and modifies the edge weights of
    graph according to Lemma 2.
    The function does not return any values.
    """
    def find_min_value(mapping):
        """
        Finds the minimum value in mapping
        """
        min_value = None
        for value in mapping.values():
            if min_value is None or value < min_value:
                min_value = value
        return min_value
    for node, neighbors in rgraph.items():
        if node==root:
            rgraph[node]={}
        else:
            smallest_weight=find_min_value(neighbors)
            for neighbor, weight in neighbors.items():
                #iterate through each edge and subtract the min edge weight
                neighbors[neighbor]=weight-smallest_weight
    #print(reverse_digraph_representation(rgraph))
def compute_rdst_candidate(rgraph: dict, root: int) -> dict:
    """
    Inputs:
        rgraph (dict): a weighted digraph in reversed representation.
        root (int): integer representing a root node in that inputted graph.
    Returns:
        dict: a weighted digraph in reversed rerpesentation that is a potential
rdst candidate
    """
    candidate = {node : {} for node in rgraph}
    candidate[root] = {}
    for node in rgraph:
        smallest_edge = 0
        smallest_neighbor = 0
        neighbors = rgraph[node]
        if root != node:
            #get correspdonding node that the smallest edge comes from
            smallest_neighbor = min(neighbors, key = neighbors.get)
            #get smallest edge the node receives
            smallest_edge = neighbors[smallest_neighbor]
            #build candidate graph
            candidate[node][smallest_neighbor] = smallest_edge
    return candidate
def compute_cycle(rdst_candidate: dict) -> tuple:
    """
    Inputs:
        rdst_candidate (dict): weighted digraph in reversed representation
    Returns:
        tuple: a tuple of integers that represents nodes in a cycle within the
inputted graph
    """
    visited = set()
    stack = set()
    cycle = []
    def dfs_cycle(node, parent):

            visited.add(node)
            stack.add(node)
            for neighbor in rdst_candidate.get(node, []):
                #iterates over neighbors of node and recursively calls dfs
                if neighbor not in visited:
                    parent[neighbor] = node
                    if dfs_cycle(neighbor, parent):
                        return True
                #if neighbor has already been visited and is still in stack then a cycle is found,
                elif neighbor in stack:
                    cycle.append(neighbor)
                    #adds all the nodes that make up the cycle
                    while node != neighbor:
                        cycle.append(node)
                        node = parent[node]
                    cycle.append(neighbor)
                    return True
            #no cycle found
            stack.remove(node)
            return False
    parent = {}
    #iterates over each node and runs the recursive dfs until dfs returns true
    for node in rdst_candidate:
        if node not in visited:
            parent[node] = None
            if dfs_cycle(node, parent):
                return tuple(set(cycle))
    #no cycle
    return None
def contract_cycle(graph: dict, cycle: tuple) -> Tuple[dict, int]:
    """
    Args:
        graph (dict): a weighed digraph in the standard representation
        cycle (tuple): a tuple of integers that represents nodes that make up a
cycle within that graph
    Returns:
        Tuple[dict, int]: a tuple that contains a weighted digraph that results
from contracting the cycle
        and a integer that represents the new node used in place of the cycle in
the graph """
    #contracts the cycle, by taking the cycle computed from compute_cycle
    #first, create a new graph
    graph_star={}
    #first, add all the nodes not in cycle
    for node in graph:
        if node not in cycle:
            graph_star[node]={}
    #find cstar
    c_star=-1
    for node in graph:
        if node>c_star:
            c_star=node
        c_star+=1
    graph_star[c_star]={}
    #print("the cstar is",c_star)
    #for every edge in graph, do the following cases
    for node, nbr_list in graph.items():
        for nbr in nbr_list.keys():
            #check if u is not in C and v is not in C
            if node not in cycle and nbr not in cycle :

                graph_star[node][nbr]=graph[node][nbr]
            elif node not in cycle and nbr in cycle:
                if c_star not in graph_star[node].keys():
                    graph_star[node][c_star]=graph[node][nbr]
            #if (node,c_star) does not already exist
                else:
                    if graph[node][nbr]<graph_star[node][c_star]:
                        graph_star[node][c_star]=graph[node][nbr]
            #when there are parallel edges, weight resets when new weight
            elif node in cycle and nbr not in cycle:
                if nbr not in graph_star[c_star].keys():
                    graph_star[c_star][nbr]=graph[node][nbr]
                    #when (c_star,nbr) doesnt exist
                else:
                    if graph[node][nbr]<graph_star[c_star][nbr]:
                        graph_star[c_star][nbr]=graph[node][nbr]
                    #when there are parallel edges, reset weight if new weight is smaller
    return tuple([graph_star,c_star])
def expand_graph(graph: dict, rdst_candidate: dict, cycle: tuple, cstar: int) ->dict:
    """ Input:
            graph (dict): graph in standard representation whose cycle is contracted
            rdst_candidate (dict): weighted digraph in standard representation,
    computed on contracted version
            of graph
            cycle (tuple): tuple of nodes on cycle that was contracted
            cstar (int): integer that labels the node that replaces the contracted
    cycle.
        Returns:
            dict: a weighted digraph that results from expanding the cycle in
    rdst_candidate
        """
    reverse_candidate = reverse_digraph_representation(rdst_candidate)
    rgraph = reverse_digraph_representation(graph)
    result = deepcopy(rdst_candidate)
    # Instantiate values need to find
    for node in cycle:
        result[node] = {}

# Find vstar
    u_node = list(reverse_candidate[cstar].keys())[0]
    weight = float("inf")
    v_node = -1
    for node in cycle:
        if node in graph[u_node]:
            if graph[u_node][node] < weight:
                weight = graph[u_node][node]
                v_node = node
    result[u_node][v_node] = weight
    # Get ustar,v edges
    for node in result[cstar]:
        weight = float("inf")
        u_star = 0
        for nbr in rgraph[node]:
            if nbr in cycle:
                if rgraph[node][nbr] < weight:
                    weight = rgraph[node][nbr]
                    u_star = nbr
        result[u_star][node] = weight
    # Remove cstar from nodes
    result.pop(cstar)
    # Remove cstar from edges
    result[u_node].pop(cstar)
    index = cycle.index(v_node)
    clist = list(cycle)
    clist = clist[index:] + clist[:index]
    clist = clist[1:] + clist[:1]
    # Create graph with last nodes
    for idx in range(len(clist) - 1):
        result[clist[idx + 1]][clist[idx]] = rgraph[clist[idx]][clist[idx + 1]]
    return result






def compute_genetic_distance(sequence1, sequence2):
    """
    Inputs:
        sequence1: First sequence representing a bacterial genome.
        sequence2:: Second sequence representing a baceterial genome
    Returns:
        distance (int): a integer representing the hamming distance between the two
genomes """
    distance = 0
    for i, val in enumerate(sequence1):
        if val != sequence2[i]:
            distance += 1
    return distance
def construct_complete_weighted_digraph(genetic_data, epi_data):
    """
    Inputs:
        genetic_data: file with patient sequences
        epi_data: file with the patient traces for the epidemiological data
    Returns:
        graph (dict): a complete, weighed, digraph whose ndoes are the patients
        and whose edge weights are based on Equation (1)
    """
    genetic_data = read_patient_sequences(genetic_data)
    graph = {}
    #get patient ids from genetic_data
    for tup in genetic_data:
        graph[tup[0]] = {}
    epi_distances= read_patient_traces(epi_data)
    genetic_distances = compute_pairwise_gen_distances(genetic_data,compute_genetic_distance)
    #print(genetic_distances)
    #determine max of matrix E
    max=-1
    for i in epi_distances.values():
        for dist in i.values():
            if dist > max:
                max = dist
    #find edge weights

    for patient in genetic_distances.keys():
        distances = {}
        for nbr in genetic_distances[patient].keys():
            genetic_pair_dist = genetic_distances[patient][nbr]
            epi_pair_dist = epi_distances[patient][nbr]
            distances[nbr] = genetic_pair_dist+(999*(epi_pair_dist/max)/100000)
        graph[patient]=distances
    return graph
#print(construct_complete_weighted_digraph("patient_sequences.txt",'patient_traces.txt'))
rdmst=infer_transmap("patient_sequences.txt",'patient_traces.txt',1)
#print(rdmst)
graph=rdmst[0]
#print(graph)
#graph={1: {2: 1.1, 6: 1.0, 4: 0.0}, 2: {6: 1.0}, 6: {}, 4: {2: 16.1, 6: 16.0}}
g = Digraph('G', filename='graph')
for node in graph:
    g.node(str(node))
for node, edges in graph.items():
    for edge, weight in edges.items():
        g.edge(str(node), str(edge), label=str(weight))
g.view()
#print(compute_genetic_distance('00101','10100'))
#print(compute_genetic_distance('00101','00101'))
#print(compute_genetic_distance('11010','00101'))