from collections import defaultdict
from networkx.utils import generate_unique_node
from copy import copy
from heapq import heappush, heappop

def negative_edge_cycle(G, weight='weight'):
    """
    If there is a negative edge cycle anywhere in G, returns True.
    Also returns the total weight of the cycle and the nodes in the cycle.
    Parameters
    ----------
    G : NetworkX graph
    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight
    Returns
    -------
    length : numeric
        Length of a negative edge cycle if one exists, otherwise None.
    edges: list
        Edges in a negative edge cycle (in order) if one exists,
        otherwise None.
    negative_cycle : bool
        True if a negative edge cycle exists, otherwise False.
    Examples
    --------
    >>> import networkx as nx
    >>> import bellmanford as bf
    >>> G = nx.cycle_graph(5, create_using = nx.DiGraph())
    >>> print(bf.negative_edge_cycle(G))
    (None, [], False)
    >>> G[1][2]['weight'] = -7
    >>> print(bf.negative_edge_cycle(G))
    (-3, [(1, 2), (2, 3), (3, 4), (4, 0), (0, 1)], True)
    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.
    This algorithm uses bellman_ford() but finds negative cycles
    on any component by first adding a new node connected to
    every node, and starting bellman_ford on that node.  It then
    removes that extra node.
    """
    newnode = generate_unique_node()
    G.add_edges_from([(newnode, n) for n in G])

    try:
        pred, _, negative_cycle_end = _bellman_ford_relaxation(G, newnode, G.number_of_nodes(), weight)

        if negative_cycle_end:
            edges = []
            negative_cycle = True
            v = negative_cycle_end
            while True:
                u = pred[G.number_of_nodes()][v]
                edges.insert(0, (u, v))
                if edges.count((u, v)) > 1:
                    end_index = edges[1:].index((u, v))
                    edges = edges[:end_index + 1]
                    break
                v = u
            length = sum(G[u][v].get(weight, 1) for (u, v) in edges)
        else:
            edges = None
            negative_cycle = False
            length = None

        return length, edges, negative_cycle
    finally:
        G.remove_node(newnode)

def bellman_ford(G, source, target, weight='weight'):
    """
    Compute shortest path and shortest path lengths between a source node
    and target node in weighted graphs using the Bellman-Ford algorithm.
    Parameters
    ----------
    G : NetworkX graph
    pred: dict
        Keyed by node to predecessor in the path
    dist: dict
        Keyed by node to the distance from the source
    source: node label
        Source node
    target: node label
        Target node
    weight: string
       Edge data key corresponding to the edge weight
    Returns
    -------
    length : numeric
        Length of a negative cycle if one exists.
        Otherwise, length of a shortest path.
        Length is inf if source and target are not connected.
    edges: list
        List of edges in a negative edge cycle (in order) if one exists.
        Otherwise, list of nodes in a shortest path.
        List is empty if source and target are not connected.
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5, create_using = nx.DiGraph())
    >>> bf.bellman_ford(G, source=0, target=4)
    (4.0, [(0, 1), (1, 2), (2, 3), (3, 4)])
    """
    return bellman_ford_k(G, source, target, G.number_of_nodes() - 1, weight)

def bellman_ford_k(G, source, target, k, weight='weight'):
    # Get shortest path tree
    pred, dist, _ = _bellman_ford_relaxation(G, source, k, weight)

    edges = []

    v = target
    for i in range(k, 0, -1):
        u = pred[i][v]
        edges.insert(0, (u, v))
        if u == source or u is None:
            break
        v = u

    length = dist[-1][target]  # should even return length...?

    return length, edges  # not returning negative cycle i guess?


def _bellman_ford_relaxation(G, source, k, weight):
    """
    Relaxation loop for Bellman–Ford algorithm
    Parameters
    ----------
    G : NetworkX graph
    pred: list(defaultdict(node label))
        Keyed by node to predecessor in the path
    dist: list(defaultdict(float))
        Keyed by node to the distance from the source
    source: list
        List of source nodes
    weight: string
       Edge data key corresponding to the edge weight
    Returns
    -------
    pred, dist : dict
        Returns two dictionaries keyed by node to predecessor in the
        path and to the distance from the source respectively.
    negative_cycle_end : node label
        Backtrack from this node using pred to find a negative cycle, if
        one exists; otherwise None.
    """
    pred = [defaultdict(lambda: None)]
    dist = [defaultdict(lambda: float('inf'))]
    dist[0][source] = 0.0

    if G.is_multigraph():
        def get_weight(edge_dict):
            return min(eattr.get(weight, 1) for eattr in edge_dict.values())
    else:
        def get_weight(edge_dict):
            return edge_dict.get(weight, 1)

    for i in range(1, k + 1):
        pred.append(copy(pred[-1]))
        dist.append(copy(dist[-1]))
        for u, v in G.edges():
            dist_v = dist[i - 1][u] + get_weight(G[u][v])
            if dist_v < dist[i - 1][v]:
                dist[i][v] = dist_v
                pred[i][v] = u

    for u, v in G.edges():  # sketch to keep using i here...!
        if dist[-1][v] > dist[-1][u] + get_weight(G[u][v]):
            negative_cycle_end = v
            return pred, dist, negative_cycle_end

    negative_cycle_end = None
    return pred, dist, negative_cycle_end

def dijkstra(G, source, target, weight='weight'):
    paths = {source: []}
    dist = {source: 0.0}  # dictionary of NOT final distances
    seen = set()
    # fringe is heapq with 2-tuples (distance, node)
    fringe = []
    heappush(fringe, (0, source))
    while fringe:
        (d, u) = heappop(fringe)
        if u in seen:  # should only happen for duplicates
            continue  # already searched this node.
        seen.add(u)
        if u == target:
            break
        for u, v in G.out_edges(u):
            cost = G[u][v][weight]
            uv_dist = dist[u] + cost
            if v not in dist or uv_dist < dist[v]:
                dist[v] = uv_dist
                heappush(fringe, (uv_dist, v))  # duplicate entries
                paths[v] = paths[u] + [(u, v)]

    return dist[target], paths[target]
