import bellmanford as bf
from collections import Counter

def negative_edge_cycles(G, weight='weight'):
    return _all_nc(G, [], [], weight)

def _all_nc(G, F, R, weight='weight'):  # F and R should be EDGE iterables of G i guess.
    if not F:
        z_upper, C = _an_nc(G, R, weight)
        pi_upper = C
        if not C:
            return
    else:  # step 2
        G_F_R = G.copy()
        v = [e[0] for e in F] + [F[-1][1]]
        G_F_R.remove_nodes_from(v[1:-1])
        G_F_R.remove_edges_from(R)
        G_F_R.remove_edges_from(list(G_F_R.out_edges(v[0])))  # errors if don't cast to list...!
        G_F_R.remove_edges_from(list(G_F_R.in_edges(v[-1])))  # i think probably an upstream bug.

        d_lower, pi_lower = bf.bellman_ford(G_F_R, v[-1], v[0], weight)
        counts = Counter(e[0] for e in pi_lower)
        c, count = max(counts.items(), key=lambda item: item[1])
        if count > 1:
            G_0_F_R = G_F_R.copy()
            G_0_F_R.remove_node(c)
            d_0, _ = bf.bellman_ford_k(G_0_F_R, v[-1], v[0], G_0_F_R.number_of_nodes(), weight)  # the tacit assumed - 1 needs to be overridden.
            G_1_F_R = G_F_R.copy()
            G_1_F_R.remove_edges_from(list(G_1_F_R.out_edges(c)))
            G_2_F_R = G_F_R.copy()
            G_2_F_R.remove_edges_from(list(G_2_F_R.in_edges(c)))
            d_1 = min(bf.bellman_ford_k(G_1_F_R, v[-1], c, k, weight)[0] + bf.bellman_ford_k(G_2_F_R, c, v[0], G_F_R.number_of_nodes() - 1 - k, weight)[0] for k in range(1, G_F_R.number_of_nodes() - 1))
            z_lower = sum(G.edges[u, v][weight] for u, v in F) + min(d_0, d_1)
        else:
            z_lower = sum(G.edges[u, v][weight] for u, v in F) + d_lower
        if z_lower >= 0:
            return
        d_upper, pi_upper = bf.dijkstra(G_F_R, v[-1], v[0])  # pi_upper is EDGES
        z_upper = sum(G.edges[u, v][weight] for u, v in F) + d_upper
        if z_upper < 0:
            C = F + pi_upper
        else:  # step 5
            for e in G_F_R.out_edges(v[-1]):
                yield from _all_nc(G, F + [e], R)
            return
    yield z_upper, C  # step 4
    for i, e in enumerate(pi_upper):
        yield from _all_nc(G, F + pi_upper[:i], R + [e])

def _an_nc(G, R, weight='weight'):  # suppress the argument F
    G = G.copy()  # revisit
    G.remove_edges_from(R)
    length, nodes, _ = bf.negative_edge_cycle(G, weight)
    return length, nodes
