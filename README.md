This repository implements Yamada and Kinoshita's algorithm for negative edge-weight cycle enumeration, [published in](https://www.sciencedirect.com/science/article/pii/S0166218X01002013) Discrete Applied Mathematics.

The standard Bellman–Ford algorithm only _detects_ a negative cycle; with some modification (implemented in [bellmanford.py](bellmanford.py)), it can be made to report one. Yet this reported cycle is arbitrary among the graph's negative cycles; in particular, it need not be the "most negative" one. Meanwhile, enumerating _all_ cycles (by say [Johnson's algorithm](https://www.cs.tufts.edu/comp/150GA/homeworks/hw1/Johnson%2075.PDF)) and filtering for negative ones is inefficient.

Yamada and Kinoshita's important algorithm fills this gap.

### Thanks:

Our [implementation of](bellmanford.py) Bellman–Ford (used as a subroutine in Yamada–Kinoshita) is based on [Nelson Uhan's](https://github.com/nelsonuhan/bellmanford) work.

## Usage:
```python
import networkx as nx
import yamadakinoshita as yk
G = nx.DiGraph()
G.add_edge('a', 'b', weight=1)
G.add_edge('b', 'c', weight=6)
G.add_edge('c', 'a', weight=4)
G.add_edge('b', 'd', weight=5)
G.add_edge('d', 'a', weight=4)
list(yk.negative_edge_cycles(G))
# []

G['b']['c']['weight'] = -6
list(yk.negative_edge_cycles(G))
# [(-1, [('c', 'a'), ('a', 'b'), ('b', 'c')])]

G['d']['a']['weight'] = -10
list(yk.negative_edge_cycles(G))
# [(-4, [('a', 'b'), ('b', 'd'), ('d', 'a')]), (-1.0, [('a', 'b'), ('b', 'c'), ('c', 'a')])]
```
