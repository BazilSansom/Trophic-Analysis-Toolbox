
**Contents**
-
- **Examples** provides some examples that illustrate methods/test functions running.

**Core functions:**
- **incoherence**   returns ***trophic levels*** (provides a wrapper for *levels*) and ***trophic incoherence*** statistics as per [1]
- **levels**        returns *trophic levels* as per [1]

**Digraph opperations:**
- **adj2edgelist**  converts adjacency matrix to edgelist
- **edge_diff**     takes diferences over each edge (i->j or j->i) for vector of node values

**Auxiliary functions:**
- **parseArgs**    is a borrowed function that use extensively to allow optional inputs in many functions (from Aslak Grinsted [here](https://uk.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot)).

**References:**
-

[1] MacKay, Johnson and Sansom (2020) How directed is a directed network?

Paper available [here](https://arxiv.org/pdf/2001.05173.pdf) and as Rebuilding Macroeconomic working paper [here](https://www.rebuildingmacroeconomics.ac.uk/how-directed-is-a-directed-network).
