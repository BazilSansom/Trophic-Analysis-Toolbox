
**Contents**
-
- **Examples** provides some examples that illustrate methods/test functions running.

**Core functions:**
- **incoherence**   returns ***trophic levels*** (provides a wrapper for *levels*) and ***trophic incoherence*** statistics as per [1]
- **levels**        returns *trophic levels* as per [1]
- **standard_levels** returns trophic levels with fixed basal nodes as per [4]
- **HHD**           returns network Helmholtz-Hodge Decomposition based on trophic levels (see appendix of [1])

**Digraph opperations:**
- **adj2edgelist**  converts adjacency matrix to edgelist
- **edge_diff**     takes diferences over each edge (i->j or j->i) for vector of node values

**Auxiliary functions:**
- **parseArgs**    is a borrowed function that use extensively to allow optional inputs in many functions (from Aslak Grinsted [here](https://uk.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot)).

**References:**
-

[1] MacKay RS, Johnson S, Sansom B. 2020 How directed is a directed network? *R. Soc. Open Sci.* **7**: 201138.
http://dx.doi.org/10.1098/rsos.201138 (Paper available [here](https://doi.org/10.1098/rsos.201138)).

[2] Johnson S, Domínguez-García V, Donetti L and Muñoz MA (2014) "Trophic coherence determines food-web stability", PNAS, 111 (50), 17923-17928. https://doi.org/10.1073/pnas.1409077111

[3] Sansom,B, Johnson S, and MacKay RS (2021) "Trophic incoherence drives sytemic risk in financial exposure networks", Rebuilding Macroeconomics Working Paper Series, Working Paper No.39 (Paper available [here](https://warwick.ac.uk/fac/sci/maths/people/staff/sansom/sansom_2021.pdf)).

[4] Levine S. (1980) Several measures of trophic structure applicable to complex food webs. J. Theor. Biol. 83, 195–207. (doi:10.1016/0022-5193(80)90288-X)
