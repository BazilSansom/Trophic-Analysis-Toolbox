from scipy.sparse.linalg import spsolve, cg
from scipy.sparse import diags, lil_matrix
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np


#The functions below are able to handle sparse matrices and so they are a much more efficient memory implementation of the previous functions you had coded up
# Thus these are able to handle much larger networks then previously. 


def get_exact_trophic_levels(G, weight=None):
    G2 = G.to_undirected(reciprocal=False,as_view=True)
    if nx.is_connected(G2):
        B = nx.adj_matrix(G, weight=weight)
        in_deg = B.sum(axis=0).A1
        out_deg = B.sum(axis=1).A1
        v = in_deg - out_deg
        w = in_deg + out_deg
        L = diags(w, 0) - (B + B.transpose())
        L[0,0 ]= 0
        h = spsolve(L, v)
        h = h - h.min()
        return h
    else:
        print('Network must be weakly connected')
        
        
def get_exact_trophic_coherence(G, weight=None):
    h = get_trophic_levels(G, weight=weight)
    B = nx.adj_matrix(G, weight=weight)
    
    H = lil_matrix(B.shape, dtype=float)    
    for i, j in zip(B.nonzero()[0], B.nonzero()[1]):
        H[i,j] = h[j] - h[i] - 1
    
    H2 = (H.tocsr()).power(2)
    F_0 = (B.multiply(H2)).sum() / B.sum()
    
    return F_0, h




#Using one of the iterative methods for linear equation solvers instead of spsolve is much more efficient but not as accurate in terms of the exact solution but the overall rnaking is still correct

def get_trophic_levels(G, weight=None):
    G2 = G.to_undirected(reciprocal=False,as_view=True)
    if nx.is_connected(G2):
        B = nx.adj_matrix(G, weight=weight)
        in_deg = B.sum(axis=0).A1
        out_deg = B.sum(axis=1).A1
        v = in_deg - out_deg
        w = in_deg + out_deg
        L = diags(w, 0) - (B + B.transpose())
        L[0,0 ]= 0
        h = cg(L, v)
        h = h - h.min()
        return h
    else:
        print('Network must be weakly connected')
        
        
        
def get_trophic_coherence(G, weight=None):
    h = get_trophic_levels(G, weight=weight)
    B = nx.adj_matrix(G, weight=weight)
    
    H = lil_matrix(B.shape, dtype=float)    
    for i, j in zip(B.nonzero()[0], B.nonzero()[1]):
        H[i,j] = h[j] - h[i] - 1
    
    H2 = (H.tocsr()).power(2)
    F_0 = (B.multiply(H2)).sum() / B.sum()
    
    return F_0, h



G_semi_coherent=nx.DiGraph()
G_semi_coherent.add_edges_from([(1,4),(1,5),(2,5),(3,5),(3,6),(4,7),(5,7),(5,8),(5,9),
                  (6,9),(1,7),(2,8),(3,9),(1,8),(2,7),(2,9)])
G_coherent=nx.DiGraph()
G_coherent.add_edges_from([(1,4),(1,5),(2,5),(3,5),(3,6),(4,7),(5,7),(5,8),(5,9),(6,9)])

G_incoherent=nx.DiGraph()
G_incoherent.add_edges_from([(1,2), (2,3),(3,1)])



# Calculating MacKay Trophic levels for coherent network and then visualising the network with the trophic level controlling the coordinates of the nodes in the plot
G=G_coherent
F,h=get_trophic_coherence(G)
pos = nx.spring_layout(G)

trophic_levels = {}
for i, node in enumerate(G.nodes):
    trophic_levels[node] = h[i]
for node in G_coherent.nodes:
    pos[node][1] = trophic_levels[node]
    
plt.figure(figsize = (5,5))
nx.draw(G, pos=pos);
nx.draw_networkx_labels(G, pos=pos);
print(f'F_0 coherent: {round(F,4)}')
print(f'h: {h}')



G=G_semi_coherent
F,h=get_trophic_coherence(G)

pos = nx.spring_layout(G)
trophic_levels = {}
for i, node in enumerate(G.nodes):
    trophic_levels[node] = h[i]
for node in G_coherent.nodes:
    pos[node][1] = trophic_levels[node]
    
plt.figure(figsize = (5,5))
nx.draw(G, pos=pos);
nx.draw_networkx_labels(G, pos=pos);
print(f'F_0 semi-coherent: {round(F,4)}')
print(f'h: {round(h,3)}')



G=G_incoherent
F,h=get_trophic_coherence(G)
pos = nx.spring_layout(G)
plt.figure(figsize = (5,5))
nx.draw(G, pos=pos);
nx.draw_networkx_labels(G, pos=pos);
print(f'F_0 incoherent: {round(F,4)}')
print(f'h: {h}')



#Visualisation of graphs where h is the trophic levels

troph_positions = {}
for i in range (len(h)):
    
    troph_positions[i]= [np.random.random(),h[i]]
 

pos= troph_positions


fig, ax = plt.subplots(figsize=(5, 10))

#networkx drawing call 
nx.draw(G, pos, node_size=40, node_color='b', ax=ax)

# turn the axis on
ax.set_axis_on()
ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
plt.ylabel('Trophic Level')
plt.xlabel('Trophic incoherence = ' + "{:.2f}".format(F))
ax.xaxis.set_major_formatter(plt.NullFormatter())
plt.xticks([], [])



