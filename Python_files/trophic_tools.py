'''
This Python module is for conducting trophic analysis as introduced by MacKay, Johnson and Sansom in [1], and also provides tools developed for network visualisation using these methods.
Equation numbers refered to in annotation are from [1].

If you make use of code provided here please site [1].

Functions included are:
- trophic_levels : returns trophic levels as per [1]
- trophic_incoherence : returns trophic incoherence as per [1] (where trophic coherence is 1-incoherence)
- trophic_layout : returns a layout where y-possitions are given by trophic levels, and x-possitions based on a modified force-directed graph drawing algorithm to spread nodes out on the x-axis. For reproducibility, user can save and specify seed.
- trophic_plot : plots network according to trophic_layout

Refs:
 [1] MacKay, Johnson & Sansom (2020), "How directed is a directed network", Royal Society Open Science

Written by: Bazil Sansom
Contact: bazil.sansom@warwick.ac.uk

Contributions welcome!

'''

# DEPENDENCIES for pip install
#from scipy.sparse.linalg import spsolve, cg
#from scipy.sparse import diags, lil_matrix
#import networkx as nx 
#import matplotlib.pyplot as plt #is this needed?
#import numpy as np
#import random2
#np.random.seed(4812)

# NOTES: this is written to work with networkx, but can easily be modified otherwise.

import numpy as np
import networkx as nx
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
from scipy.sparse import diags, lil_matrix
import matplotlib.pyplot as plt

#### --- TROPHIC ANALYSIS --- #####

# Define function to obtaion trophic levels

def trophic_levels(G):
    '''
    This function takes networkx graph object as input, and calculates trophic levels
    INPUTS:
        G networkx graph object
    OUTPUTS:
        h array of trophic levels
    '''
    # FUNCTION:

    # Checks:
    if not isinstance(G, nx.Graph):
        msg = "This function takes networkx graph object."
        raise ValueError(msg)
        
    # Check network weakly connected
    G2 = G.to_undirected(reciprocal=False,as_view=True)
    if nx.is_connected(G2):
        W = nx.adj_matrix(G)                  # weighted adjacency matrix as Compressed Sparse Row format sparse matrix of type '<class 'numpy.longlong'>'
        in_weight = W.sum(axis=0).A1          # Eq.2.1
        out_weight = W.sum(axis=1).A1         # Eq.2.1
        u = in_weight + out_weight            # (total) weight Eq.2.2
        v = in_weight - out_weight            # imbalance Eq.2.3 (the difference between the flow into and out of the node)
        L = diags(u, 0) - (W + W.transpose()) # (weighted) graph-Laplacian operator Eq.2.5
        L[0,0 ]= 0
        h = spsolve(L, v)                     # solve Eq.2.6 for h using sparse solver spsolve
        h = np.round(h - h.min(),decimals=10) # put lowest node at zero and remove numerical error
        return h
    else:
        # Should extend to identify components and obtain trohpic levels 
        # for each of these. However layout function currently only takes single 
        # connected component so this needs extending first.
        msg = 'Network must be weakly connected.'
        raise ValueError(msg)
        

# Define function to obtaion trophic incoherence

def trophic_incoherence(G):
    ''' 
    This function takes networkx graph object as input, and returns 
    trophic levels and coherence
    INPUTS
      G networkx graph object
    OUTPUTS:
      F_0 trophic incoherence
      h   array of trophic levels
    '''
    # FUNCTION
    
    h = trophic_levels(G)
    W = nx.adj_matrix(G)
    hj, hi = np.meshgrid(h, h)
    H=np.power([hj-hi-1],2)
    F_0 = (W.multiply(H)).sum() / W.sum()
    
    return F_0, h


# Define function to obtaion network Helmholtz-Hodge Decomposition from trophic levels

def HHD(G):
    ''' 
    This function takes networkx graph object as input, and returns 
    network HHD: the HHD enables us to break down the flow on a directed network into 
    two flow components: potential flow (also directonal or gradient flow) and 
    circular flow. The potential flow between a pair of nodes is given by the difference 
    of their potentials obtained by the Helmholtz-Hodge decomposition [ref] or 
    equivalently trophic-levels (shown in [1]) and runs from a node with higher potential 
    to a node with lower potential. On the other hand, the circular flow component is 
    balanced so contributes nothing to node imbalance vector w_in-w_out (Eq.2.3 (the 
    difference between the flow into and out of each node)).
    
    INPUTS
      G networkx graph object
    OUTPUTS:
      F_c circular flow network
      F_p potential flow network
    '''
    # FUNCION
    
    # Get trophic levels
    h = trophic_levels(G)
    
    # Obtain flow decompossition F=F_c+F_p
    W = nx.adj_matrix(G).todense()
    hj, hi = np.meshgrid(h, h)
    H=hj-hi
    F_p = np.multiply(W ,H)
    F_c = W - F_p
    
    # Directed network with possitive sign ordered pairs
    for i in range(len(W)):
        for j in range(len(W)):
            if F_p[i,j]<0:
                F_p[j,i]=F_p[j,i]+abs(F_p[i,j])
                F_p[i,j]=0
    for i in range(len(W)):
        for j in range(len(W)):
            if F_c[i,j]<0:
                F_c[j,i]=F_c[j,i]+abs(F_c[i,j])
                F_c[i,j]=0
                
    F_c=F_c+0
    F_p=F_p+0 
    
    return F_c, F_p

### -- VISUALISATION  -- ####

# Define a function to spread nodes on x-axis taking y-possitions as given.

def trophic_layout(G, 
                   k=None, 
                   ypos=None, 
                   iterations=50, 
                   seed=None,    
                   threshold=1e-4):
    ''' 
    This function position nodes in network G using modified Fruchterman-Reingold layout. 
    The layou spreads the nodes on the x-axis taking y-possitions as given. 
    By default the function uses for y-possitions trophic levels as defined in [1], 
    but may also be passed any y-possitions user chooses to define.
    
    REQUIRED INPUTS:
    G   : networkx graph object.
         Positions will be assigned to every node in G.
    
    OPTIONAL INPUTS
    k   : integer or None (default=None). If None the distance is set to 1/sqrt(nnodes) 
        where nnodes is the number of nodes.  Increase this value to spread nodes farther apart on x-axis .
    ypos: array or None (default=None). Initial y-positions for nodes. If None, then use
        trophic levels as defined by [1]. Alternatively user can pass any desired possitions as ypos. 
    iterations : integer (default=50)
        Maximum number of iterations taken.
    seed : integer, tuple, or None (default=None). For reproducible layouts, specify seed for random number generator.
           This can be an integer; or the output seedState obtained from a previous run of trophic_layout or trophic_plot in order to reporduce
           the same layout result.
           If you run:
               pos1, seedState = trophic_layout(G)
               pos2, _ = trophic_layout(G,seed=seedState)
            then pos1==pos2
    threshold: float (default = 1e-4)
        Threshold for relative error in node position changes.
        The iteration stops if the error is below this threshold.
    
    OUTPUTS
        pos a dictionary of possitions for each node in G
        seedState : (tuple) the seedState needed to reproduce layout obtained (e.g. if you run:
                        pos1, seedState = trophic_layout(G)
                        pos2, _ = trophic_layout(G,seed=seedState)
                  then pos2==pos1
    
    '''
    
    if seed is None or isinstance(seed,int):
        np.random.seed(seed) # sets seed using default (None) or user specified (int) seed
        seedState = np.random.get_state() # saves seed state to be returned (for reproducibility of result)
    elif isinstance(seed,tuple): # allows user to pass seedState obtained from previous run to reproduce same layout
        np.random.seed(None)
        np.random.set_state(seed)
        seedState=seed
    else:
        msg = '"Seed should be None (default); integer or tuple (use seedState which is output of trophic_layout).'
        raise ValueError(msg)
        
    import networkx as nx
    
    # Check networkx graph object
    if not isinstance(G, nx.Graph):
        msg='This function takes networkx graph object'
        raise ValueError(msg)
    
    # Check network weakly connected
    G2 = G.to_undirected(reciprocal=False,as_view=True)
    if not nx.is_connected(G2):
        msg='Network must be weakly connected'
        raise ValueError(msg)
        
    A = nx.to_numpy_array(G)
    dim=2
    nnodes, _ = A.shape
    
    A=A+np.transpose(A) # symmetrise for layout algorithm

    if ypos is None:
        h=trophic_levels(G)
        pos = np.asarray(np.random.rand(nnodes, dim), dtype=A.dtype)
        pos[:,1]=h
        # random initial positions
        #pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
        pos = np.asarray(np.random.rand(nnodes, dim), dtype=A.dtype)
        pos[:,1]=ypos
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)
        
    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)
    
    # the initial "temperature"  is about .1 of domain area
    # this is the largest step allowed in the dynamics.
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
    
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations + 1)
    
    for iteration in range(iterations):
        
        # matrix of difference between points
        delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
           # This is nnodes*nnodes*2 array giving delta_ij1=x_i-x_j delta_ij2=y_i-y_j
        
        # distance between points 
        distance = np.linalg.norm(delta, axis=-1) 
             # This is the nnodes*nnodes euclidian distance matrix with elements
             # d_ij=sqrt((x_i-x_j)^2+(y_i-y_j)^2)
                
        # enforce minimum distance of 0.01
        np.clip(distance, 0.01, None, out=distance)
        
        # displacement "force"
        displacement = np.einsum('ijk,ij->ik',
                                 delta,
                                 (k * k / distance**2 - A * distance / k))
        
        # update positions
        length = np.linalg.norm(displacement, axis=-1) # returns the Euclidean norm of each row in displacement (this norm is also called the 2-norm, or Euclidean LENGTH).
        length = np.where(length < 0.01, 0.01, length)  # enforce: where length<0.01, replace with 0.01, otherwise keep length
        delta_pos = np.einsum('ij,i->ij', displacement, t / length)
        delta_pos[:,1]=0 
        pos += delta_pos 
        
        # cool temperature
        t -= dt
        
        # check convergence
        err = np.linalg.norm(delta_pos) / nnodes
        if err < threshold:
            break
            
    pos = dict(zip(G, pos))
    return pos, seedState


# Define a function to plot network using trophic_layout

def trophic_plot(G,
                 k=None, 
                 ypos=None,
                 iterations=50, 
                 seed=None,
                 threshold=1e-4,
                 title='',
                 node_color=[[0, 0.451, 0.7412]]):
    '''
    This is just a wrapper for trophic_layout that automates some plotting decissions.
    Has same options and defaults as trophic_layout.
    Will add some plotting options.
    
    REQUIRED INPUTS:
    
    G   : networkx graph object.
         Positions will be assigned to every node in G.
    
    OPTIONAL INPUTS
    
    Layout options:
    k   : integer or None (default=None). If None the distance is set to 1/sqrt(nnodes) 
        where nnodes is the number of nodes.  Increase this value to spread nodes farther apart on x-axis .
    ypos: array or None (default=None). Initial y-positions for nodes. If None, then use
        trophic levels as defined by [1]. Alternatively user can pass any desired possitions as ypos. 
    iterations : integer (default=50)
        Maximum number of iterations taken.
    seed : integer, tuple, or None (default=None). For reproducible layouts, specify seed for random number generator.
           This can be an integer; or the output seedState obtained from a previous run of trophic_layout or trophic_plot in order to reporduce
           the same layout result.
           If you run:
               pos1, seedState = trophic_layout(G)
               pos2, _ = trophic_layout(G,seed=seedState)
            then pos1==pos2
    threshold: float (default = 1e-4)
        Threshold for relative error in node position changes.
        The iteration stops if the error is below this threshold.
    
    Plotting options (a selection of draw_networkx options):
        title : (str) otionally provide a title for the chart
        node_color : color string, or array of floats, (default RGB triplet [0, 0.451, 0.7412])
                     Node color. Can be a single color format string, or a  sequence of colors with
                     the same length as nodelist. If numeric values are specified they will be mapped 
                     to colors using the cmap and vmin,vmax parameters.  See matplotlib.scatter for 
                     more details.
    
    OUTPUTS
        plots input network G
        seedState : (tuple) the seedState needed to reproduce layout obtained (e.g. if you run:
                        seedState = trophic_plot(G)
                        then this exact plot can be replicated by running
                        _ = trophic_plot(G,seed=seedState)
                  
     '''
    
    F_0,_ = trophic_incoherence(G)
    pos, seedState = trophic_layout(G,k=k,ypos=ypos,iterations=iterations,threshold=threshold, seed=seed)
    fig, ax = plt.subplots(figsize = (10,10))
    nnodes=G.number_of_nodes()
    scale=1/nnodes
    nx.draw_networkx(G, with_labels=False,pos=pos,node_size=scale*2000,arrowsize=scale*150,width=scale*10, ax=ax, node_color=node_color);
    limits=plt.axis('on') # turns on axis
    ax.tick_params(left=True, labelleft=True)
    plt.ylabel('Trophic Levels')
    plt.xlabel('Trophic coherence = ' + "{:.2f}".format(1-F_0))
    font = { "fontweight": "bold", "fontsize": 16}
    ax.set_title(title,font)
    
    return seedState


