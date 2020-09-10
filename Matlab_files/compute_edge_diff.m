function [dif_ji,abs_dif] = compute_edge_diff(levels,DG,format)

% This function computes differences on edges in a directed network (DG)
% based on a vector of node 'heights' given by 'levels'.

% INPUTS:
%  levels   is a vector assigning a 'height' or magnitude to each node in
%           DG. This may be output obtained using standard_level or
%           generalised_levels.
% DG        a directed network (as adjacency matrix or edgelist)
% format    specify as 'adj' or 'edgelist' depending on input format

% Function takes DG as adj or edgelist, but requires edgelist
if strcmp(format,'adj')
%if ismatrix(DG) && (size(DG, 1) == size(DG, 2))
    edgelist=adj2edgelist(DG);
elseif strcmp(format,'edgelist')
    edgelist=DG;
else
    error('Must specify DG format as "adj" or "edgelist"')
end

% compute the difference over each edge

%edgelist=adj2edge(adj);                 % convert adjacency to edgelist (indices)

nEdges=size(edgelist,1);               % simple number of edges
abs_dif=zeros(size(edgelist,1),1); % prealocate vector for pairwise differences on each edge
dif_ji=zeros(size(edgelist,1),1); % prealocate vector for pairwise differences on each edge

for edge=1:nEdges
    abs_dif(edge,1)=abs(levels(edgelist(edge,2),1)-levels(edgelist(edge,1),1));
    dif_ji(edge,1)=(levels(edgelist(edge,2),1)-levels(edgelist(edge,1),1));
end

end

