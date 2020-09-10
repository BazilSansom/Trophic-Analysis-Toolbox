function [dif] = edge_diff(levels,DG,format,dir)

% This function computes differences on edges in a directed netwrok (DG)
% based on a vector of magnitudes given by 'levels'. This can be thought of
% as a gradient or potential difference accross each edge.

% INPUTS:
%  levels   is a vector assigning a 'height' or magnitude to each node in
%           DG. This may be output obtained using standard_level or
%           generalised_levels.
% DG        a directed network (as adjacency matrix or edgelist - need to specify as 'format')
% format    specify as 'adj' or 'edgelist' depending on input format.
% dir       need to specify 'ji' or 'ij'

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

nEdges=size(edgelist,1);          % simple number of edges
dif=zeros(size(edgelist,1),1); % prealocate vector for pairwise differences on each edge

if strcmp(dir,'ji')
    for edge=1:nEdges
        dif(edge,1)=(levels(edgelist(edge,2),1)-levels(edgelist(edge,1),1));
    end
elseif strcmp(dir,'ij')
    for edge=1:nEdges
        dif(edge,1)=(levels(edgelist(edge,1),1)-levels(edgelist(edge,2),1));
    end
else
    error('Must specify dir as "ji" or "ij"')
end


end

