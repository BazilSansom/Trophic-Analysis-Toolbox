function el=adj2edgelist(adj)

% Inputs
% - adj  adjacency matrix (weighted or unweighted)
% Outputs
% - Edge list (colums one and two of el are s->t)
% - Edge weights (column three of el is weights as specified in adj)

n=length(adj); % number of nodes
edges=find(adj>0); % indices of all edges (only allows possitive edge weights)
%edges=find(adj~=0); % indices of all edges (also allows negative edge weights)

el=[];
for e=1:length(edges)
  [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  (Convert linear indices to subscripts)
  el=[el; i j adj(i,j)];
end
