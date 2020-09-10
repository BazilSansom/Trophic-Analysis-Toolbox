function [x_out] = ewmean(x_in,A)

% This opperator takes average over all edges of a graph using possibly
% weighted edges.

% This corresponds to Eq.4 [1]

% Inputs:
%  - A     adjacency matrix (weighted or simple)
%  - x_in  the edge variable to be averaged (i.e. this needs to be an n*1 
%          vector where n is the simple number of edges.

% References:
%  - [1] (Our paper)

% FUNCTION

L=sum(sum(A));
E=adj2edge(A);
x_out=L^-1*sum(E(:,3).*x_in);

end

% SUBFUNCTIONS

function el=adj2edge(adj)

% Inputs
% - Adjacency matrix
% Outputs
% - Edge list (colums one and two are s->t)
% - Edge weights (column three is weights as specified in adj)

n=length(adj); % number of nodes
edges=find(adj>0); % indices of all edges

el=[];
for e=1:length(edges)
  [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  (Convert linear indices to subscripts)
  el=[el; i j adj(i,j)];
end

end

