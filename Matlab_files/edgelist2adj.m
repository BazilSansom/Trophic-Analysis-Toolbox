function [adj]=edgelist2adj(EdgeList)

% This function takes an edgelist (numeric or string IDs) and converts to a
% numeric adjacency matrix format.

% INPUTS
% EdgeList
% 
% DEPENDENCIES
% tickers2numbers (TC toolbox)

if ~isnumeric(EdgeList)
    [~,EdgeList]=tickers2numbers(EdgeList);
end

nVar=max(max(EdgeList(:,1),max(EdgeList(:,2))));
adj=zeros(nVar,nVar);

if size(EdgeList,2)==2
    for i=1:size(EdgeList,1)  
        adj(EdgeList(i,1),EdgeList(i,2))=1;
    end
elseif size(EdgeList,2)==3
    for i=1:size(EdgeList,1)  
        adj(EdgeList(i,1),EdgeList(i,2))=EdgeList(i,3);
    end
else
    error('Input does not seem to be edge list')
end

 
end

