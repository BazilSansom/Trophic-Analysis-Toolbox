function [nodeIDs,edgeNums] = tickers2numbers(EdgeList,NodeList)
% This function takes a network with string based IDs and returns node list
% and edge list with numberic IDs. If a nodelist is provided, then IDs are
% assigned according to order of this nodelist, otherwise nodelist created
% from edgelist.

% To do: Should have flexibility to take weighted edgelist as input
% Should be called something 2 numeric..?

% If nodelist not provided creates one from EdgeList
if nargin == 0
    error('At least edgelist required')
elseif nargin == 1
    NodeList=uniqueRowsCA([EdgeList(:,1);EdgeList(:,2)]);
end

nodeIDs=(1:size(NodeList,1))';
edgeNums=zeros(size(EdgeList,1),size(EdgeList,2));

for i=1:size(EdgeList,1)
  for j=1:2  
   ind=strcmp(EdgeList(i,j),NodeList(:,1));
   edgeNums(i,j)=nodeIDs(ind);
  end 
end
    
end

