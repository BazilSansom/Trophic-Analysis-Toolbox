function [s] = standard_levels(adj_mtrx)
% This function calculates standard trophic levels (Eq.2 [1]) 
% according to [2].
%
% INPUTS:
%  A  a (possibley weighted) adjacency matrix
%
% OUTPUTS:
%  s   a vector of trophic levels according to standard definition
%      (corresponding to Eq.2 [1]).

% REFERENCES
% [1] Mackay, Johnson and Sansom ...
% [2] Levine (1980) JTB Vol 83
%
% Function beggins

% Construct v
nNodes=size(adj_mtrx,1); % number of nodes
k_in=sum(adj_mtrx,1);    % vector of in-degrees

%v=zeros(nNodes,1);% pre-allocating
%for i=1:nNodes 
%    v(i)=max([k_in(i),1]);
%end

v=max([k_in',ones(nNodes,1)],[],2); % basal nodes have level 1 by convention
%Lambda=diag(v)-transpose(adj_mtrx);
Lambda=diag(v)-adj_mtrx';
%Lambda=diag(v)-adj_mtrx;

if det(Lambda)==0
    %fprintf('Error: Lambda is a singular matrix.\n  Network must have at lease one basal node to obtain standard trophic levels.\n Consider using new levels.')
    error('Error: Lambda is a singular matrix.\n  Network must have at lease one basal node to obtain standard trophic levels.\n Consider using new levels.')
    s=[];
    % 
    return
else
    % Solve Ls=v for s
    %s=(Lambda^-1)*v; 
    s=Lambda\v;
end

end

