function [s] = standard_levels(W)
% This function calculates standard trophic levels (Eq.2 [1]) 
% according to [2].
%
% INPUTS:
%  W  a (possibly weighted) adjacency matrix
%
% OUTPUTS:
%  s   a vector of trophic levels according to standard definition
%      (corresponding to solution to Eq.4.2 [1]).

% REFERENCES
% [1] MacKay RS, Johnson S, Sansom B. 2020 How directed is a directed network? 
%       R. Soc. Open Sci. 7: 201138. http://dx.doi.org/10.1098/rsos.201138
% [2] Levine (1980) JTB Vol 83
%
% Function beggins

% Construct v
nNodes=size(W,1); % number of nodes
k_in=sum(W,1);    % vector of in-weights

%v=zeros(nNodes,1);% pre-allocating
%for i=1:nNodes 
%    v(i)=max([k_in(i),1]);
%end

v=max([k_in',ones(nNodes,1)],[],2); % basal nodes have level 1 by convention
Lambda=diag(v)-W'; % Laplacian

if det(Lambda)==0
    %fprintf('Error: Lambda is a singular matrix.\n  Network must have at lease one basal node to obtain standard trophic levels.\n Consider using new levels.')
    error('Error: Lambda is a singular matrix.\n  Network must have at lease one basal node to obtain standard trophic levels.\n Consider using improved levels.')
    s=[];
    % 
    return
else
    % Solve Ls=v for s (Eq.4.2 [1])
    %s=(Lambda^-1)*v; 
    s=Lambda\v;
end

end

