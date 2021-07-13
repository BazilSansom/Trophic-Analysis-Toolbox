function [F_c, F_p] = mcHHD(W,varargin)
% This function is a wrapper for the function HHD (the network Helmholtz-Hodge decompossition 
% (see e.g. [2] and appendix in [1])) in order to provide flexibility to 
% compute HHD for a network with multiple connectedcomponents.
%
% Inputs:
%  W : a (weighted) adjacency matrix
%
% Outputs:
%  F_c : this is the circular flow network (flows are berfectly balanced)
%  F_g : this is the 'gradient' flow network (net equivalent with W and a directed asyclic graph)
%
% Dependencies:
%  HHD (TC toolbox)
%    |_ levels (TC toolbox)
%    |_ parseArgs (TC toolbox)
%
% REFERENCES
% - [1] MacKay RS, Johnson S, Sansom B. (2020) How directed is a directed network? 
%          R. Soc. Open Sci. 7: 201138. http://dx.doi.org/10.1098/rsos.201138
% - [2] Kichikawa et al. (2019) Community structure based on circular
%                   flow in a large-scale transaction network, 8, Applied Network Science

% Contact: bazil.sansom@warwick.ac.uk

% ------- FUNCTION BEGGINS ------------------------------

% Default options:
opts=struct('flows','pos');     % negative flows interpreted as reversal 
                                %   of flow direction (all edges possitive). Alternative option is 'basic' (in which case this adjustment not made).
opts=parseArgs(varargin,opts);  % parse optional inputs

% Identify weakly connected components (WCC):
G=digraph(W);
[bin,binsize] = conncomp(G,'Type','weak');
nodeID=1:size(W,1);

F_c=zeros(size(W)); % preallocate
F_p=zeros(size(W)); % preallocate

% Loop over all WCCs
for i=1:max(bin)
    idx = bin == i;
    id=nodeID(idx);
    if binsize(i)>1
        SG = subgraph(G, id);               % obtain WCC id as subgraph 
        A=full(adjacency(SG,'weighted'));   % create adjacency matrix representation
        [c,p]=HHD(A,'flows','basic');       % obtain HHD of WCC id
        
        % Write circular flows for WCC id to F_c
        for ki=1:size(A)
            for kj=1:size(A)
                if c(ki,kj)~=0 % why do I have this condition?
                   F_c(id(ki),id(kj))= c(ki,kj);
                end
            end
        end
        
        % Write to F_p
        for ki=1:size(A)
            for kj=1:size(A)
                if p(ki,kj)~=0  % why do I have this condition?
                   F_p(id(ki),id(kj))= p(ki,kj);
                end
            end
        end
        
    end
    
end

% Having obtained mcHHD, if negative edges to be interpreted as reversal of
% direction, make this adjustment

if strcmp(opts.flows,'pos')
    
    % Interprete negative flows as reversal of edge direction
        % potential flow matrix
        for i=1:size(W,1)
            for j=1:size(W,1)
                if F_p(i,j)<0
                    F_p(j,i)=abs(F_p(i,j))+F_p(j,i);
                    F_p(i,j)=0;
                end
            end
        end
        
        % circular flow matrix
        for i=1:size(W,1)
            for j=1:size(W,1)
                if F_c(i,j)<0
                    F_c(j,i)=F_c(j,i)+abs(F_c(i,j));
                    F_c(i,j)=0;
                end
            end
        end
        
        F_p=round(F_p,5); F_c=round(F_c,5); % rounding avoids very small numerical error based non-zero values

elseif strcmp(opts.flows,'basic')
    %F_p=F_p; F_c=F_c;

else
    F_p=[]; F_c=[];
    disp('ERROR: "flows" option must be specified as "basic" or "pos"')
end

end

