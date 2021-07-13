function [F_c, F_p] = HHD(W,varargin)
% This function returns the decomposition of a the edge weights (considered 
% as flows) of a weighted directed network (W) into a directed part and a 
% circulating part: flows in the circulating part are perfectly balanced 
% (in and out-flow are balanced at every node); the directed part (the 
% 'potential' or 'gradient' flow network) preserves the node imbalances (i.e. 
% node imbalance, our Eq.3 in [1]) and is a directed ascylic graph. Neither part 
% introduces any edges not existing in the original network (W). Note however 
% that  while new edges are not introduced some edge flows may be negative 
% (negative flows can also be interpret as reversal of direction), thus 
% others may have C(i,j)>w(i,j) or P(i,j)>w(i,j).
%
% This is the network Helmholtz-Hodge decompossition (see e.g. [3] and appendix in [1]).
%
% For true circular/directional flow decomposition, see [2] and function CDFD().
%
% Inputs:
%  W : an adjacency matrix (can be weighted)
%
% Outputs:
%  F_c : this is the circular flow network (flows are berfectly balanced)
%  F_g : this is the 'gradient' flow network (net equivalent with W and a directed asyclic graph)
%
% Dependencies
%   - levels (TC toolbox)
%   - parseArgs (TC toolbox)
%
% REFERENCES
% - [1] MacKay RS, Johnson S, Sansom B. (2020) How directed is a directed network? 
%          R. Soc. Open Sci. 7: 201138. http://dx.doi.org/10.1098/rsos.201138
% - [2] Sansom, B, MacKay, RS, Johnson, S, Zhou, Y (2021) True circular and directional
%                  flow decomposition (forethcoming).
% - [3] Kichikawa et al. (2019) Community structure based on circular
%                   flow in a large-scale transaction network, 8, Applied Network Science

% Contact: bazil.sansom@warwick.ac.uk

% ------- FUNCTION BEGGINS ------------------------------

% Default options:
opts=struct('flows','pos');     % negative flows interpreted as reversal 
                                %   of flow direction (all edges possitive). Alternative option is 'basic' (in which case this adjustment not made).
opts=parseArgs(varargin,opts);  % parse optional inputs

h=levels(W);      % obtain levels (equivalent to Hodge potentials, see [1])
h_dif=meshgrid(h)-meshgrid(h)'; % obtain level differences (this is all to all, but we will only consider h_dif over edges in the nework)

if strcmp(opts.flows,'basic') || strcmp(opts.flows,'pos')
    F_p=W.*h_dif; % obtain gradient flow network based on trophic level differences
    F_c=W-F_p; % circular flow network is then the difference between total flow 
           % and gradient flow
    
    % Function will either return the raw or "basic" decomposition, or
    % interprete negative flows as reversal of edge direction
    
    if strcmp(opts.flows,'basic')
        F_p=round(F_p,5); F_c=round(F_c,5); % rounding avoids very small numerical error based non-zero values
    
    elseif strcmp(opts.flows,'pos')
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
    end
else
        disp('ERROR: "flows" option must be specified as "basic" or "pos"')
end


end

