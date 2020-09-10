function [trophic_levels,q,q_hat] = incoherence_stand(DG,levels)

% This function computes the standard trophic levels according to [1] (Eq.2 in [3]) and 
% trophic coherence statistic (q) introduced by [2] (Eq.3 in [3]) as well
% as q relative to null expectation.
%
% INPUTS:
%    A        a (possibly weighted) adjacency matrix.
%    levels  (Optional input) vector of heights to be used.
% OUTPUTS:
%    trophic_levels   Is a vector n*1 of trophic levels according to
%                     Eq.2 [3] where n is the number of nodes.
%    q                Is 'trophic incoherence' statistic Eq.3 [3]. q is 0 
%                     for perfectly coherent networks, and positive otherwise. 
%                     We therefore refer to q as 'trophic incoherence' [2].
%    q_hat            This is q/q_tilda where q_tilda (Eq.6 [4]) provides
%                     null expectation. q_hat thus close to one when q
%                     close to random expectaion. Values lower than 1 reveal 
%                     coherent networks, while values greater than 1 
%                     reveal incoherent ones (see ).
%
% Dependencies:
%     standard_levels   (TC toolbox)
%     compute_edge_diff (TC toolbox)
%     ewmean            (TC toolbox)
%
% REFERENCES:
% [1] Levine (1980) JTB Vol 83
% [2] Johnson, S., Domínguez-García, V., Donetti, L., & Muñoz, M. A. (2014). 
%     Trophic coherence determines food-web stability. 
%     Proceedings of the National Academy of Sciences, 111(50), 17923?17928.
% [3] Our paper....
% [4] Johnson and Jones (2017) '?Looplessness in networks is linked to
%       trophic coherence'
% [5] Pagani et al. (2019) '?Resilience or robustness : identifying
%         topological vulnerabilities in rail networks'
%
% Example:
% [provide simple examples]

% CHECKS and setup

% Set up to take DG as adjacency matrix or as edgelist

if ismatrix(DG) && (size(DG, 1) == size(DG, 2))
    adj_mtrx=DG;
else
    adj_mtrx=edgelist2adj(DG);
end

% CALCULATE

if nargin==1
    try
        % Compute trophic levels from adjacency matrix
        [trophic_levels] = standard_levels(adj_mtrx);
    catch
        return
    end
elseif nargin==2
    trophic_levels=levels;
end

% Compute trophic differences from trophic levels
[~,trophic_diff] = compute_edge_diff(trophic_levels,adj_mtrx,'adj');

% Compute coherence (q) from trophic differneces
q=sqrt(ewmean(trophic_diff.^2,adj_mtrx)-1);

% Compute q relative to null model (basal ensemble expectation q_tilda)
% q/q_tilda
[q_tilda]=basal_ensemble_expectation(adj_mtrx); % Obtain basel ensemble expectation q_tilda
q_hat=q/q_tilda;                                % calculate q/q_tilda

end

function [q_tilda] = basal_ensemble_expectation(A)

L=sum(sum(A));         % total weighted number of edges
k_out=sum(A,2);        % weighted out degree
k_in=sum(A,1)';        % weighted in degree
source=k_in==0;        % find source nodes
Li=sum(k_out(source)); % Count number of edges connected to sources ('basal nodes')

q_tilda=sqrt(L/Li-1);  % corresponds to Eq.6 [4].

end
