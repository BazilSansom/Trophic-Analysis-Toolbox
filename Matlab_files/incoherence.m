function [h,F0] = incoherence(W,varargin)

% This function calculates generalised heights and incoherence 
%  statistics as introduced in ref [1] for any weakly connected directed
%  network.

% INPUTS:
% Required
% W      adjacency matrix for directed graph (may be weighted or unweighted
%        and can include self-loops. Allows at most one edge from a node 
%        m to a node n (function assumes no multi-edges but might check for 
%        this and amalgamate in future)). This graph should be weakly connected
%        (in the sense that it would be connected if undirected).
% Optional:
% 'h0'  (optional) allows specification of how to set zero height. Default is to set
%       height of lowest node h=0; can also specify 'h0', 'wm' in order to
%       set the (weighted) degree weighted mean level to zero; or 'sm' to
%       centre around simple mean. Fmin is invarient to any shift by an
%       arbitrary constant, but this choice of zero level may be relevant
%       for some comparisons, or presentational issues when it comes to
%       heights.
%
% Outputs:
% - h    This is vector of heights for each node in W obtained based on 
%        improved notion of trophic level introduced by [1] (solution to (Eq.6)).
%        Standardised around some zero level (defined by 'h0' or default by
%        which lowest node has level zero). These heights minimize 
%        cost function ('trophic confusion') F(h) Eq.8 [1].
% - F0   This is the cost function F(h) (Eq.8 [1]) evaluated at h. This is
%        our measure of incoherence. F0 = 0 iff all the level differences 
%        zmn = hn − hm are 1; and F0 = 1 iff all the level differences are 0.
%        Trophic coherence is thus defined as 1-F0.
% - eta  the ratio of the edge-weighted mean height difference (which
%        evaluates to 1-F0) over the standard deviatio of the height
%        differences, which evaluates to sqrt(F0/(1-F0)). ([1] Eq.16).
%
% DEPENDENCIES:
% - levels              (TC toolbox)
% - adj2edgelist        (TC toolbox)
% - edge_diff           (TC toolbox)
%                       (have made this stand alone function as used accross 
%                       other functions and may be useful general opperation for DG).
% - parseArgs           (TC toolbox)
%
% References
% - [1] MacKay, Johnson and Sansom (2020) How directed is a directed
%         network? (Paper available at: https://arxiv.org/pdf/2001.05173.pdf).
% 
%Contact: bazil.sansom@warwick.ac.uk

%%%%%%%%%%%%%%   FUNCTION BEGGINS     %%%%%%%%%%%%%%%

% Optional inputs 
%(this approach gives some flexibility going forward)

opts=struct('h0','min');      % default options
opts=parseArgs(varargin,opts);% optional inputs

% - CHECKS -

% Test weakly connected
if ~isConnected(((W+W')>0))
    disp('ERROR: Heights only defined for a weakly connected network. Consider using GC option?')
    return 
end


% Test if directed
if isSymetric(W) 
    disp('WARNING: Symetric adjacecnyc matrix - newtwork may be un-directed?')
end

% - COMPUTE STATISTICS -

h=levels(W,'h0',opts.h0);            % obtain trophic levels ([1] Eq.6)
H=(meshgrid(h)-meshgrid(h)'-1).^2;
F0=sum(sum((W.*H)))/sum(sum(W));     % F0=F(h) ([1] Eq.7)
eta=(F0/(1-F0))^(1/2);               % ([1] Eq.16) 

end


% - SUBFUNCTIONS - 

% isConnected: tests whether 'weakly connected' defined as:
%              weakly connected iff connected as an undirected graph
% isSymetric:  tests whether (adjacecy) matrix symetric (as possible
%              indication of undirected graph)

function S = isConnected(adj)
% Check whether nework connected (for undirected newtork)
    
    if not(isempty(find((sum(adj)==0),1)))
        S = false; 
        return
    end
    
    n = length(adj);
    x = [1; zeros(n-1,1)]; % [1,0,...0] nx1 vector 
    
    while 1
        y = x;
        x = adj*x + x;
        x = x>0;
        
        if x==y
            break
        end
    end
    S = true;
    if sum(x)<n
        S = false;
    end
end


function S=isSymetric(adj)
% Check whether (adjacency) matrix is symetric
    S = false;
    if adj==transpose(adj)
        S = true; 
    end
end
