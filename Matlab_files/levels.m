function [h] = levels(W,varargin)
% Obtains generalised heights as introduced in [1] (Eq.6). 
% All equation numbers refer to our paper [1].

% INPUTS
% Required inputs
% -  W    adjacency matrix for directed graph (may be weighted or unweighted
%         and can include self-loops. Allows at most one edge from a node 
%         m to a node n (function assumes no multi-edges but might check for this)).
% Optional inputs (name value pairs):
% - 'h0' allows specification of how to set zero height. Default is to set
%        height of lowest node h=0; can also specify 'h0', 'wm' in order to
%        set the (weighted) degree weighted mean level to zero; or 'h0','sm' to
%        centre around simple mean.
%
% Outputs
% - h    This is vector of heights for each node in W obtained based on 
%        improved notion of trophic level introduced by [1] (solution to (Eq.6)).
%        Standardised around some zero level (defined by 'h0' or default by
%        which lowest node has level zero).

% CHECKS

% (1) Weak connectedness (defined as weakly connecged iff connected when
%     treated as undirected).
% (2) Directed network (can't actually be checked for, but symetric matrix
%     taken as possible indication of undirected network.

% DEPENDENCIES
% - parseArgs (TC Toolbox).

% REFERENCES
% - [1] MacKay, Johnson and Sansom (2020) How directed is a directed
%         network? (Paper available at: https://arxiv.org/pdf/2001.05173.pdf).
% 
%Contact: bazil.sansom@warwick.ac.uk

% ------- FUNCTION BEGGINS ------------------------------

% Parse optional inputs

opts=struct('h0','min');       % default options
opts=parseArgs(varargin,opts); % optional inputs

% Checks
% (1) Weak connectedness
if ~isConnected(((W+W')>0))
    disp('ERROR: Heights only defined for a weakly connected network. Consider using GC option?')
    return 
end

% (2) Directed network
if isSymetric(W) 
    disp('WARNING: Symetric adjacecnyc matrix - newtwork may be un-directed?')
end

% Calculations

% (1) Setup
k_out=sum(W,2); k_in=sum(W,1)'; % weighted out and in degrees ([1] Eq.1)
u=k_in+k_out;                   % node weight   ([1] Eq.2)
v=k_in-k_out;                   % node imblance ([1] Eq.3)
lambda=diag(u)-W-W';            % our laplacian ([1] Eq.6)

% (2) Solve lambda*h=v for h
lambda(1,1)=0;   % remove arbitrariness of levels
%h=(lambda^-1)*v;
h=linsolve(lambda, v);

% (3) After solving this system subtract min(h) from all h(i), 
%     yielding the solution h which has min(h)=0, or define other zero
%     level.

if strcmp(opts.h0,'min')
    h=h-min(h);
    h=round(h,100);
elseif strcmp(opts.h0,'wm')
    h = h-sum(u.*h,1)./sum(u,1); % where u total weighted degree
    h=round(h,100);
elseif strcmp(opts.h0,'sm')
    h = h-mean(h);
    h=round(h,100);
else
    disp('ERROR: h0 incorrectly specified')
end

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

