function [A,dim] = adj2adjm(adjobj)
%
% Function to convert an adjacency object matrix in SYSDYNET format to the 
% adjacency matrix of the extended graph.
%
% Variables:
% adjobj: adjacency matrix object of a network;
%      
% A:   adjacency matrix of the extended graph. 
% dim: dimensions of the different signals [L p K]
%
%   Author:  Paul Van den Hof
%            Eindhoven University of Technology.
%   Version: 1.1 
%   Date:    14-2-2021, 15-2-2021
%            1-8-2021: adapted to ordering [G H R] and corresponding
%            dimensions [L p K]
%
    if nargin > 1
      err('number of input arguments is not correct')
    end  
%
    L = size(adjobj.G,1); % number of nodes
    p = size(adjobj.H,2); % number of white noises
    K = size(adjobj.R,2); % number of external excitations
    n = L+K+p;
%
% Build the extended graph of the network. 
% If no external excitation signals are present, or no noise signals are
% present then the corresponding graphs are not included.
    if K > 0 && p > 0
        A = [adjobj.G adjobj.H adjobj.R; zeros(K+p,K+p+L)];
    elseif K > 0
        A = [adjobj.G adjobj.R; zeros(K+p,K+p+L)];
    elseif p > 0
        A = [adjobj.G adjobj.H; zeros(K+p,K+p+L)];
    else
        A = adjobj.G; 
    end
%
    dim = [L p K];
%

