function [Test] = modidfanalysis_partial(varargin,moduleindex,maxiter)
% Verify generic identifiability of a single module with partial
% measurement using a sufficient condition (Proposition 2 in the reference)
%
%
% [Test] = modidfanalysis_partial(varargin,moduleindex,maxiter)
% Function created on Matlab 2021a
%
% Note that this function only works when set Xj and C have small cardinality,
% since it possibly requires to evaluate over all the subsets of Xj and C.
% We can use 'maxiter' input to restrict the maximum iterations allowed for
% the subset search
%
% Inputs: varargin: 1-dimensional or 2-dimensional cell array, where 
%                   1. varargin{1} contains the adjacency matrix object of a network and must be specified; 
%                           a. varargin{1} is a structure consisting of
%                           fields G, R, H, C, which are binary matrices
%                           specifying the network topology
%                           b. R=[] or H =[] if r or e does not exist
%                   2. varargin{2} contains the adjacency matrix object of the fixed modules:
%                           a. If not specified, then the defaul is that
%                           only the entries in 'R' are fixed
%                           b. If specified, defined similarly as
%                           varargin{1} with fields G,R,H
%                               b1: since R is a binary matrix
%                                   in the reference paper, 'R' in varargin{2} should be
%                                   identical to 'R' in varargin{1}, i.e.,actuators
%                                   dynamics are always known
%                               b2: H in varargin{1} and varargin{2} must
%                                   have the same dimension
%         moduleindex: 1x 2 row vector that contains the indices of target module, where 
%                      a. row vector [i j] denotes a directed edge from w_i to w_j and thus the module G_{ji}  
%                      b. Input and the output of the target module must be
%                      measured
%         maxiter: Maximum iteration allowed for searching for subsets;
%                  If not specified, the default is maxiter=infinity
% Outputs: Test: 1 if the target module is identifiable,
%                0 if the sufficient condition is not satisfied after evaluating all possible subsets,
%                NaN if sufficient condition not yet sastisfied when the maximum iteration is reached
%          
% Reference:
% (a) S. Shi, X. Cheng and P.M.J. Van den Hof, "Single module identifiability in linear dynamic networks with partial excitation and measurement",  
%     arXiv preprint arXiv:2012.11414,2020.
%
%   Author:  Shengling Shi
%            Control Systems Group
%            Eindhoven University of Technology.
%   Version: 1.0
%   Date:    23- Sep-2021
%
%   Note: Rely on the functions "mindismaxpath", "adj2adjm", and "PowerSet"

%% Check inputs
if ~((nargin == 2)||(nargin == 3) )
  error('number of input arguments is not correct')
elseif nargin == 2
    maxiter=inf;
end  
 if ~(size(moduleindex,2) == 2)
  error('the dimension of the input moduleindex is not correct')
end  
adjobj = varargin{1};   
testmatrix = adjobj.G;    
for n=1:size(moduleindex,1)
    testindex = moduleindex(n,:);
    if ~(testmatrix(testindex(2),testindex(1))==1)
          error('the specified target modules do not exist in the network')  
    end
end

L = size(adjobj.G,1); % number of nodes
if isfield(adjobj,'C')&&(~isempty(adjobj.C))
    if size(adjobj.C,2)~=L
         error('Incorrect column dimension for C matrix')  
    end
else
    error('Measurement matrix C is not specified in varargin{1}')
end
[~,initalmeasure] = find(adjobj.C); % The measured internal signals C

if ~(all(ismember(moduleindex,initalmeasure) )) % Check whether input and output are measured
    error('Input or output of the target module is not measured, so the sufficient identifiability condition is trivially not satisfied')
end

%% Basic build up
% Build the extended graph of the network

[A,~] = adj2adjm(adjobj);
K = size(adjobj.R,2); % number of external excitations
p = size(adjobj.H,2); % number of white noises

if length(varargin)<2
    adjobjfix.G=zeros(L,L);
    if p==0
        adjobjfix.H=[];
    else
        adjobjfix.H=zeros(L,p);
    end
    if K==0
        adjobjfix.R=[];
    else
        adjobjfix.R=zeros(L,K);
    end
else
    adjobjfix = varargin{2};
    if ~(all(size(adjobjfix.R)==size(adjobj.R)) && all(size(adjobjfix.H)==size(adjobj.H)))  
          error('The specified R and H field in varargin{2} should have the same dimension as the R and H field in varargin{1}, respectively')  
    end  
end


%% verify identifiability of a single module
Gmax = adjobj.G;
Gmaxfix = adjobjfix.G;
Rmax = adjobj.R;
Rmaxfix = adjobjfix.R;
Hmax = adjobj.H;
Hmaxfix = adjobjfix.H;
RHmatrix = [Hmax Rmax]+[Hmaxfix Rmaxfix];

% Set up basic sets
    wj=moduleindex(2);
    inneighbor = find(Gmax(wj,:));
    inneighborfix = find(Gmaxfix(wj,:));
    Nj = setdiff(inneighbor,intersect(initalmeasure,inneighborfix )); % Definition of Nj set in Proposition 2, excluding measured signals that have known edge to wj
    Wjbar = setdiff(moduleindex(1),inneighborfix ); % check whether the target module is already known or not
    
% Start verification  
    if ~isempty(Wjbar)   % If the target module is not already known      
        if K+p==0          
           Xj=[];       
        else  
            X=1:1:K+p;
            Xex = find(RHmatrix(wj,:)==1);
            Xj = setdiff( X,Xex);      
            Xj = Xj+L;      % Indices of intially present external signals, which have no unknown edge to wj, in extended graph 
        end     % Indices of external signals, which have no unknown edge to wj, in extended graph       
        Cfull = setdiff(initalmeasure,moduleindex(1));
        
        PPX=PowerSet(Xj);
        PPC=PowerSet(Cfull);
        % Search for a subset of measured signals
        iter=0;
        for nx=length(PPX):-1:1         % Search subsets
            for nc = length(PPC):-1:1
                iter=iter+1;
                Cbar = PPC{nc};    
                Xbar = PPX{nx};
                UnionSet=union(Nj,Cbar);
                [b1,~,~]= mindismaxpath(A',Xbar ,UnionSet);
                [b2,~,~]= mindismaxpath(A',Xbar ,setdiff(UnionSet,moduleindex(1))); 
                [b3,~,~]= mindismaxpath(A',Xbar ,Cbar); 
                if (b1==b2+1) && (b2 == b3)     % Identifiability conditions in Proposition 2
                    Test=1;
                    return
                elseif maxiter< iter
                    Test= nan;
                    return        
                else              % check which module in the MISO is not identifiable
                    Test=0;
                end
            end
        end
        
    else
        Test=1;  
    end

end

