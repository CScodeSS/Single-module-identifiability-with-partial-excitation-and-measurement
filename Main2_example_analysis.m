clear all
clc
%% Example of testing generic identifiability of a single module
% Formulate network adjacency matrices
adj.G = [0 0 0 0; %1
         1 0 0 0; %2
         0 1 0 0; %3
         1 1 1 0;]; %4
adj.H = [1 0 0;0 1 0;0 0 0;0 0 1]; % Canonical noise model
adj.R = [];
adj.C =[1 0 0 0;
    0 1 0 0;
    0 0 0 1];
varargin{1}=adj;

% Specify the known transfer functions in the network
adj2.G=[0 0 0 0; %1
         0 0 0 0; %2
         0 0 0 0; %3
         0 1 0 0;]; % G42 is known
adj2.R=[];
adj2.H=zeros(4,3);
varargin{2}=adj2;

% Specify the target moule and then conduct analysis
targetmodule=[1 4]; % Edge from w1 to w4
[Test] = modidfanalysis_partial(varargin,targetmodule);

% Plot the network and show the result
ExtendGraph = adj2adjm(adj);
Graph=digraph((ExtendGraph)',{'w1','w2','w3','w4','e1','e2','e3'});
H=plot(Graph,'Layout','force');
if Test ==1         
 highlight(H,targetmodule,'EdgeColor','g')   % Green if identifiable
else
     highlight(H,targetmodule,'EdgeColor','r') % Red if sufficient condition not satisfied
end
% Highlight the measured internal nodes in green
C = adj.C;
[~,sensor]=find(C);
highlight(H,sensor,'NodeColor','g')
