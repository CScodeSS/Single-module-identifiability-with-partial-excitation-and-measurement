clear all
clc
% Example of allocating actuators and sensors to achieve generic identifiability

%% Specifiy the network
% Specify network topology
adj.G = [0 0 0 0; %1
         1 0 0 0; %2
         0 1 0 0; %3
         1 1 1 0;]; %4
adj.H = [1;0;0;0];
adj.R = [];
adj.C = [1 0 0 0;
    0 0 0 1];
varargin{1}=adj;

% Specify the known modules in the network
adjfix.G = [0 0 0 0; %1
         0 0 0 0; %2
         0 0 0 0; %3
         0 1 1 0;]; %4
adjfix.H = [0;0;0;0];     % fix.H and R must have the same dimenion as H and R in varargin{1}, respectively
adjfix.R = [];
varargin{2}=adjfix;

% Plot the network 
 ExtendGraph = adj2adjm(adj);
 Graph=digraph((ExtendGraph)',{'w1','w2','w3','w4','e'});
 H=plot(Graph,'Layout','force');
 % Highlight the measured internal nodes in green
[~,sensor]=find(adj.C);
highlight(H,sensor,'NodeColor','g')
 
 % Veify whether the target module is already identifiable
 targetmodule=[1 4];
 [Test] = modidfanalysis_partial(varargin,targetmodule);
if Test ==1         
 highlight(H,targetmodule,'EdgeColor','g')   % Green if identifiable
else
     highlight(H,targetmodule,'EdgeColor','r') % Red if sufficient condition not satisfied
end
 

%%  Conduct actuator and sensor allocation
[vararginOut] = modidfsynthesis_partial(varargin,targetmodule);

% Plot new network
adjnew=vararginOut{1};

% Since function 'adj2adjm' only works without field 'C', we build a new
% network without the C matrix in 'adjnew' 
adjDya.G=adjnew.G;
adjDya.H=adjnew.H;
adjDya.R=adjnew.R;

% Draw the new network
figure
GraphNew=digraph((adj2adjm(adjDya))',{'w1','w2','w3','w4','e','r'});
H2=plot(GraphNew,'Layout','force');

% Highlight the measured internal nodes in green
C = adjnew.C;
[~,sensor]=find(C);
highlight(H2,sensor,'NodeColor','g')

% Verify identifiability in the new network
[Test] = modidfanalysis_partial(vararginOut,targetmodule);
if Test ==1         
 highlight(H2,targetmodule,'EdgeColor','g')   % Green if identifiable
else
     highlight(H2,targetmodule,'EdgeColor','r') % Red if sufficient condition not satisfied
end

