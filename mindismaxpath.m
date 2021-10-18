function [b,D,P] = mindismaxpath(adjacencymat,V1,V2)
% In a directed graph, compute a minimum disconnecting set and a set of maximum number of vertex
% disjoint paths from a vertex set V1 to set V2. This function relies on
% the "maxflow" function in Matlab
% 
% [b,D,P] = mindismaxpath(adjacencymat,V1,V2)
% Function created on Matlab 2021a
%
% Inputs: adjacencymat: Binary adjacency matrix of a directed graph ( NOTE: non-zero (i,j) entry denoes a directed edge from i to j )
%         V1: vector of the indices of nodes in V1      
%         V2: vector of the indices of nodes in V2
%
% Outputs: b: The maximum number of vertex disjoint paths from V1 to V2 ( = minimum size of V1-V2 disconnecting set)
%          D：Indices of the nodes in a minimum disconnecting set
%          P: The adjacency matrix of the subgraph that has the same vertices but only contains the edges in a set of maximum
%             number of vertex disjoint paths from V1 to V2
%
%   Author:  Shengling Shi
%            Control Systems Group
%            Eindhoven University of Technology.
%   Version: 1.1 
%   Date:    05- Aug-2021
%
%   Note: Algorithm follows from Chapter 9 of book "combinatorial optimization polyhedra and efficiency (volumn A)"



%% Check inputs
% Check the input arguments
    if nargin ~= 3
        error('incorrect number of input arguments');
    end
% Node indices in V1 and V2 should not exceed the dimension of adjacency
% matrix
L = length(adjacencymat);

    if isempty(V1) || isempty(V2)
        b=0;
        D = [];
        P = [];
        return
    end

    if (max(V1)>L) || (max(V2)>L)
        error('Node indices in V1 and V2 should not exceed the dimension of adjacency matrix');
    end

%% Reformulation of the graph by adding source and sink and then node splitting
% add source s that has edges to all nodes in V1
rows=zeros(1,L+1);
rows(V1)=1;
As = vertcat(horzcat(adjacencymat,zeros(L,1)), rows);
% add sink t to which all nodes in V2 have edges
columnt = zeros(L+1,1);
columnt(V2)=1;
A = vertcat(horzcat(As,columnt),zeros(1,L+2)); % Now s corrsponds to the seconed last column, and t corresponds the last one

% Node splitting starts: Any vertex v is splitted into v' and v''
 Asplit=[];
for i=1:(L+2)
    column = A(:,i);
    for k=1:2
       columnVs = [];
        for j = 1:(L+2)
            if k ==1
               if column(j)==0
                  columnVs = vertcat(columnVs,zeros(2,1)); 
                              
               else
                   columnVs = vertcat(columnVs,[0 ;1]);         
               end    
            else
                columnVs = zeros(2*(L+2),1);
                columnVs((i-1)*2+1,1)=1;
            end
            
        end
     Asplit = [Asplit columnVs];  % Adjacency matrix with splitted nodes
    end     
end
G= digraph(Asplit);
%% Compute the outputs
[b,GF,cs,~]= maxflow(G,(L+2)*2-2,(L+2)*2-1);       % b = the maximum flow for the graph with splitted nodes and weight 1 for the edges

% computer indices of nodes in disconnecting set in the original graph
set=[];
for i = 1: length(cs)
    set = [set ;successors(G,cs(i))];                         
end
indexvector = setdiff(set,cs);

D=[];
for i=1: length(indexvector)
    if rem(indexvector(i),2) ==0   % If it is even
       D=[D indexvector(i)/2]; 

    else
        D=[D (indexvector(i)+1)/2];     
    end
end

% Compute the output P
P1 = adjacency(GF);
P2 = P1(1:(2*L),1:(2*L));       % Remove the source s and the sink t
% Starts merging the splitted nodes
P =[];
for i=1:L
    column = P2(:,i*2-1);
    columnP=zeros(L,1);
    for j =1:L
       if column(j*2)==1
        
            columnP(j) = 1;        
       end
    end
    P = [P columnP];
end

end

