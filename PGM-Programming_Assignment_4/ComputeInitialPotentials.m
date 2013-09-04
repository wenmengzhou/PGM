%COMPUTEINITIALPOTENTIALS Sets up the cliques in the clique tree that is
%passed in as a parameter.
%
%   P = COMPUTEINITIALPOTENTIALS(C) Takes the clique tree skeleton C which is a
%   struct with three fields:
%   - nodes: cell array representing the cliques in the tree.
%   - edges: represents the adjacency matrix of the tree.
%   - factorList: represents the list of factors that were used to build
%   the tree. 
%   
%   It returns the standard form of a clique tree P that we will use through 
%   the rest of the assigment. P is struct with two fields:
%   - cliqueList: represents an array of cliques with appropriate factors 
%   from factorList assigned to each clique. Where the .val of each clique
%   is initialized to the initial potential of that clique.
%   - edges: represents the adjacency matrix of the tree. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function P = ComputeInitialPotentials(C)

% number of cliques
N = length(C.nodes);

% initialize cluster potentials 
P.cliqueList = repmat(struct('var', [], 'card', [], 'val', []), N, 1);
P.edges = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% First, compute an assignment of factors from factorList to cliques. 
% Then use that assignment to initialize the cliques in cliqueList to 
% their initial potentials. 

% C.nodes is a list of cliques.
% So in your code, you should start with: P.cliqueList(i).var = C.nodes{i};
% Print out C to get a better understanding of its structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = length(C.factorList);
assigned = zeros(1, M);

for i = 1:N
    P.cliqueList(i).var = C.nodes{i};
    P.cliqueList(i).card = zeros(1, length(P.cliqueList(i).var));
    for j = 1:M
            [tr idx] = ismember(C.factorList(j).var, P.cliqueList(i).var);
            if all(tr)
                P.cliqueList(i).card(idx) = C.factorList(j).card;  %尽可能的获得变量的card
                 if ~assigned(j)
                    assigned(j) = i;
                 end
            end
    end
    
    for k = 1:length(P.cliqueList(i).var)               %factor并没有填满所有的C_i， 导致有的C_i中变量的card为0， 我们需要把这些变量card获得，否则用ones语句就无法正确初始化变量概率为1(line 70).
        if P.cliqueList(i).card(k) == 0
            for j = 1:M
                [tr idx] = ismember( P.cliqueList(i).var(k), C.factorList(j).var);
                if all(tr)
                    P.cliqueList(i).card(k) = C.factorList(j).card(idx);
                    break;
                end
            end
        end
    end
    
    P.cliqueList(i).val = ones(1, prod(P.cliqueList(i).card));
end

for i = 1:M
    idx = assigned(i);
    P.cliqueList(idx) = FactorProduct(P.cliqueList(idx),  C.factorList(i));
end


P.cliqueList = StandardizeFactors(P.cliqueList); 

P.edges = C.edges;
end

