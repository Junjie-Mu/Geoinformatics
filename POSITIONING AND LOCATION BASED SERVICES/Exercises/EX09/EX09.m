%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EX09: Dijkstra
% POS&LBS A.A. 2023/2024
% MU JUNJIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
close all
clear all
node = {'London', 'Zurich', 'Milan', 'Paris', 'Wien', 'Istanbul', 'Kyoto'};
% arc columns = start node, end node, cost
arc = [1 2 3; ...
       2 1 3; ...
       1 4 10; ... 
       4 1 10; ...
       2 3 1; ...
       3 2 1; ... 
       3 4 13; ... 
       4 3 13; ... 
       3 5 8; ...
       5 3 8; ...
       4 5 3; ...
       5 4 3; ...
       5 6 6; ...
       6 5 6; ...
       6 4 2; ...
       4 6 2; ...
       6 7 4; ...
       7 6 4; ...
       2 7 18; ... 
       7 2 18; ];
 % Departure node 
first_id = find(strcmp(node,'Milan'));
 
 % Arrival node
final_id = find(strcmp(node,'Kyoto'));   
 
 % Costs: Initialize by setting all to infinity
ttn = inf(size(node));
% but the departure: cost 0
ttn(first_id) = 0;
 % prev_id: for each node the previous in the minimum path: set NaN
prev_id = nan(size(node));
 % but the departure: reachable from itself
prev_id(first_id) = first_id; 
 % visited_id: already explored nodes vector: set empty
visited_id = [];
 % id_to_visit: nodes that can be (costs already estimated) visited / explored: at the beginning only the departure
id_to_visit = [ first_id ];
while ~ismember(final_id,visited_id)
    % something exists to check: for example until the arrival has been visited / explored
    % presently visited /explored node: find minimum cost node in the set of nodes that can be visited
    mincost = min(ttn(id_to_visit));
    cuurent_id = find(ttn == mincost);
    % remove the presently visited node from the set of nodes that can be visited
    id_to_visit = id_to_visit(id_to_visit ~= cuurent_id);
    % ... and add the presently visited node to the set of already visited nodes
    visited_id = [visited_id, cuurent_id];
    % search outgoing arcs from visited node and check reachable nodes from currently visited node
    reachable_nodes = arc(arc(:, 1) == cuurent_id,:);
    % check: do not go back 
    reachable_nodes = reachable_nodes(~ismember(reachable_nodes(:,2),visited_id),:);    
    % Cycle For each reachable node
    for i = 1 : size(reachable_nodes,1)
            % check cumulated cost of the node; in case and existing cost exists compare new and existing; set the minimum
            % check and update the vector of the previous in the minimum paths
            item = reachable_nodes(i,:);
            currentEsstimate = item(1);
            index = item(2);
            cost = item(3) + mincost;
            if(cost < ttn(index)) 
                ttn(index) = cost;
                prev_id(index) = currentEsstimate;
            end
        	% if not yet in the set of nodes to visit, add it
            if(~ismember(index,id_to_visit))
                id_to_visit = [id_to_visit,index];
            end         
    end
end
index = final_id;
fprintf('best path: \n')
while(index ~= first_id)
    disp(node(index))
    index = prev_id(index);
end
disp(node(index))



