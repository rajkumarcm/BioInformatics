clear;
clc;
close all;
% Adjency matrix...
adj_mat = zeros(9,9);
adj_mat(1,2:3) = 1;
adj_mat(2,[1,3]) = 1;
adj_mat(3,[1,2,4,5]) = 1;
adj_mat(4,[3,6]) = 1;
adj_mat(5,[3,6]) = 1;
adj_mat(6,[4,5,7,8,9]) = 1;
adj_mat(7,[6,8,9]) = 1;
adj_mat(8,[6,7,9]) = 1;
adj_mat(9,[6,7,8]) = 1;

degree_vec = sum(adj_mat,2);

% Degree distribution
deg = unique(degree_vec);
deg_vec = zeros(1,length(deg));
for i = 1:length(deg)
    n = sum(degree_vec==deg(i));
    deg_vec(i) = n/9;
end
 bar(deg_vec,'b');

% Avergage clustering coefficient and its spectrum
cl_c = zeros(1,9);
for i = 1:9
    
    nodes = find(adj_mat(i,:) ~= 0);
    n = length(nodes);
    max_edges = factorial(n)/(factorial(2)*factorial(n-2));
    m = get_edges(adj_mat,nodes);
    % Since undirected, m = 2xm
    cl_c(i) = (m)/max_edges;
end
avg_clcoeff = mean(cl_c);
figure

% plot cl_c
cl_dis = zeros(1,length(deg));
for i = 1:length(deg)
    node_indices = find(degree_vec==deg(i));
    n = length(nodes);
    cl_dis(i) = sum(cl_c(node_indices))/n;
end
plot(deg,cl_dis,'o- b','LineWidth',2,'MarkerSize',10)

dis_mat = zeros(9,9);
dist_mat(1,2) = 1;
dist_mat(1,3) = 1;
dist_mat(1,4) = 2;
dist_mat(1,5) = 2;
dist_mat(1,6) = 3;
dist_mat(1,7) = 4;
dist_mat(1,8) = 4;
dist_mat(1,9) = 4;

dist_mat(2,3) = 1;
dist_mat(2,4) = 2;
dist_mat(2,5) = 2;
dist_mat(2,6) = 3;
dist_mat(2,7) = 4;
dist_mat(2,8) = 4;
dist_mat(2,9) = 4;

dist_mat(3,4) = 1;
dist_mat(3,5) = 1;
dist_mat(3,6) = 2;
dist_mat(3,7) = 3;
dist_mat(3,8) = 3;
dist_mat(3,9) = 3;

dist_mat(4,5) = 2;
dist_mat(4,6) = 1;
dist_mat(4,7) = 2;
dist_mat(4,8) = 2;
dist_mat(4,9) = 2;

dist_mat(5,6) = 1;
dist_mat(5,7) = 2;
dist_mat(5,8) = 2;
dist_mat(5,9) = 2;

dist_mat(6,7) = 1;
dist_mat(6,8) = 1;
dist_mat(6,9) = 1;

dist_mat(7,8) = 1;
dist_mat(7,9) = 1;

dist_mat(8,9) = 1;

% get unique distances
u_dist = unique(dist_mat);
u_dist = u_dist(2:end);
nodes_dist = zeros(1,length(u_dist));
for i = 1:length(u_dist)    
    nodes_dist(i) = length(find(dist_mat==u_dist(i)));
end

total_dist = sum(nodes_dist);
for i = 1:length(nodes_dist)
    nodes_dist(i) = nodes_dist(i)/total_dist;
end
figure
title('Shortest distance spectrum');
bar(u_dist,nodes_dist,'r')

% Betweenness centrality
centrality = zeros(1,9);
for i = 1:9
    for j = 1:9
        if i ~= j && j >= i
            
        end
    end
end