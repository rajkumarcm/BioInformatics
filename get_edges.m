function edges = get_edges(adj_mat,nodes)

edges = 0;
for i = 1:length(nodes)
    for j = 1:length(nodes)
        if i ~= j && j >= i
            edges = edges + adj_mat(nodes(i),nodes(j));
        end
    end
end

end