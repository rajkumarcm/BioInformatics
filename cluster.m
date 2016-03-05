function [ clusters ] = cluster( X )
% Name: Rajkumar Conjeevaram Mohan
% Tutorial: Introduction to Bioinformatics Tutorial2
% Clustering genes into two clusters using Manhattan
% distance and single-linkage method

% Initializing vars....
clc;
[m,~] = size(X);
d_mat = zeros(size(X));
clusters = cell(1,m);
for i = 1:m
    temp = zeros(1,m);
    temp(1) = i;
    clusters{i} = temp;
end

% Initialize similarity matrix
for d1 = 1:m
    for d2 = 1:m
        if d1 ~= d2
            % Manhattan distance
%             d_mat(d1,d2) = sum(abs(X(d1,:)-X(d2,:)));
            % Euclidean distance
              d_mat(d1,d2) = sqrt(sum((X(d1,:)-X(d2,:)).^2));
        end
    end
end
display('Similarity matrix');
display(d_mat);
%----------------------
formed_clusters = length(clusters);
while formed_clusters ~= 3
    % Iterate over clusters to group them...
    cluster_sim = Inf(formed_clusters);
    for c1 = 1:formed_clusters
        for c2 = 1:formed_clusters
            if c1 ~= c2
                % Samples in cluster c1
                cl1 = clusters{c1};
                dt_pts1 = cl1(cl1~=0);
                cl2 = clusters{c2};
                dt_pts2 = cl2(cl2~=0);
                
                if ~isempty(dt_pts1) && ~isempty(dt_pts2)
                    % distance matrix between dt_pts in two clusters
                    temp_dmat = Inf(length(dt_pts1),length(dt_pts2));
                    for i = 1:length(dt_pts1)
                        for j = 1:length(dt_pts2)
                            temp_dmat(i,j) = d_mat(dt_pts1(i),dt_pts2(j));
                        end
                    end

                    % get the most min value index
                    temp = min(min(temp_dmat));
                    index = randi([1,length(temp)]);
                    cluster_sim(c1,c2) = temp(index);                    
                end
            end
        end
    end
    indx = find(cluster_sim == min(min(cluster_sim)));
    indx = indx(randi([1,length(indx)]));
    cl2 = mod(indx,size(cluster_sim,2));
    if cl2 == 0
        cl2 = size(cluster_sim,2);
    end
    cl1 = ceil(indx/size(cluster_sim,1));
    
    % Merge clusters
    % get empty positions in cl1
    temp_cl1 = clusters{cl1};
    try
        temp_cl2 = clusters{cl2};
    catch
        display('catch me');
    end
    
    % m_cl -> merged cluster
    m_cl = horzcat(temp_cl1(temp_cl1~=0),temp_cl2(temp_cl2~=0));
    m_cl = horzcat(m_cl,zeros(1,m-length(m_cl)));
    clusters{cl1} = m_cl;
    clusters{cl2} = zeros(1,m);
    
    formed_clusters = sum(sum(vertcat(clusters{:}),2)~=0);
    
    new_clusters = cell(1,formed_clusters);
    j = 1;
    % Lets just get rid of clusters with no data points inside
    for i = 1:length(clusters)
        dt_pts = clusters{i};
        if ~isempty(dt_pts(dt_pts~=0))
            new_clusters{j} = dt_pts;
            j = j+1;
        end
    end
    clusters = new_clusters;
    clear new_clusters;
    
end


end

