function [ clusters ] = cluster_v2( X )
% Name: Rajkumar Conjeevaram Mohan
% Coursework: Introduction to Bioinformatics
% Clustering genes into three clusters using Euclidean distance

% Initializing vars....
clc;
[m,~] = size(X);

clusters = cell(1,3);
% Initialising medoids
medoids = [1,2,3];

% Initialize similarity matrix
d_mat = get_sim_matrix(X);
% display('Similarity matrix');
% display(d_mat);

% Step 1: Initialize datapoints closest to medoids
sim_ = zeros(size(medoids));
for i = 1:m
    if ~sum(medoids==i)
        for m_i = 1:length(medoids)
            sim_(m_i) = d_mat(medoids(m_i),i);
        end
        [~,I] = min(sim_);
        indx = randi([1 length(I)]);
        % In case if there is more than one same minimal value
        I = I(indx);
        dt_pts = clusters{I};
        dt_pts(end+1) = i;
        clusters{I} = dt_pts;
    end
    
end

% Iteration beings
max_epochs = 3;
% Step 2
for epoch = 1:max_epochs
    
    % Step 2.1 
    % Recompute the medoids in each cluster
    for m_i = 1:length(medoids)
        temp = clusters{m_i};
        % Include also the medoid of this cluster in 
        % its dataset while re-estimating medoid
        temp(end+1) = medoids(m_i);
        new_medoid = estimate_medoid(temp);
       
        % if existing medoid is the same as estimated medoid
        % then make no change else swap them.
        if medoids(m_i) ~= new_medoid
            clusters{m_i} = temp(temp~=new_medoid);
            medoids(m_i) = new_medoid;
        end
    end
    
    % Step 2.2 Compute distance between medoid of one cluster and 
    % non-medoids of other clusters
    total_dtpts = 1:m;
    non_medoids = setdiff(total_dtpts,medoids);
    
    for nm_i = 1:length(non_medoids)
        cost = zeros(1,3);
        for m_i = 1:length(medoids)
            cost(m_i) = d_mat(m_i,non_medoids(nm_i));
        end
        
        [~,ob_medoid] = min(cost);
        ex_medoid = get_medoid(clusters,non_medoids(nm_i));
        
        % if non-medoid dt_pt's existing cluster is the same
        % as obtained, then make no change, else make a swap
        if ob_medoid ~= ex_medoid
            
            % Now assign the non_medoid dt_pt to new cluster
            % ob_medoid
            temp = clusters{ob_medoid};
            temp(end+1) = non_medoids(nm_i);
            clusters{ob_medoid} = temp;
            
            % Then delete the non_medoid dt_pt from old cluster
            % to avoid redundancy
            temp = clusters{ex_medoid};
            clusters{ex_medoid} = temp(temp~=non_medoids(nm_i));
            
        end
    end
end

% Merge medoids with its corresponding cluster bins
for m_i = 1:length(medoids)
    temp = clusters{m_i};
    temp(end+1) = medoids(m_i);
    clusters{m_i} = temp;
end


end

function medoid = get_medoid(clusters,dt_pt)
    for i = 1:length(clusters)
        temp = clusters{i};
        if sum(temp==dt_pt)==1
            medoid = i;
            break;
        end
    end
end

function medoid = estimate_medoid(cluster)
    d_mat = get_sim_matrix(cluster);
    d_mat = sum(d_mat,2);
    [~,I] = min(d_mat);
    % In case if there is more than one same minimal value
    indices = randi([1 length(I)]);
    medoid = cluster(I(indices));
end

function d_mat = get_sim_matrix(X)
    m = size(X,1);
    d_mat = zeros(m,m);
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
end
