function filename = create_sphere_rips_demo( num_points, dimension )
% Example: construct the dense distance matrix for sphere
filename = ['sphere_' num2str( dimension ) '_' num2str( num_points ) '.DDM'];
RandStream.setGlobalStream(RandStream('mt19937ar','seed',pi));

%% create actual data
points = zeros( dimension, num_points );
cur_num_points = 0;
while cur_num_points < num_points
    random_point = 2 * (rand( dimension, 1 ) - 0.5 );
    if norm(random_point) > 0.75 && norm(random_point) < 1
        cur_num_points = cur_num_points + 1;
        points(:, cur_num_points) = random_point;
    end
end 

%% compute distance matrix
distance_matrix = zeros( num_points );
for i = 1:num_points
    for j = 1:num_points
        distance_matrix(i,j) = norm( points(:, i) - points(:, j));
    end
end
% scatter(points(1,:), points(2,:));

%% save to disk in DIPHA format
Save_Dense_Distance_Matrix( points, distance_matrix, filename );
end