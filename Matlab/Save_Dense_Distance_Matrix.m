function Save_Dense_Distance_Matrix( points, distance_matrix, filename )
%% Open file for writing
fid = fopen( filename, 'w' );

%% File type identifier - Dense Distance Matrix 1
fwrite( fid, 1, 'uint' );

%% Total number of input points n
fwrite( fid, size( distance_matrix, 2 ), 'uint' );

%% Write dimension of each point
fwrite( fid, size( points, 1), 'uint' );

%% Write point positions
fwrite( fid, points(:), 'double' );

%% Floating point values for the coordinates of the points
fwrite( fid, distance_matrix(:), 'double' );

%% Close file
fclose(fid);
end