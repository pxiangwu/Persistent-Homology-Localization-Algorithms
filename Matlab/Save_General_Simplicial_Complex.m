function Save_General_Simplicial_Complex( points, values, organizations, filename )
%% Open file for writing
fid = fopen( filename, 'w' );

%% File type identifier - General simplicial complex 2
fwrite( fid, 2, 'uint' );

%% Maximum dimension of simplex
fwrite( fid, size(organizations, 2), 'uint');

%% Total number of input points n
fwrite( fid, size( points, 2 ), 'uint' );

%% Write dimension of each point
fwrite( fid, size( points, 1), 'uint' );

%% Write point positions and values
fwrite( fid, points(:), 'double' );
fwrite( fid, values(:), 'double' );

%% Write indices of edges, faces, etc ...
num = size(organizations, 2);
for dim=1:num
    org = organizations{dim};
    N = size(org, 2);
    fwrite( fid, dim, 'uint' );
    fwrite( fid, N, 'uint' );
    fwrite( fid, org(:), 'uint' );
end

%% Close file
fclose(fid);
end