function [points, distMatrix, pers_list, red_list, bd_list] = Read_Pers_Results_FullRips(dataName)
%% read point position data
dataFileName = dataName;
pers_fname = [dataFileName,'.pers'];
red_fname = [dataFileName,'.red'];
bd_fname = [dataFileName, '.bnd'];

data_fid = fopen(dataFileName, 'r');
[~] = fread(data_fid, 1, 'uint' ); % read file type identifier

numPoints = fread(data_fid, 1, 'uint'); % read the number of points
pointDim = fread(data_fid, 1, 'uint'); % read point dimension

points = fread(data_fid, numPoints*pointDim, 'double'); % read point positions
points = reshape(points, [pointDim, numPoints]);

distMatrix = fread(data_fid, numPoints*numPoints, 'double');
distMatrix = reshape(distMatrix, [numPoints, numPoints]);

fclose(data_fid);

%% read persistence data
pers_list = {};
pers_fid = fopen(pers_fname, 'r');
maxDim = fread(pers_fid, 1, 'uint'); % read the maximum dimension
persSize = fread(pers_fid, maxDim, 'uint');

for i=1:maxDim
    tmp = fread(pers_fid, persSize(i)*2*2, 'uint');
    tmp = reshape(tmp, [2*2, persSize(i)]);
    pers_list{i} = tmp;
end
fclose(pers_fid);

%% read reduction data
red_list = {};
for i=1:maxDim
    red_fname_withdim = [red_fname, '.', int2str(i)];
    red_fid = fopen(red_fname_withdim, 'r'); 
    [~] = fread(red_fid, 1, 'uint'); % trash
    redSize = fread(red_fid, 1, 'uint'); % read the number of reductions
    
    red_list{i} = {};
    for j=1:redSize
        listSize = fread(red_fid, 1, 'uint'); % read the size for this dimension
        tmp = fread(red_fid, listSize, 'uint');
        red_list{i}{j} = tmp;
    end
    fclose(red_fid);
end

%% read boundary data
bd_list = {};
for i=1:maxDim
    bd_fname_withdim = [bd_fname, '.', int2str(i)];
    bd_fid = fopen(bd_fname_withdim, 'r');
    [~] = fread(bd_fid, 1, 'uint'); % trash
    bdSize = fread(bd_fid, 1, 'uint'); % read the number of reductions
    
    bd_list{i} = {};
    for j=1:bdSize
        listSize = fread(bd_fid, 1, 'uint');
        tmp = fread(bd_fid, listSize, 'uint');
        bd_list{i}{j} = tmp;
    end
    fclose(bd_fid);
end

end