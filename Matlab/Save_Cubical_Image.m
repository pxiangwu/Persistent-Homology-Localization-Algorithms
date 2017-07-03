function Save_Cubical_Image( phi, dataname )
% Generate the data file compatible with the persistence homology program
fid = fopen( dataname, 'w');

% File type identifier - Cubical Image Data 0
fwrite(fid, 0, 'uint');

d = length(size(phi));
tmpM = [d; flipud(size(phi)')];
fwrite(fid, tmpM, 'uint');

D = double(phi(:));
fwrite( fid, D, 'double');
fclose(fid);

end

