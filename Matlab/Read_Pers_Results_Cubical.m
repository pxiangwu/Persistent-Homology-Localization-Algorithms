function [ phi, pers_list, red_list, bd_list ] = Read_Pers_Results_Cubical( dataname ) 
%% read filter function data
dataFileName = dataname;
pers_fname = [dataname,'.pers'];
red_fname = [dataname,'.red'];
bd_fname = [dataname, '.bnd'];

fid = fopen( dataFileName, 'r' );
[~] = fread( fid, 1, 'uint' );

d = fread( fid, 1, 'uint' );
phi_dim = fread( fid, d, 'uint' );
phi_dim = flipud(phi_dim);
datasize = 1;
for i = 1:d
    datasize = datasize*phi_dim(i);
end
phi = fread( fid, datasize, 'double' );
phi = reshape( phi, phi_dim' );
fclose( fid );

%% read persistence data
pers_list = {};
pers_fid = fopen( pers_fname, 'r' );
tmpd = fread( pers_fid, 1, 'uint' );
assert( tmpd == d );
pers_dim = fread( pers_fid, d, 'uint' );
for i = 1:d 
    tmpM = fread( pers_fid, pers_dim(i)*2*d, 'uint' );
    tmpM = reshape( tmpM, [ 2*d, pers_dim(i) ] );
    pers_list{i} = tmpM;
end
fclose( fid );

%% read reduction data
red_list = {};
for i = 1:d 
    red_fname_withdim = [red_fname, '.', int2str(i)];
    red_fid = fopen( red_fname_withdim, 'r' );
    tmpd = fread( red_fid, 1, 'uint' );
    assert( tmpd == d );
    red_dim = fread( red_fid, d, 'uint' );

    red_list{i} = {};
    for j = 1:red_dim(1)
        tmpDataInfo = fread( red_fid, d, 'uint' );
        assert( tmpDataInfo(1) > 0 );
        assert( max(tmpDataInfo(2:end))==0 );
    
        tmpM = fread( red_fid, tmpDataInfo(1)*d, 'uint' );
        tmpM = reshape( tmpM, [ d, tmpDataInfo(1) ] );
    
        red_list{i}{j} = tmpM;
    end
    
    fclose( red_fid );
end

%% read boundary data
bd_list = {};
for i = 1:d 
    bd_fname_withdim = [bd_fname, '.', int2str(i)];
    bd_fid = fopen( bd_fname_withdim, 'r' );
    tmpd = fread( bd_fid, 1, 'uint' );
    assert( tmpd == d );
    bd_dim = fread( bd_fid, d, 'uint' );

    bd_list{i} = {};
    for j = 1:bd_dim(1)
        tmpDataInfo = fread( bd_fid, d, 'uint' );
        assert( tmpDataInfo(1) > 0 );
        assert( max(tmpDataInfo(2:end))==0 );
    
        tmpM = fread( bd_fid, tmpDataInfo(1)*d, 'uint' );
        tmpM = reshape( tmpM, [ d, tmpDataInfo(1) ] );
    
        bd_list{i}{j} = tmpM;
    end
    
    fclose( bd_fid );
end

end

