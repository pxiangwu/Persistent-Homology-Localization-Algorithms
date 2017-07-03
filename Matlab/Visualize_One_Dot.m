function Visualize_One_Dot( phi, pers_list, red_list, bd_list, d, i )

cpts = pers_list{d}(:,i);
red_verts = red_list{d}{i};
bd_verts = bd_list{d}{i};

assert( (length(size(phi)) == 2) | (length(size(phi)) == 3) );
if( length(size(phi)) == 2 )
    figure, imshow(phi,[]);
    hold on; plot(red_verts(2,:), red_verts(1,:), 'g.');
    hold on; 
    plot(bd_verts(2,:), bd_verts(1,:), 'cs'); 
    plot(cpts(2:2:end), cpts(1:2:end), 'r*'); 
    axis([1, size(phi,2), 1, size(phi,1)]); grid on; hold off;
else
    figure, plot3(red_verts(1,:), red_verts(2,:), red_verts(3,:), 'g.');
    hold on; 
    plot3(bd_verts(1,:), bd_verts(2,:), bd_verts(3,:), 'bs'); 
    plot3(cpts(1:3:end), cpts(2:3:end), cpts(3:3:end), 'r*'); 
    axis([1, size(phi,1), 1, size(phi,2), 1, size(phi,3)]); grid on; hold off;
end

end