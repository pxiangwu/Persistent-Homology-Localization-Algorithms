%% Demo for full rips-complex
close all;
clear;
clc;

%% plot the results
dataName = 'rips_data/sphere_2_50.DDM';
[points, distMatrix, pers_list, red_list, bd_list] = Read_Pers_Results_FullRips(dataName);

plot([0, max(distMatrix(:))], [0, max(distMatrix(:))], 'g');
hold on;

% plot the persistence diagram for d dimension
d = 2;
for i=1:size(pers_list{d},2)
    cpts = pers_list{d}(:,i);
    birth = distMatrix(cpts(1), cpts(2));
    death = distMatrix(cpts(3), cpts(4));
    plot(birth, death, 'r.');
end
axis equal;
axis tight;
grid on;
title('Persistence Diagram');
hold off;

%% plot the boundaries
for i=1:length(bd_list{d})
    figure;
    scatter(points(1,:), points(2,:)); hold on;
    ptsIdx = int32(bd_list{d}{i});
    scatter(points(1, ptsIdx), points(2, ptsIdx), [], 'r', 'filled');
end
title('Optimal cycle (indicated by red vertices)');
hold off;

