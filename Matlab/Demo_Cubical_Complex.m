%% Demo for cubical complex
close all;
clear;
clc;

%% plot the results
dataName = 'cubical_data/peaks_2d_demo.dat';
[phi, pers_list, red_list, bd_list] = Read_Pers_Results_Cubical(dataName);

plot([0, max(phi(:))], [0, max(phi(:))], 'g');
hold on;

% plot the persistence diagram for d dimension
d = 2;
for i=1:size(pers_list{d},2)
    cpts = pers_list{d}(:,i);
    birth = phi(cpts(1), cpts(2));
    death = phi(cpts(3), cpts(4));
    plot(birth, death, 'r.');
end
axis equal;
axis tight;
grid on;
title('Persistence Diagram');
hold off;

%% plot the cycles
for d=2
    for i=1:length(red_list{d})
        Visualize_One_Dot(phi, pers_list, red_list, bd_list, d, i);
    end
end