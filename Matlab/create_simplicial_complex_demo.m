function create_simplicial_complex_demo()
points = [];
delta = [0, 0.5, 1, 1.5, 2, 2.5];
for y=1:6
    for x=delta(y)+1:6-delta(y)
        points = [points; x, y;];
    end
end
points = points';

pointsValue = [1.1, 2, 2.3, 1.3, 3.3, 5.2, 2.8, 5, 10, 3, 1.1, 2.2,...
    10, 10, 7.1, 2.7, 8, 9.1, 2, 2, 3];

edges = [1 2; 2 3; 3 4; 4 5; 5 6; 7 8; 8 9; 9 10; 10 11; 12 13; 13 14; 14 15; 16 17; 17 18; 19 20;
    1 7; 7 12; 12 16; 16 19; 19 21; 2 8; 8 13; 13 17; 17 20; 3 9; 9 14; 14 18; 4 10; 10 15; 5 11;
    2 7; 3 8; 8 12; 4 9; 9 13; 13 16; 5 10; 10 14; 14 17; 17 19; 6 11; 11 15; 15 18; 18 20; 20 21;];
edges = edges' - 1; % index starting from 0

faces = [1 2 7; 2 7 8; 2 3 8; 3 8 9; 3 4 9; 4 9 10; 4 5 10; 5 10 11; 5 6 11;
    7 8 12; 8 12 13; 8 9 13; 9 13 14; 9 10 14; 10 14 15; 10 11 15;
    12 13 16; 13 16 17; 13 14 17; 14 17 18; 14 15 18; 16 17 19; 17 19 20; 17 18 20; 19 20 21];
faces = faces' - 1; % index starting from 0

pointsX = points(1,:);
pointsY = points(2,:);

organizations = {};
organizations{1} = edges;
organizations{2} = faces;

Save_General_Simplicial_Complex(points, pointsValue, organizations, 'sim_complex_demo.dat');

%% plot the results ...
edges = edges + 1;
faces = faces + 1;

figure;
hold on;
for i=1:size(edges, 2)
    plot([pointsX(edges(1, i)), pointsX(edges(2, i))], [pointsY(edges(1, i)), pointsY(edges(2, i))]);
end
title('Dimension 1');
hold off;

figure;
hold on;
for i=1:size(faces, 2)
    patch([pointsX(faces(1, i)), pointsX(faces(2, i)), pointsX(faces(3, i))],...
        [pointsY(faces(1, i)), pointsY(faces(2, i)), pointsY(faces(3, i))], 'green');
end
title('Dimension 2');
hold off;
end

