close all
clear
clc

Mesh.Order = 4;
% Mesh.Nel = [10 15];
% dim = [0 10; 0 15];

% BB1=[2 3; 0 4];
% BB2=[2 3; 10 20];

% Mesh = Mesher(Mesh,dim);
% remove elements
%[~,ElList] = BB(Mesh,BB1);
%Mesh.Connectivity(ElList,:)= [];
%[~,ElList] = BB(Mesh,BB2);
%Mesh.Connectivity(ElList,:)= [];
[nn,ee] = readGmshQuadMesh("Rooms_quad.msh");
Mesh.Connectivity = [];
Mesh.Nodes = [];
Mesh.Connectivity(:,2:5) = ee;
Mesh.Connectivity(:,1) = 1:1:size(ee,1);
Mesh.Nodes = nn(:,1:2);

if 1
    figure;
    hold on;
    axis equal;

    % Plot elements
    for el = 1:size(Mesh.Connectivity,1)
        nodeIDs = Mesh.Connectivity(el,2:5);
        x_coords = Mesh.Nodes(nodeIDs,1);
        y_coords = Mesh.Nodes(nodeIDs,2);

        % Close the quadrilateral
        x_coords(5) = x_coords(1);
        y_coords(5) = y_coords(1);

        plot(x_coords, y_coords, 'b-');
    end

end

Mesh = MakeSEMMesh(Mesh);

% Plot nodes
plot(Mesh.Nodes(:,1), Mesh.Nodes(:,2), 'ro', 'MarkerFaceColor', 'r');


kappa=1;
rho = 1;
[M,K,M_el,K_el] = makeMK(Mesh,kappa,rho);

% save("MK_Full_Room1.mat");