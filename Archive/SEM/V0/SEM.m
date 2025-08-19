close all
clear
clc

Mesh.Order = 4;
Mesh.Nel = [15 10];
% dim = [0 10; 0 20];

dim = [0 15; 0 10];

% BB1=[2 3; 0 4];
% BB2=[2 3; 10 20];

BB1=[6 7; 0 3];
BB2=[6 7; 6 10];

Mesh = Mesher(Mesh,dim);
% remove elements
[~,ElList] = BB(Mesh,BB1);
Mesh.Connectivity(ElList,:)= [];
[~,ElList] = BB(Mesh,BB2);
Mesh.Connectivity(ElList,:)= [];

% change coords here

idd = BB(Mesh,[6 7; 4 5]);
Mesh.Nodes(idd,2) = [4.25;4.25;4.75;4.75];

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
[M,K] = makeMK(Mesh,kappa,rho);

% save("MK_Full_Room.mat","M","K","Mesh");

% %%
% for el = 1:size(Mesh.Connectivity,1)
%         nodeIDs = Mesh.Connectivity(el,2:5);
%         x_coords = Mesh.Nodes(nodeIDs,1);
%         y_coords = Mesh.Nodes(nodeIDs,2);
%         if (min(x_coords)==2) && (min(y_coords)==5)
%             source_ID = nodeIDs;
%             source_x = x_coords;
%             source_y = y_coords;
%         end
% end
% 
% Smin_x = min(source_x); Smin_y = min(source_y);
% Smax_x = max(source_x); Smax_y = max(source_y);
% 
% Spt_ind = ((Mesh.Nodes(:,1)>=Smin_x)&(Mesh.Nodes(:,1)<=Smax_x)&(Mesh.Nodes(:,2)>=Smin_y)&(Mesh.Nodes(:,2)<=Smax_y));
% 
% figure(21),
% plot(Mesh.Nodes(Spt_ind,1),Mesh.Nodes(Spt_ind,2),"*r");

