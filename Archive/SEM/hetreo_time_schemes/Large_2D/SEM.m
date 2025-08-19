% close all
clear
clc

Mesh.Order = 4;

% [nn,ee] = readGmshQuadMesh("Rooms_V5.msh");
[nn,ee] = readGmshQuadMesh("Rooms_test_III.msh");        % Rooms_V7_2507.msh
            % % FOR INP FILES (FROM LUCIO)
            % [nn, ee, eType] = readinp("Job-2.inp");
            % plotMesh(nn, ee, eType)

            % For Rooms_V5
            % el_inverted = [51:199 317:331 356:382];
            % for kk=1:length(el_inverted)
            %     temp1 = ee(el_inverted(kk),:);
            %     ee(el_inverted(kk),2) = temp1(4);
            %     ee(el_inverted(kk),4) = temp1(2);
            %     clear temp1;
            % end

            % [nn,ee] = readGmshQuadMesh("Rooms_quad_V2.msh");
Mesh.Connectivity = [];
Mesh.Nodes = [];
                % % FOR INP FILES (FROM LUCIO)
                % Mesh.Connectivity(:,2:5) = ee(:,2:5);
Mesh.Connectivity(:,2:5) = ee(:,1:4);
Mesh.Connectivity(:,1) = 1:1:size(ee,1);
                % % FOR INP FILES (FROM LUCIO)
                % Mesh.Nodes = nn(:,2:3);
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

        plot(x_coords, y_coords, 'b-');hold on;
            % text(mean(x_coords),mean(y_coords),num2str(el),"FontSize",14)
    end

end

Mesh = MakeSEMMesh(Mesh);

% Plot nodes
plot(Mesh.Nodes(:,1), Mesh.Nodes(:,2), 'ro', 'MarkerFaceColor', 'r');


kappa=1;
rho = 1;
[M,K,M_el,K_el] = makeMK(Mesh,kappa,rho);

% save("MK_Full_Room1.mat");
%% testing for connectivity issues
AA_M = zeros(size(Mesh.Nodes,1));
AA_K = zeros(size(Mesh.Nodes,1));
for ll=1:el
    coooord = Mesh.Connectivity(ll,2:end);
    AA_M(coooord,coooord) = AA_M(coooord,coooord)+M_el(:,:,ll);
    AA_K(coooord,coooord) = AA_K(coooord,coooord)+K_el(:,:,ll);
    clear coooord
end
