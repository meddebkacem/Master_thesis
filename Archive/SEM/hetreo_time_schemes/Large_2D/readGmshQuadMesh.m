function [nodes, quads] = readGmshQuadMesh(filename)
% READGMSHQUADMESH Read nodes and quadrilateral elements from GMSH .msh file
%   [nodes, quads] = readGmshQuadMesh(filename) reads a GMSH 2.2 format file
%   and returns:
%       nodes - Nx3 matrix of node coordinates [x,y,z]
%       quads - Mx4 matrix of quadrilateral element connectivity

    nodes = [];
    quads = [];
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Check file format
    line = fgetl(fid);
    if ~strcmp(line, '$MeshFormat')
        fclose(fid);
        error('Not a valid GMSH file format');
    end
    
    line = fgetl(fid);
    version = sscanf(line, '%f %d %d');
    if version(1) < 2.2
        fclose(fid);
        error('This function only supports GMSH 2.2 format');
    end
    
    % Skip to nodes section
    while ~feof(fid)
        line = fgetl(fid);
        if strcmp(line, '$Nodes')
            break;
        end
    end
    
    % Read nodes
    numNodes = str2double(fgetl(fid));
    nodes = zeros(numNodes, 3);
    
    for i = 1:numNodes
        line = fgetl(fid);
        data = sscanf(line, '%f %f %f %f');
        nodes(i,:) = data(2:4)'; % Skip node number
    end
    
    % Skip to elements section
    while ~feof(fid)
        line = fgetl(fid);
        if strcmp(line, '$Elements')
            break;
        end
    end
    
    % Read elements
    numElements = str2double(fgetl(fid));
    quadCount = 0;
    quads = zeros(1000, 4); % Preallocate (adjust size as needed)
    
    for i = 1:numElements
        line = fgetl(fid);
        data = sscanf(line, '%d %d %d %d %d %d %d %d %d %d %d');
        
        % Check if it's a quad element (type 3 in GMSH 2.2)
        if data(2) == 3
            quadCount = quadCount + 1;
            quads(quadCount,:) = data(end-3:end)'; % Last 4 numbers are node tags
        end
    end
    
    % Trim unused preallocated space
    quads = quads(1:quadCount,:);
    
    fclose(fid);
    
    % Display summary
    fprintf('Read %d nodes and %d quadrilateral elements\n', size(nodes,1), size(quads,1));
end