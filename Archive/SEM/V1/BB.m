function [NList,ElList] = BB(obj,BoundBox)
%% MeshClass.BB
% Found the bounding box and return all the elements inside of
% the BB
%
% Syntax:
%   obj = obj.BB( def );
%
% Inputs:
%   BoundBox : BB with [xmin xmax; ymin ymax;zmin zmax]
%
% Outputs:
%   ElList : Logical list of the elements touched by the BB
%   NList : Logical list of the nodes inside the BB
%
% See also,
NList = logical((obj.Nodes(:,1)>=BoundBox(1,1)).*(obj.Nodes(:,1)<=BoundBox(1,2)) ...
    .*(obj.Nodes(:,2)>=BoundBox(2,1)).*(obj.Nodes(:,2)<=BoundBox(2,2)) );
nodeName= obj.RefNodes(NList,1);
clear aux
aux = obj.Connectivity(:,2:end);
aux = reshape(aux',1, 4*size(obj.Connectivity,1))';
El = zeros(numel(aux),1);
for inode = 1:numel(nodeName)
    El = (aux == nodeName(inode)) + El;
end
ElList = El>0;
clear El
ElList = reshape(ElList,[4 size(obj.Connectivity,1)]);
ElList = sum(ElList,1)>0;
ElList = ElList.';
end