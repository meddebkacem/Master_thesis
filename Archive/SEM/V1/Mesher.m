function obj = Mesher(obj,Dimensions)

%nodes
nnodes = obj.Nel + 1;
[X, Y] = ndgrid(0:(nnodes(1)-1), 0:(nnodes(2)-1));  
indices = (1:prod(double(nnodes)))';                
obj.RefNodes = [indices, X(:), Y(:)];              
index = 1;
obj.Connectivity = zeros(prod(double(obj.Nel)),5);
for j = 1 : obj.Nel(2)
    for i = 1 :obj.Nel(1)
        obj.Connectivity(index,:) = [i + (j-1)*obj.Nel(1)...
            (i+0) + (j+0-1)*nnodes(1)...
            (i+1) + (j+0-1)*nnodes(1)...
            (i+1) + (j+1-1)*nnodes(1)...
            (i+0) + (j+1-1)*nnodes(1)];
        index = index + 1 ;
    end
end
%coords
x = (Dimensions(1,2)-Dimensions(1,1))*double(obj.RefNodes(:,2))/double(max(obj.RefNodes(:,2)))+Dimensions(1,1);
y = (Dimensions(2,2)-Dimensions(2,1))*double(obj.RefNodes(:,3))/double(max(obj.RefNodes(:,3)))+Dimensions(2,1);

obj.Nodes = [x y];

end