function obj = MakeSEMMesh(obj)
obj.Dim = 2;
%obj.Order = order;
order = obj.Order
dx = 1e50*ones(1,obj.Dim);

[x,~,~] = lglnodes(order);
[y,~,~] = lglnodes(order);

out = zeros(size(obj.Connectivity,1),numel(x),numel(y),2);
for iel = 1:size(obj.Connectivity,1)
    elcoords = obj.Nodes(obj.Connectivity(iel,2:end),:);
    dx(1) = min([dx(1); abs(elcoords(1,1) - elcoords(2,1)); abs(elcoords(3,1) - elcoords(4,1))]);
    dx(2) = min([dx(2); abs(elcoords(1,2) - elcoords(4,2)); abs(elcoords(3,2) - elcoords(2,2))]);
    for ix = 1:numel(x)
        for iy = 1:numel(y)
            out(iel,ix,iy,:) = LinearMap([x(ix) y(iy)],elcoords);
        end
    end
end

dx = dx;
dx=[0.1 0.1];
bb = [min(obj.Nodes); max(obj.Nodes)];
L = abs(bb(1,:) - bb(2,:));
expo = [max(ceil(log10(L))) min(floor(log10(dx)))];

expFactor = 6;
nodes_new = floor(10^(-expo(2)+expFactor)*out);
out = floor(10^(-expo(2)+expFactor)*out);

aux = zeros(size(out,1)*(order+1)^2, 2);
ino = 1;
for iel = 1:size(out,1)
    for iy = 1 : size(out,3)
        for ix = 1 : size(out,2)
            aux(ino,:) = squeeze(nodes_new(iel,ix,iy,:));
            ino = ino + 1;
        end
    end
end
nodes_new = unique(aux,'rows');
ElementList = zeros(size(out,1),(order+1)^2);
org = zeros(4,1);
for iel = 1:size(out,1)
    ino = 1;
    for iy = 1 : size(out,3)
        for ix = 1 : size(out,2)
            ind = find(sum(squeeze(out(iel,ix,iy,:))'==nodes_new,2)==obj.Dim);
            ElementList(iel,ino) = ind;

            if ix == 1 && iy == 1
                org(1) = ino;
            elseif ix == size(out,2) && iy == 1
                org(2) = ino;
            elseif ix == size(out,2) && iy == size(out,3)
                org(3) = ino;
            elseif ix == 1 && iy == size(out,3)
                org(4) = ino;
            end

            ino = ino + 1;
        end
    end
end
P = 10^(expo(2)-expFactor)*nodes_new;
nel = linspace(1,size(obj.Connectivity,1),size(obj.Connectivity,1))';
obj.Nodes = [];
obj.Connectivity = [];
obj.Nodes = P;
ind = linspace(1,(order+1)^2,(order+1)^2)';
ind(org) = [];
ind = cat(1,org,ind);
obj.Connectivity =cat(2,nel,ElementList(:,ind));
obj.ElOrder = ind;

end
function MapCoord = LinearMap(ref,coords)

N(:,1) = (1-ref(:,1)).*(1-ref(:,2))*0.25;
N(:,2) = (1+ref(:,1)).*(1-ref(:,2))*0.25;
N(:,3) = (1+ref(:,1)).*(1+ref(:,2))*0.25;
N(:,4) = (1-ref(:,1)).*(1+ref(:,2))*0.25;

MapCoord = N*coords;
end

function [x,w,P]=lglnodes(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, 'Spectral Methods
%   in Fluid Dynamics,' Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation + 1
N1=N+1;
% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';
% The Legendre Vandermonde Matrix
P=zeros(N1,N1);
% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.
xold=2;
while max(abs(x-xold))>eps
    xold=x;
    P(:,1)=1;
    P(:,2)=x;
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end
w=flip(2./(N*N1*P(:,N1).^2));
x = flip(x');
w = w';
end