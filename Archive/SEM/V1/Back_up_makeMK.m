function [M,K,mtemp,ktemp] = makeMK(Mesh,kappa,rho)

[x,w] = lglnodes(Mesh.Order);
%calc the function derivatives
[N,dN] = ShapeFunctionLagrange(x);

nodes = zeros(1,numel(Mesh.Connectivity(1,2:end)));
%K = sparse(zeros(2*size(Mesh.Nodes,1),2*size(Mesh.Nodes,1)));
%M = sparse(zeros(2*size(Mesh.Nodes,1),2*size(Mesh.Nodes,1)));

K = sparse(zeros(size(Mesh.Nodes,1),size(Mesh.Nodes,1)));
M = sparse(zeros(size(Mesh.Nodes,1),size(Mesh.Nodes,1)));

% build a (NxN)x2xM matrix for phi
dphi = zeros(size(N,2),size(N,2),size(N,1)*size(N,1),2);%
for iwx = 1 : size(N,2)
    for iwy = 1 : size(N,2)
        count = 1;
        for in2 = 1 : size(N,1)
            for in1 = 1 : size(N,1)
                phi(iwx,iwy,count) = N(in1,iwx)*N(in2,iwy);
                dphi(iwx,iwy,count,1) = dN(in1,iwx)*N(in2,iwy); %(ponto de integracao epsilon, ponto d integracao xsi, funcoes de forma, derivada em eps ou xsi) derivada em relacao a epsilon
                dphi(iwx,iwy,count,2) = N(in1,iwx)*dN(in2,iwy); % derivada em relacao a xsi.
                count = count + 1 ;
            end
        end
    end
end

for iel = 1:size(Mesh.Connectivity,1)
    clear ktemp mtemp
    %ktemp = zeros(2*(Mesh.Order+1)^2,2*(Mesh.Order+1)^2);
    %mtemp = zeros(2*(Mesh.Order+1)^2,2*(Mesh.Order+1)^2);
    ktemp = zeros((Mesh.Order+1)^2, (Mesh.Order+1)^2);
    mtemp = zeros((Mesh.Order+1)^2, (Mesh.Order+1)^2);
    nodes(Mesh.ElOrder) = Mesh.Connectivity(iel, 2:end);
    coords = Mesh.Nodes(nodes,:);
    
    % Material properties 
    matC = kappa;  
    matRho = rho;  
    
    for iwy = 1:numel(x)
        for iwx = 1:numel(x)
            % Calculate Jacobian
            J = squeeze([dphi(iwx,iwy,:,1) dphi(iwx,iwy,:,2)]) * coords; 
            if det(J) <= 0
                error('Negative Jacobian, please check the mesh');
            end
            
            gradN = (J' \ squeeze([dphi(iwx,iwy,:,1) dphi(iwx,iwy,:,2)]));
            
            % Build B 
            B = zeros(2, numel(nodes));
            %build B
            %B = zeros(3,2*numel(nodes));
            %for i = 1 : numel(nodes)
            %    B(:,(i*2-1):2*i) = [gradN(i,1) 0;0 gradN(i,2);gradN(i,2) gradN(i,1)]; % aqui tirar os termos de cizalahamnto. assim a matriz B fica 2x2.
            %end

            % NN = zeros(2,2*numel(nodes));
            % c = 1;
            % for i = 2 : 2 : 2*numel(nodes)
            %     NN(1,i-1)   = phi(iwx,iwy,c);
            %     NN(2,i)     = phi(iwx,iwy,c);
            %     c = c+1;
            % end
            % 
            % matC = [kappa 0 0; 0 kappa 0; 0 0 kappa];
            % matRho = [rho 0; 0 rho];
            dNx = gradN(1,:);
            dNy = gradN(2,:);
            for i = 1:numel(nodes)
                B(:,i) = gradN(:,i)';
             end
            ktemp = ktemp + w(iwx)*w(iwy)*det(J)*(B'*matC*B);
            NN = squeeze(phi(iwx,iwy,:))';              
            mtemp = mtemp + w(iwx)*w(iwy)*det(J)*(NN'*matRho*NN);
        end
    end
    
    % Assembly into global matrices
    K = AssemblyGlobal(K, ktemp, nodes);  % Stiffness matrix
    M = AssemblyGlobal(M, mtemp, nodes);  % Mass matrix 
end
end

function Global = AssemblyGlobal(Global,Local,nodes)
dim = 1;
Pos = zeros(1,dim*numel(nodes));
Pos = nodes;
% for i = 1 : dim
%     Pos(i:dim:end) = (nodes-1)*2 +i;
% end

Global(Pos,Pos) = Global(Pos,Pos) + Local;
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

function [phi,dphi] = ShapeFunctionLagrange(x)
suport = x;
phi   = LagrangePolynomial(x,suport)';
dphi  = DifferentialLagrangePolynomial(x,suport)';
end

function [ LagrangePolynomialEvaluations ] = DifferentialLagrangePolynomial( EvaluationPoints,InterpolationPoints )
% This function computes the the derivative of the
% Lagrange polynomial basis functions
% based on the interpolation points defined in the COLUMN
% vector 'InterpolationPoints'. The evaluation points defined in
% 'EvaluationPoints' can be either row- or column formated.
[nrow, ncol] = size(EvaluationPoints);
if (ncol > nrow)
    EvaluationPoints = EvaluationPoints'; % get column vector
end
n = numel(EvaluationPoints);
m = numel(InterpolationPoints);
LagrangePolynomialEvaluations = zeros(n,m);
PolynomialConstants = ones(m,1);
% Normalizing Constants
for k = 1:m
    for i = 1:m
        if i ~= k
            PolynomialConstants(k) = ...
                PolynomialConstants(k)/(InterpolationPoints(k) - InterpolationPoints(i));
        end
    end
end
% Polynomial Derivative Evaluation
for k = 1:m
    for j = 1:m
        temp = ones(n,1);

        for i = 1:m
            if i ~= j && i ~= k
                temp = temp.*(EvaluationPoints - InterpolationPoints(i));
            end
        end

        if j ~= k
            LagrangePolynomialEvaluations(:,k) = ...
                LagrangePolynomialEvaluations(:,k) + ...
                temp*PolynomialConstants(k);
        end
    end
end
end

function [ LagrangePolynomialEvaluations ] = LagrangePolynomial( EvaluationPoints,InterpolationPoints )
%The function computes a matrix of Lagrange polynomials with respect to the nodal basis 'InterpolationPoints'.
%The polynomials are evaluated at the points 'EvaluationPoints'.
%The polynomials are stacked columnwise, with the rows corresponding to the evaluation points.
%The input is assumed to be one dimensional arrays.
[nrow, ncol] = size(EvaluationPoints);
if (ncol > nrow)
    EvaluationPoints = EvaluationPoints'; % get column vector
end
n = numel(EvaluationPoints);
m = numel(InterpolationPoints);
LagrangePolynomialEvaluations = ones(n,m);
PolynomialConstants = ones(m,1);
% Normalizing Constants
for k = 1:m
    for i = 1:m
        if i ~= k
            PolynomialConstants(k) = ...
                PolynomialConstants(k)/(InterpolationPoints(k) - InterpolationPoints(i));
        end
    end
end
% Polynomial Evaluation
for k = 1:m
    for i = 1:m
        if i ~= k
            LagrangePolynomialEvaluations(:,k) = ...
                LagrangePolynomialEvaluations(:,k).*(EvaluationPoints - InterpolationPoints(i));
        end
    end
    LagrangePolynomialEvaluations(:,k)=LagrangePolynomialEvaluations(:,k)*PolynomialConstants(k);
end
end