clear variables
close all
clc

c = 1; L = 1; TM = 0.6;
Lambda = L/10;

r = 5;                              % order of interp r-1
[x_gll,w_gll] = lglnodes(r);        % positions and weights in [0;1]
n_gll = length(x_gll);              % length of the corresponding vector
                               
ml = 1:1:20;                                 % order of taylor approx

h = 0.01;                           % for Δt:cte
alp = 0.01:0.005:1.1;
dTList = alp*h;

N0=round(2*L/h+1);
X = linspace(-L,L,(N0-1)*(n_gll-2)+N0);     % position (interp pts included)
M = zeros(1,length(X));                     % Mass matrix diagonal
for ij=1:(n_gll-1):length(X)-n_gll+1
    X(ij:ij+n_gll-1)=x_gll*(X(ij+n_gll-1)-X(ij))+X(ij);     % positions for plotting purpose
    M(ij:ij+n_gll-1)=M(ij:ij+n_gll-1) + (X(ij+n_gll-1)-X(ij))*w_gll';
end
M = sparse(diag(M));    % making sparse mass matrix
D = zeros(r+1,r+1);     % elementary derivatives with gauss lobatto interp
for i=1:r+1
    for j=1:r+1
        D(i,j) =alternative_dl(i,x_gll,x_gll(j));
    end
end
Ke = zeros(r+1);        % elementary stiffness matrix
for ii=1:r+1
    for jj=1:r+1
        for kk=1:r+1
            Ke(ii,jj)=Ke(ii,jj)+w_gll(kk)*D(ii,kk)*D(jj,kk);
        end
    end
end
K = zeros(length(X));   % global stiffness matrix (matching the Mass matrix)
for il=1:(n_gll-1):length(X)-n_gll+1
    K(il:il+n_gll-1,il:il+n_gll-1) = K(il:il+n_gll-1,il:il+n_gll-1) + Ke;
end
K = K./h;

err_M = zeros(length(dTList),length(ml));

for jt=1:length(dTList)
    dT = dTList(jt);
    err_par = zeros(1,length(ml));      % ADDED FOR //
    parfor jm=1:length(ml)
        m = ml(jm);
        T = 0:dT:TM;
        KK = sparse(zeros(size(K)));
        NN = inv(M)*K; %#ok<*MINV>
        for oo=2:m
            KK = KK + (((-1)^oo)*(c*dT)^(2*oo)/factorial(2*oo))*((NN)^(oo));
        end
        U0 = exp(-(X/Lambda).^2);               % Initial wave
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        for kt=2:length(T)
            Unp = 2*Un-Unm-(((c*dT)^2)*(NN)+2*KK)*Un;
            Unm = Un;
            Un = Unp;
        end
        % final theoretical solution + normalized error
        theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
        
        error_val = trapz(X,(theo_fn'-Unp).^2)/trapz(X,theo_fn.^2);
        err_par(jm) = error_val;
    end
    clear Unp Unm Un U0; 
    err_M(jt,:) = err_par;
end
%%
err_M(err_M>10^2) = 1;

max_t = zeros(1,length(ml));
for jk=1:length(ml)
    M_err = find(err_M(:,jk)<0.001);
    max_t(jk) = dTList(M_err(end));
end

          
xlabel("\Delta t"),ylabel("Error");
% loglog(alp,err_M,"LineWidth",1),
% xlabel("CFL (\alpha)"),ylabel("Error");figure()
loglog(dTlist,err_M,"LineWidth",1),

% figure()
% loglog(h,err_M,"LineWidth",1),
% xlabel("\Delta x"),ylabel("Error");
% xlim([8*10^(-4) 3*10^(-3)])
% ylim([10^(-6) 10^(1)])

figure()
plot(ml,max_t/h)
xlabel("order"),ylabel("CFL")



% % % function: interpolation derivatives for gll
function y = alternative_dl(j,x,z)
y = 0;
n = length(x);
for l=1:n
    if not(l==j)
        k = 1/(x(j)-x(l));
        for m=1:n
            if not(m==j) && not(m==l)
                k = k*(z-x(m))/(x(j)-x(m));
            end
        end
        y = y + k;
    end
end
end