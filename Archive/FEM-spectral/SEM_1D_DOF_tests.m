clear variables
% close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks out for r = 1:5 & m = 1;
% CFL page 197/355 of PDF >>> case of normal SEM + leap-frog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for r =2 ; 3 & m = 2 >> CFL page 199/355 
% max CFL sim = α_m/2 for both r=2 & r=3 :/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 1; L = 1; TM = 0.6;
Lambda = L/10;
r = 2;                                  % order of interp r-1
m = 1;                                  % order of taylor approx

alpha_L = [0.005:0.005:0.42 0.6 0.8 1];
% hList = 0.0009:0.00005:0.003;
dTList =  0.001;
hList = c*dTList./alpha_L;

[x_gll,w_gll] = lglnodes(r);        % positions and weights in [0;1]
n_gll = length(x_gll);              % length of the corresponding vector

err_M = zeros(length(dTList),length(hList));
DOF_L = zeros(length(hList),1);
N0_L = zeros(length(hList),1);
parfor jx=1:length(hList)
    for jt=1:length(dTList)
        h = hList(jx);
        dT = dTList(jt);
        disp(c*dT/h)
        T = 0:dT:TM;
        N0=round(2*L/h+1);
        X = linspace(-L,L,(N0-1)*(n_gll-2)+N0);     % position (interp pts included)
        DOF_L(jx) = length(X); N0_L(jx) = N0;
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
                            KK = sparse(zeros(size(K)));
                            NN = inv(M)*K; %#ok<*MINV>
                            for oo=2:m
                                KK = KK + (((-1)^oo)*(c*dT)^(2*oo)/factorial(2*oo))*((NN)^(oo));
                            end
        U0 = exp(-(X/Lambda).^2);               % Initial wave
        Unm = U0';
        Un = U0';
        Unp = 0;
        for kt=2:length(T)
            Unp = 2*Un-Unm-(((c*dT)^2)*(NN)+2*KK)*Un;
            Unm = Un;
            Un = Unp;
        end
        % final theoretical solution
        TF = (round(TM/2))*2-TM;
        theo_fn = 0.5*exp(-((X-c*TF)/Lambda).^2) + 0.5*exp(-((X+c*TF)/Lambda).^2);

        error_val = trapz(X,(theo_fn'-Unp).^2)/trapz(X,theo_fn.^2);
        err_M(:,jx) = error_val;
    end
end
%%
err_M(err_M>10^2) = 1;


% DOF_L_r5=DOF_L;
% N0_L_r5=N0_L;
% err_M_r5=err_M;
% save("data_r5.mat","DOF_L_r5","N0_L_r5","err_M_r5")

figure()
loglog(N0_L,err_M,"LineWidth",1),
xlabel("\Delta t"),ylabel("Error");


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