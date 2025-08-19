clear variables
% close all; 
clc

c = 1; L = 1; TM = 1;
Lambda = L/10;
h = 0.01;
N_ord = 16;

% dT = [0.00001,0.00002,0.000025,0.00004,0.00005,0.0001,0.0002,0.00025,0.0004,0.0005,...
%     0.001,0.002,0.0025,0.004,0.005,0.01,0.02,0.025,0.04,0.05,0.1,0.2,0.25,0.4,0.5];
dT = [0.00001,0.00002,0.00004,0.00005,0.00008,0.0001,0.00016,0.0002,0.00025,0.00032,...
    0.0004,0.0005,0.0008,0.001,0.00125,0.0016,0.002,0.0025,0.004,0.005,0.00625,...
    0.008,0.01,0.0125,0.02,0.025,0.03125,0.04,0.05,0.0625,0.1];


N = round(2*L/h+1);
X = -L:h:L;

K = 0.5*diag(ones(N-1,1),1)-0.5*diag(ones(N-1,1),-1);
K(1,1:2) = [-1,1];
K(N,N-1:N) = [-1,1];
K = sparse(K);

U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),N_ord);
ord = 1:16;

for oo=1:N_ord
    for kk=1:length(dT)
        alfa = c*dT(kk)/h;
        T = 0:dT(kk):TM;
        Un = U0';       
        UL1 = sparse(ord(1)*(-alfa*K));
        UL2 = sparse(ord(2)*(1/factorial(2))*(-alfa*K)^2);
        UL3 = sparse(ord(3)*(1/factorial(3))*(-alfa*K)^3); 
        UL4 = sparse(ord(4)*(1/factorial(4))*(-alfa*K)^4);
        UL5 = sparse(ord(5)*(1/factorial(5))*(-alfa*K)^5);
        UL6 = sparse(ord(6)*(1/factorial(6))*(-alfa*K)^6);
        UL7 = sparse(ord(7)*(1/factorial(7))*(-alfa*K)^7);
        UL8 = sparse(ord(8)*(1/factorial(8))*(-alfa*K)^8);
        UL9 = sparse(ord(9)*(1/factorial(9))*(-alfa*K)^9);
        UL10 = sparse(ord(10)*(1/factorial(10))*(-alfa*K)^10);
        UL11 = sparse(ord(11)*(1/factorial(11))*(-alfa*K)^11);
        UL12 = sparse(ord(12)*(1/factorial(12))*(-alfa*K)^12);
        UL13 = sparse(ord(13)*(1/factorial(13))*(-alfa*K)^13);
        UL14 = sparse(ord(14)*(1/factorial(14))*(-alfa*K)^14);
        UL15 = sparse(ord(15)*(1/factorial(15))*(-alfa*K)^15);
        UL16 = sparse(ord(16)*(1/factorial(16))*(-alfa*K)^16);
        for tt=1:length(T)
            Unn = (speye(N,N)+UL1+UL2+UL3+UL4+UL5+UL6+UL7+UL8+UL9...
                +UL10+UL11+UL12+UL13+UL14+UL15+UL16)*Un;
            Un = Unn;
        end
        theo_fn = exp(-((X-c*TM)/Lambda).^2);
        error_val = trapz(X,(theo_fn'-Un).^2)/trapz(theo_fn'.^2);
        err_M(kk,oo) = error_val;
        clear alfa  T error_val
    end
end

err_M1 = log10(err_M);
max_err = zeros(1,8);
for jj=1:N_ord
    M_err = find(err_M1(:,jj)<0);
    max_err(jj) = dT(M_err(end));
end
%%
ordo = 1:N_ord;
figure(1)
surf(ordo,dT,err_M1),
xlabel("order"); ylabel("\Delta t"),zlabel("Log10(error)");
colorbar

figure(2)
plot(ordo,max_err,"LineWidth",1.4),hold on,grid on;
plot(ordo,max_err,"*","LineWidth",1.4)
xlabel("order"),ylabel("\Deltat")
title("max stable \Deltat")

