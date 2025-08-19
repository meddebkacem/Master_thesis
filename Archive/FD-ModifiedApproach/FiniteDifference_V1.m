clear variables
close all; 
clc
                                tic
c = 1; L = 1; TM = 0.6;
Lambda = L/10;
h = 0.01;
N_ord = 20;%12;

% dT = [0.001,0.00125,0.0016,0.002,0.0025,0.004,0.005,0.00625,...
%     0.008,0.01,0.0125,0.02,0.025,0.03125,0.04,0.05,0.0625,0.1];
dT = [0.001:0.0002:0.12 0.15]; %0.1
nt = unique(round(TM./dT));
dT = flip(TM./nt);

N = round(2*L/h+1);
X = -L:h:L;

K = diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);
K(1,N) = -1;
K(N,1) = 1;

% K(1,1:2) = [-1,1];
% K(N,N-1:N) = [-1,1];
K = sparse(0.5*K);

U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),N_ord);
% eigen_UL = zeros(length(dT),N_ord);

for oo=1:N_ord
    parfor kk=1:length(dT)
        alfa = c*dT(kk)/h;
        T = 0:dT(kk):TM;
        Un = U0';
        UL = sparse(N,N);
        for lk=1:oo
            UL = UL + sparse((1/factorial(lk))*(-alfa*K)^lk);
        end
        % figure()
        for tt=1:length(T)
            Unn = (speye(N,N)+UL)*Un;
            Un = Unn;
            % plot(X,Unn), ylim([-.5 1.5]), grid on
            % drawnow
        end
        theo_fn = exp(-((X-c*TM)/Lambda).^2);
        error_val = trapz(X,(theo_fn'-Un).^2)/trapz(theo_fn'.^2);
        err_M(kk,oo) = error_val;
        % eigen_UL(kk,oo) = max(abs(eigs(UL)));
    end
end
                                    toc
%%
err_M(err_M>1)= 1;
ordo = 1:N_ord;
figure(1)
surf(ordo,dT,err_M),shading flat; zscale("log")
xlabel("order"); ylabel("\Delta t"),zlabel("error");
colorbar
clim([10^(-4) 10^(0)]);
set(gca,'ColorScale','log')
%%
max_err = zeros(1,8);
for jj=1:N_ord
    M_err = find(err_M(:,jj)<0.008);
    max_err(jj) = dT(M_err(end));
end
figure(2)
plot(ordo,c*max_err/h,"LineWidth",1.4),hold on,grid on;
plot(ordo,c*max_err/h,"*","LineWidth",1.4)
xlabel("order"), ylabel("CFL number \alpha",'Interpreter','tex')

