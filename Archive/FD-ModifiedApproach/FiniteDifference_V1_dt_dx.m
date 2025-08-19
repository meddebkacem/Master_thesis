clear variables
close all; 
clc

c = 1; L = 1; TM = 0.6;
Lambda = L/10;
N_ord = 4;

hL = 0.0005:0.0005:0.05;
% dTL = 0.0005:0.0005:0.1;

dTL = [0.001:0.0002:0.12 0.15]; %0.1
nt = unique(round(TM./dTL));
dTL = flip(TM./nt);


uu = ((mod(2,hL)*2000));
uu(uu~=0) = true;
hL = hL(uu==0);

err_M = zeros(length(dTL),length(hL));
alfa = zeros(length(dTL),length(hL));

for jx=1:length(hL)
    h = hL(jx);
    parfor jt=1:length(dTL)
        dT = dTL(jt);
        N = round(2*L/h+1);
        X = -L:h:L;
        K = diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);
        K(1,N) = -1;
        K(N,1) = 1;
        K = sparse(0.5*K);
        alfa(jt,jx) = c*dT/h;
        T = 0:dT:TM;
        U0 = exp(-(X/Lambda).^2);
        Un = U0';
        UL = sparse(N,N);
        for lk=1:N_ord
            UL = UL + sparse((1/factorial(lk))*(-alfa(jt,jx)*K)^lk);
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
        err_M(jt,jx) = error_val;
        % clear alfa  T error_val UL U0 Un Unn
    end
end

%%

err_M1 = err_M;
err_M1(err_M>1) = 1;
err_M1(isnan(err_M1)) = 1;


ordo = 1:N_ord;
figure(1)
surf(hL,dTL,err_M1),shading flat, zscale("log"), xscale("log"), yscale("log");
xlabel("h"); ylabel("\Delta t"),zlabel("(error)");
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
