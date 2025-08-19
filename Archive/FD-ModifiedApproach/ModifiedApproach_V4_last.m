clear variables
close all; 
clc
%%%%%%%%%%%%%%%%%%%%%%%%
% Stability P 69/355 PDF
% α_m = 1 for all orders (also mentioned in 100/355 
% [6] says it should decrease, 
%           but also =1 over all for 2m-2m scheme
%%%%%%%%%%%%%%%%%%%%%%%%
c = 1; L = 1; TM = .6;
Lambda = L/10;
h = 0.01;
N = round(2*L/h+1);
X = -L:h:L;
dT = 0.001:0.0005:0.04;%

K = -2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
K(1,N) = 1; 
K(N,1) = 1;
K = sparse((1/h^2)*K);

U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),1);
alfa = zeros(length(dT),1);

N_ord = 10;              % 2== modified equation
U_all_dt = zeros(length(dT),N,N_ord);

for m=1:N_ord
    for kk=1:length(dT)
        alfa(kk) = c*dT(kk)/h;
        T = 0:dT(kk):TM;
        Un = U0';
        Unm = U0';
        KK = sparse(N,N);
        for oo=1:m
            KK = KK + ((c*dT(kk))^(2*oo)/factorial(2*oo))*K^(oo);
        end
        % figure()    
        for tt=2:length(T)
            Unp = (2*speye(N,N)+2*KK)*Un-Unm;
            Unm = Un;
            Un = Unp;
            % plot(X,Un), ylim([-.5 1.5]), grid on
            % drawnow
        end
        U_all_dt(kk,:,m) = Unp;
        X1 = -L:(h/100):L;
        Unp1 = interp1(X,Unp,X1);
        theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        err_M(kk,m) = error_val;
        clear T error_val
    end
end

% figure(11)
% loglog(dT,(err_M),"LineWidth",1),
% xlabel("\Delta t"),ylabel("Log10(error)");

max_err = zeros(1,N_ord);
for jj=1:N_ord
    M_err = find(err_M(:,jj)<0.001);
    max_err(jj) = dT(M_err(end));
end

ordo = 1:N_ord;
figure(1)
err_M1 = err_M;
err_M1(err_M1>0.1) = 0.1;
surf(2*ordo,dT,err_M1),shading flat;zscale log;
xlabel("order"); ylabel("\Delta t"),zlabel("error");
colorbar

figure(2)
plot(2*ordo,c*max_err/h,"LineWidth",1.4),hold on,grid on;
plot(2*ordo,c*max_err/h,"*","LineWidth",1.4)
xlabel("order"), ylabel("CFL number \alpha",'Interpreter','tex')

figure(3),loglog(dT,err_M(:,1)),hold on,loglog(dT,err_M(:,2))
% hold on,loglog(dT,err_M(:,3)),hold on,loglog(dT,err_M(:,4))

% figure(11)
% for ik=1:N_ord
%     U1 = U_all_dt(:,:,ik);
%     U1(abs(U1)>2)=0;
%     surf(X,dT,U1),shading flat,zlim([-2 2]);
%     drawnow;
%     pause(0.1)
% end

%%
beta_p = @(p) (2^(2*p-1))*((factorial(p-1))^2)/factorial(2*p);
BB = 0;
for ii=1:N_ord
    BB = BB + beta_p(ii);
end
BB = BB^(-1/2);
disp(BB)
