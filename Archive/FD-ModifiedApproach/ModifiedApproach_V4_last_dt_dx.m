clear variables
% close all; 
clc

c = 1; L = 1; TM = .6;
Lambda = L/10;

% hL = 0.0005;

dT = 0.0001:0.00001:0.1;
uu = ((mod(TM,dT)*10000));
uu(uu~=0) = true;
dT = dT(uu==0);

% hL = divisors(200)/1000;
hL = 0.0001:0.0001:0.2;
uu = ((mod(2,hL)*2000));
uu(uu~=0) = true;
hL = hL(uu==0);

% dT = 0.001;
% hL = dT./(0.01:0.01:2);
% hL = hL(100:2:174);

m = 2;

err_M = zeros(length(dT),length(hL));
alfa = zeros(length(dT),length(hL));

for jj=1:length(hL)
    h = hL(jj);
    N = round(2*L/h+1);
    X = -L:h:L;
   
    K = -2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
    K(1,N) = 1; 
    K(N,1) = 1;
    K = sparse((1/h^2)*K);

    U0 = exp(-(X/Lambda).^2);

    for kk=1:length(dT)
        alfa(kk,jj) = c*dT(kk)/h;
        T = 0:dT(kk):TM;
        Un = U0';
        Unm = U0';
        KK = sparse(N,N);
        for oo=1:m
            KK = KK + ((c*dT(kk))^(2*oo)/factorial(2*oo))*K^(oo);
        end
    
        for tt=2:length(T)
            Unp = (2*speye(N,N)+2*KK)*Un-Unm;
            Unm = Un;
            Un = Unp;
                % if mod(tt,10)==0
                %     plot(X,Unp)
                %     ylim([-0.2 1.2]);
                %     pause(0.05);
                %     drawnow
                % end
        end
        X1 = -L:(h/100):L;
        Unp1 = interp1(X,Unp,X1);
        theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        err_M(kk,jj) = error_val;
        % clear T error_val
    end
end

%%
err_M(err_M>2)=1;
figure(11)
loglog(hL,err_M(3,:),"LineWidth",1),
xlabel("\Delta x"),ylabel("(error)");
figure(12)
loglog(dT,err_M(:,5),"LineWidth",1),
xlabel("\Delta t"),ylabel("(error)");


figure()
surf(hL,dT,err_M),shading flat,xscale log,yscale log,zscale log;
xlabel("\Delta x"),ylabel("\Delta t"),zlabel("(error)");
err_M(isnan(err_M))=1;
figure()
surf(hL,dT,err_M),shading flat,xscale log,yscale log,zscale log;
xlabel("\Delta x"),ylabel("\Delta t"),zlabel("(error)");


% ordo = 0.01:0.01:2;
% ordo= ordo(100:2:174);
% figure(2)
% plot(ordo,c*dT./hL,"LineWidth",1.4),hold on,grid on;
% % plot(ordo,c*dT./hL,"*","LineWidth",1.4)
% xlabel("order"), ylabel("CFL number \alpha",'Interpreter','tex')

% figure(12)
% dT1 = dT(1):(dT(1)/100):dT(end);
% err_int = interp1(dT,err_M,dT1);
% loglog(dT1,(err_int),"LineWidth",1),
% xlabel("\Delta t"),ylabel("Log10(error)");


%%

beta_p = @(p) (2^(2*p-1))*((factorial(p-1))^2)/factorial(2*p);
BB = 0;
for ii=1:m
    BB = BB + beta_p(ii);
end
BB = BB^(-1/2);
cond_t = BB*h/c;
disp(append("max Δt: ",num2str(cond_t)))
cond_h = dT/BB;
disp(append("max Δx: ",num2str(cond_h)))