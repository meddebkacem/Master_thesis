clear variables
close all; 
clc

c = 1; L = 1; TM = .6;
Lambda = L/10;
h = 0.01;

dT = 0.001:0.0001:0.025; %0.1

N = round(2*L/h+1);
X = -L:h:L;

K11 = -2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
K11(1,N) = 1; 
K11(N,1) = 1;
K11 = sparse((4/3/h^2)*K11);

K12 = -2*diag(ones(N,1))+diag(ones(N-2,1),2)+diag(ones(N-2,1),-2);
K12(1,N-1) = 1; 
K12(N,2) = 1;
K12 = sparse((-1/3/4/h^2)*K12);

K2 = 6*diag(ones(N,1))-4*diag(ones(N-1,1),1)-4*diag(ones(N-1,1),-1)...
    +diag(ones(N-2,1),2)+diag(ones(N-2,1),-2);
K2(1,N) = -4; 
K2(N,1) = -4;
K2(1,N-1) = 1; 
K2(N,2) = 1;
K2 = sparse((1/h^4)*K2);
                      

U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),1);
alfa = zeros(length(dT),1);
U_all_dt = zeros(length(dT),N);

for kk=1:length(dT)
    alfa(kk) = c*dT(kk)/h;
    T = 0:dT(kk):TM;
    Un = U0';
    Unn = (0.5*exp(-((X-c*dT(kk))/Lambda).^2) + 0.5*exp(-((X+c*dT(kk))/Lambda).^2))';
    K = sparse(((dT(kk)^2)*c^2)*(K11+K12)+(dT(kk)^4)*((c^4)/12)*K2);

    for tt=3:length(T)
        Unnn = (2*speye(N,N)+K)*Unn-Un;
        Un = Unn;
        Unn = Unnn;
        % if mod(tt,10)==0
        %     plot(X,Unnn)
        %     ylim([-0.2 1.2]);
        %     pause(0.02);
        %     drawnow
        % end
    end
    U_all_dt(kk,:) = Unnn;
    theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
    error_val = trapz(X,(theo_fn'-Un).^2)/trapz(theo_fn'.^2);
    err_M(kk) = error_val;
    clear T error_val
end

%%
figure(1),subplot(2,1,1)
plot(dT,log10(err_M),"LineWidth",1.4),
xlabel("\Delta t"),ylabel("Log10(error)");

figure(1),subplot(2,1,2)
plot(dT,alfa,"LineWidth",1.4)
xlabel("\Delta t"),ylabel("Linear CFL");

% NN = 1:N; 
% figure(),
% surf(NN,dT,U_all_dt),shading flat
% colorbar,
% ylim([0.01 0.015])
% 
% prob = dT(log10(err_M)>0);
% disp (0.6./prob);