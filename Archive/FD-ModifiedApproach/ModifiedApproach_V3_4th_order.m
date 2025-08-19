clear variables
% close all; 
clc
%                       V3 : book eq 4.39 + 4.24 for Δ (PDF-55) and 4.41 for Δ²(PDF-59)
c = 1; L = 1; TM = .6;
Lambda = L/10;
h = 0.01;
N = round(2*L/h+1);
X = -L:h:L;
dT = 0.001:0.0001:0.025;
% dT = 0.001;

K11 = -2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
K11(1,N) = 1; 
K11 = sparse((4/3/h^2)*K11);
K12 = -2*diag(ones(N,1))+diag(ones(N-2,1),2)+diag(ones(N-2,1),-2);
K12(1,N-1) = 1; 
K12 = sparse((-1/3/4/h^2)*K12);
K1 = K11+K12;
K2 = 6*diag(ones(N,1))-4*diag(ones(N-1,1),1)-4*diag(ones(N-1,1),-1)...
    +diag(ones(N-2,1),2)+diag(ones(N-2,1),-2);
K2(1,N) = -4; 
K2(2,N) =1;
K2(1,N-1) = 1; 
K2 = sparse((1/h^4)*K2);


U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),1);
alfa = zeros(length(dT),1);
U_all_dt = zeros(length(dT),N);

for kk=1:length(dT)
    alfa(kk) = c*dT(kk)/h;
    T = 0:dT(kk):TM;
    Un = U0';%(0.5*exp(-((X-c*dT(kk))/Lambda).^2) + 0.5*exp(-((X+c*dT(kk))/Lambda).^2))';
    Unm = U0';
    KK = ((dT(kk)*c)^2)*K1+(((dT(kk)*c)^4)/12)*K1*K1;
    for tt=1:length(T)
        Unp = (2*speye(N,N)+KK)*Un-Unm;
        Unm = Un;
        Un = Unp;
            % if mod(tt,5)==0
            %     plot(X,Unp)
            %     ylim([-1.2 1.2]);
            %     pause(0.05);
            %     drawnow
            % end
    end
    U_all_dt(kk,:) = Unp;
    X1 = -L:(h/10):L;
    Unp1 = interp1(X,Unp,X1);
    theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
    error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
    err_M(kk) = error_val;
    % clear T error_val
end
% figure(10);
% plot(X,theo_fn,X,Unp),legend("theor","simul");


%%
figure(1),% subplot(2,1,1)
loglog(dT,(err_M),"LineWidth",1.4),
xlabel("\Delta t"),ylabel("Log10(error)");

% figure(1),subplot(2,1,2)
% plot(dT,alfa,"LineWidth",1.4)
% xlabel("\Delta t"),ylabel("Linear CFL");

%%
% alpha_infp = @(p) 2*((-1)^(p-1));

% alpha_pm = @(m,p) ((-1)^(p-1))*(2*factorial(m)^2)/(factorial(m-p))/(factorial(m+p));
% 
% function AA = a0m(m)
% alpha_pm = @(m,p) ((-1)^(p-1))*(2*factorial(m)^2)/(factorial(m-p))/(factorial(m+p));
% AA = 0;
% for p=1:m
%     AA = AA + alpha_pm(m,p)/(p^2);
% end
% AA = -2*AA;
% end
% 
% function Dl = m_Delta_h(m,U,j,h)
% Dl = a0m(m)*U(j);
% for q=1:m
%     Dl= Dl + (alpha_pm(m,q)/(q^2))*(U(j+q)+U(j-q));
% end
% Dl = Dl/(h^2);
% end

