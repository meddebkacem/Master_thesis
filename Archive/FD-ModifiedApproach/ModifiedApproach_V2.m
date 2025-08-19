clear variables
close all; 
clc
%           V2 : doc2 eq 1.8 (PDF-3)
c = 1; L = 1; TM = .6;
Lambda = L/10;
h = 0.01;
N = round(2*L/h+1);
X = -L:h:L;

dT = 0.001:0.0001:0.025;
% dT = 0.001;

function Kn = Lp_m(N,h,m)
alpha_pm = @(m,p) ((-1)^(p-1))*(2*factorial(m)^2)/(factorial(m-p))/(factorial(m+p));
Kn = sparse(N,N);
for p=1:m
    Kph = (-2*diag(ones(N,1))+diag(ones(N-p,1),p)+diag(ones(N-p,1),-p));
    Kph(1,N-p+1) = 1; 
    Kph(N,p) = 1;
    Kph = sparse(Kph./((p*h)^2));
    Kn = Kn + alpha_pm(m,p)*Kph;
end
end

U0 = exp(-(X/Lambda).^2);
err_M = zeros(length(dT),1);
alfa = zeros(length(dT),1);
U_all_dt = zeros(length(dT),N);

for kk=1:length(dT)
    
    alfa(kk) = c*dT(kk)/h;
    T = 0:dT(kk):TM;
    Un = U0'; %(0.5*exp(-((X-c*dT(kk))/Lambda).^2) + 0.5*exp(-((X+c*dT(kk))/Lambda).^2))';
    Unm = U0';    
    m= 12;
    % KM = sparse(N,N);
    % for j=1:m
    %     Kj = Lp_m(N,h,j);       % 2j th oroder Laplacian
    %     KM = KM + (2*((c*dT(kk))^(2*j))/factorial(2*j))*Kj;
    % end
    KM = ((c*dT(kk))^2)*Lp_m(N,h,m);
    for tt=1:length(T)
        Unp = (2*speye(N,N)+KM)*Un-Unm;
        Unm = Un;
        Un = Unp;
            % if mod(tt,2)==0
            %     plot(X,Unp)
            %     ylim([-0.2 1.2]);
            %     pause(0.02);
            %     drawnow
            % end
    end
    U_all_dt(kk,:) = Unp;
    X1 = -L:(h/100):L;
    Unp1 = interp1(X,Unp,X1);
    theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
    error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
    err_M(kk) = error_val;
    % clear T error_val
end
% figure(10);
% plot(X1,theo_fn,X1,Unp1),legend("theor","simul");


%%
figure()%,subplot(2,1,1)
loglog(dT,(err_M),"LineWidth",1.4),
xlabel("\Delta t"),ylabel("Log10(error)");

U1 = U_all_dt;
U1(U1>0.51)=0;
U1(U1<-0.51)=0;
figure(2)%,subplot(2,1,2)
surf(X,dT,U1),shading flat,
xlabel("space"),ylabel("\Deltat"),
colorbar
% xlim([-1 1]),ylim([0 0.03]),zlim([-1000 1000])

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
%%
beta_p = @(p) (2^(2*p-1))*((factorial(p-1))^2)/factorial(2*p);
BB = 0;
for ii=1:m
    BB = BB + beta_p(ii);
end
BB = BB^(-1/2);
cnd = BB*h/c;
disp(cnd)