clear variables
close all
clc

c = 1; L = 1; TM = 0.6;
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:10;                           % order of taylor approx
Nt_L = [51:5:201 201:10:2001 2201:100:5001];
dT_L = TM./Nt_L;

h = 0.01;
[M,K,X,hh,Me_L1,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L);
M = sparse(M); K = sparse(K);

Me_L = zeros(size(Me_L1));
        for kkl=2:size(Me_L1,3)-1
            Me_L(:,:,kkl) = Me_L1(:,:,kkl) ;
            Me_L(1,1,kkl) = Me_L1(1,1,kkl) + Me_L1(end,end,kkl-1);
            Me_L(end,end,kkl) = Me_L1(end,end,kkl) + Me_L1(1,1,kkl+1);
        end
        Me_L(:,:,1) = Me_L1(:,:,1);
        Me_L(end,end,1) = Me_L1(end,end,1) + Me_L1(1,1,2);
        Me_L(:,:,end) = Me_L1(:,:,end);
        Me_L(1,1,end) = Me_L1(1,1,end) + Me_L1(end,end,end-1);

U0 = exp(-(X/Lambda).^2);               % Initial wave
N0 = length(X);
NN = inv(M)*K; %#ok<*MINV>

                    tic
error_vector = zeros(length(m_L),length(dT_L));
CFL = zeros(1,length(dT_L));

for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/(h);      % c*dT/(h/10);
    error_X = zeros(1,length(m_L));

    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        NN_glob = zeros(size(K));
        KKe_L = zeros(size(Ke_L));
        for lll=1:size(Me_L,3)
            if (lll>=N1-1) && (lll<=N1+9+1)
                continue
            else
                NNe = (((c*dT)^2)/2)*inv(Me_L(:,:,lll))*Ke_L(:,:,lll);
                KKe_L(:,:,lll) = NNe;
            end
            clear NNe
        end
        for oo=1:m
            for ff=N1-1:N1+9+1
                NNe_g = inv(Me_L(:,:,ff))*Ke_L(:,:,ff);
                KKe_L(:,:,ff) = KKe_L(:,:,ff) + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NNe_g)^(oo));
                clear NNe_g
            end
        end

        for il=1:(r):length(K)-r
            NN_glob(il:il+r,il:il+r) = NN_glob(il:il+r,il:il+r) + KKe_L(:,:,(il+r-1)/r);
        end
        NN_glob = 2*sparse(NN_glob);

        for kt=2:N_t
            Unp = 2*Un - Unm + NN_glob*Un;
            Unm = Un;
            Un = Unp;
            
            % if mod(kt,5)==0
            %     plot(X,Unp),
            %     ylim([-0.1 1.1]);
            %     xlim([-L L])
            %     xlabel("Position (m)");
            %     ylabel("amplitude");
            %     grid on;
            %     % fontsize(gcf, 14, "points")
            %     drawnow;
            % end

        end
        % X1 = min(X):(h/100):max(X);
        % Unp1 = interp1(X,Unp,X1);
        % theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        % error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
        error_val = trapz(X,(theo_fn-Unp').^2)/trapz(X,theo_fn.^2);
        error_X(jm) = error_val;
    end
    error_vector(:,jt) = error_X;
end
                    toc

%% plotting: 3D

error_vector(error_vector>1000) = 10;
error_vector(isnan(error_vector)) = 1;
error_vector(isinf(error_vector)) = 1;
Ev = zeros(length(m_L)+1,length(dT_L));
Ev(1:length(m_L),1:length(dT_L)) = error_vector;
Ev(end,:) = error_vector(end,:);

figure(1)
surf(dT_L,[m_L m_L(end)+1],Ev), shading flat;
xlabel("\Delta t"), ylabel("order m"), zlabel("error"), zscale log, xscale log,% yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
view(90,-90),

%%
% error_vector = zeros(length(m_L),length(dT_L));
XBound = zeros(length(m_L),1); 
YBound = zeros(length(m_L),1);
for kk=1:length(m_L)
    idd = find(error_vector(kk,:)<1);
    if ~isempty(idd)
        YBound(kk) = idd(1);
        XBound(kk) = kk;
    end
end
XBound = nonzeros(XBound); YBound = nonzeros(YBound);

dt_f = dT_L(YBound);
cfl_f = c*dt_f/(h/10);
figure(21),hold on
plot(XBound,cfl_f,"-*",LineWidth=1.4,MarkerSize=10); grid on, hold on;
xlim([0 max(m_L)+1])
ylim([0 max(cfl_f)*1.5])
xlabel("order"), ylabel("CFL");
% fontsize(figure(21), 16, "points")


%% copied from elt_all_HO
% XBound = zeros(length(m_L),1); 
% YBound = zeros(length(m_L),1);
% err_v = error_vector(:,end:-1:1);
% for kk=1:length(m_L)
%     idd = err_v(kk,:);
%     for jkk=1:length(dT_L)
%         if idd(jkk)>=1
%             YBound(kk) = (jkk-1);
%             XBound(kk) = kk;
%             break
%         end
%     end
% end
% XBound = nonzeros(XBound); YBound = nonzeros(YBound);
% dT_L1 = dT_L(end:-1:1);
% dt_f = dT_L1(YBound);
% cfl_f = c*dt_f/(h/10);
% figure(21),hold on
% plot(XBound,cfl_f,"-*",LineWidth=1.4,MarkerSize=10); grid on, hold on;
% xlim([0 max(m_L)+1])
% ylim([0 max(cfl_f)*1.5])
% xlabel("order"), ylabel("CFL");
%%
function [M,K,X,hh,Me_L,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L)

[x_gll,w_gll] = lglnodes(r);
n_gll = length(x_gll);
Me = diag(w_gll);
D = zeros(r+1,r+1);     % elementary derivatives with gll
for iii=1:r+1
    for jjj=1:r+1
        D(iii,jjj) =alternative_dl(iii,x_gll,x_gll(jjj));
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
N_elt = round(2*L/h) +9 ;
N0 = N_elt*r+1; % round((2*L*r/h)+1);
N1 = round(0.3*(round((2*L/h)+1)))+1;

    hh1 = h.*ones(N_elt,1);
    hh1(N1:N1+9) = h/10;
hh = h.*ones(N0,1);
% N1 = round(0.3*N0)+1;
% hh(N1:N1+r-1) = h/10;

X = linspace(-L,L,N0);
M = zeros(N0); K = zeros(N0);
Me_L = zeros(n_gll,n_gll,N_elt);
Ke_L = zeros(n_gll,n_gll,N_elt);

for lk=1:N_elt
    Me_L(:,:,lk) = -hh1(lk)*Me;
    Ke_L(:,:,lk) = Ke/hh1(lk);
end


for il=1:(r):length(X)-n_gll+1
    X(il:il+r)=x_gll.*hh(il)+X(il);
    M(il:il+r,il:il+r) = M(il:il+r,il:il+r) + hh(il)*Me;
    K(il:il+r,il:il+r) = K(il:il+r,il:il+r) + Ke/hh(il);
end
M = -M;

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
end

function [x,w,P]=lglnodes(N)
% Truncation + 1
N1=N+1;
% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';
% The Legendre Vandermonde Matrix
P=zeros(N1,N1);
% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.
xold=2;
while max(abs(x-xold))>eps
    xold=x;
    P(:,1)=1;
    P(:,2)=x;
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end 
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );      
end
w=2./(N*N1*P(:,N1).^2);
% Kacem change
x=x(end:-1:1);
w=w(end:-1:1);
x = (x+1)/2;w = w/2;
end