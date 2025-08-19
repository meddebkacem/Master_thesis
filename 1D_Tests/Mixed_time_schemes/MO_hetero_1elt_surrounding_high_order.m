clear variables
close all
clc

% Initial parameters
c = 1; L = 1; TM = 0.6;
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:10;                           % List of orders of the modified equation

Nt_L = [51:5:201 201:10:2001 2201:100:5001];        % List of number of iterations
dT_L = TM./Nt_L;                                    % List of time steps

Nelts = 2*100;              % Number of elements in the mesh grid
h = 2*L./Nelts;             % Mesh size

% function to make the M & K matrices with the small element,
% the output: the M & K matrices, the nodes vector X, the Me_L1 & Ke_L1
% that are 3 matrices containing the M & K matrices per element respectively
[M,K,X,Me_L1,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L);

M = sparse(M); K = sparse(K);

MMg = M + zeros(size(M));
KKg = K + zeros(size(M));
dt_critic = 2/sqrt(max(eig(KKg,-MMg)));
disp(append("max CFL for first order:",num2str(dt_critic/(h/10))))


Me_L = zeros(size(Me_L1));
% correcting elemental mass matrices with connectivity to surrounding elements
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
% calculating the M^(-1)K matrix for each element
KKe_L = zeros(size(Ke_L));
for lll=1:size(Me_L,3)
    NNe_g = inv(Me_L(:,:,lll))*Ke_L(:,:,lll);
    KKe_L(:,:,lll) = KKe_L(:,:,lll) + NNe_g;
    clear NNe_g
end
% making the global M^(-1)K matrix with previously calvulated ones for the
% leapfrog term (the first term of the modified equation series)
KK_glob = zeros(size(K));
for il=1:(r):length(K)-r
    KK_glob(il:il+r,il:il+r) = KK_glob(il:il+r,il:il+r) + KKe_L(:,:,(il+r-1)/r);
end

% calculating element by element the local (M^(-1)K)^2 matrix around the small element
N1g = (round(size(Me_L,3)*0.3)+1);
Ne = 1;         % max power exponent - 1
ke_el = Ke_L(:,:,N1g-Ne:N1g+Ne);
me_el = Me_L(:,:,N1g-Ne:N1g+Ne);
Np = (2*Ne+1);
ke_el_1 = zeros(Np*r+1);
me_el_1 = zeros(Np*r+1);
for lk=1:r:(Np*r+1)-r
    ke_el_1(lk:lk+r,lk:lk+r) = ke_el_1(lk:lk+r,lk:lk+r)+ke_el(:,:,(lk+r-1)/r);
    me_el_1(lk:lk+r,lk:lk+r) = me_el_1(lk:lk+r,lk:lk+r)+me_el(:,:,(lk+r-1)/r);
end

[KK_part,KK_sq] = square_parts_KM(me_el,ke_el,r,Np);
KK_pow_2_el = KK_part+KK_sq;

error_vector = zeros(length(m_L),length(dT_L));
CFL = zeros(1,length(dT_L));
ll = zeros(length(Nt_L),length(m_L),size(Me_L,3));
ll_X = zeros(length(Nt_L),length(m_L));

for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/(h);
    error_X = zeros(1,length(m_L));

    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        
        NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob;
        NN_el_1 = inv(me_el_1)*ke_el_1;
        % calculating the modified equation matrix locally around the small element
        % (starting from m=2, as the m=1 is part of the global matrix previously calculated)
        for oo=2:m
            Pw_N = Nth_power(KK_pow_2_el,NN_el_1,oo);
            NN_glob(N1-Ne*r:N1+(Ne+1)*r,N1-Ne*r:N1+(Ne+1)*r) = ...
                NN_glob(N1-Ne*r:N1+(Ne+1)*r,N1-Ne*r:N1+(Ne+1)*r) ....
                + ((c*dT)^(2*oo)/factorial(2*oo))*Pw_N;
        end

        NN_glob = 2*sparse(NN_glob);
        for kt=2:N_t
            Unp = 2*Un - Unm + NN_glob*Un;
            Unm = Un;
            Un = Unp;
        end
        theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
        error_val = trapz(X,(theo_fn-Unp').^2)/trapz(X,theo_fn.^2);
        error_X(jm) = error_val;
    end
    error_vector(:,jt) = error_X;
end
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
yticks([m_L m_L(end)+1]);
ylim([1 m_L(end)+0.99]);
%% max stable CFL per order
XBound = zeros(length(m_L),1); 
YBound = zeros(length(m_L),1);
err_v = error_vector(:,end:-1:1);
for kk=1:length(m_L)
    idd = err_v(kk,:);
    for jkk=1:length(dT_L)
        if idd(jkk)>=1
            YBound(kk) = (jkk-1);
            XBound(kk) = kk;
            break
        end
    end
end
XBound = nonzeros(XBound); YBound = nonzeros(YBound);
dT_L1 = dT_L(end:-1:1);
dt_f = dT_L1(YBound);
cfl_f = c*dt_f/(h/10);
figure(21),hold on
plot(XBound,cfl_f,"-*",LineWidth=1.4,MarkerSize=10); grid on, hold on;
xlim([1 max(m_L)])
ylim([0 floor(10*max(cfl_f))*0.15])
xlabel("order"), ylabel("CFL");

%% function used (only for these cases)
function [Pw_N] = Nth_power(KK_pow_2,KK_pow_1,n_pow)
switch mod(n_pow,2)
    case 0
        M = KK_pow_2^(floor(n_pow/2));
    case 1
        M = KK_pow_2^(floor(n_pow/2))*KK_pow_1;
end
Pw_N = M;
end

function [KK_part,KK_sq] = square_parts_KM(Me_L,Ke_L,r,N_elts)
N = N_elts*r+1;
KK_part = zeros(N);
KK_sq = zeros(N);
A = inv(Me_L(:,:,1))*Ke_L(:,:,1);
KK_sq(1:r+1,1:r+1) = A^2;

for jj=2:N_elts
    indx = (jj-1)*r+1;
    A = inv(Me_L(:,:,jj-1))*Ke_L(:,:,jj-1);
    B = inv(Me_L(:,:,jj))*Ke_L(:,:,jj);
    KK_sq(indx:indx+r,indx:indx+r) = KK_sq(indx:indx+r,indx:indx+r) + B^2;

    K12= A(1:end-1,4);
    K21= A(4,1:end-1);
    K2 = A(end,end);
    
    K3 = B(1,1);
    K34= B(1,2:end);
    K43= B(2:end,1);
    
    KK_part(indx-r:indx,indx:indx+r) = KK_part(indx-r:indx,indx:indx+r) + [K12;K2]*[K3 K34];
    KK_part(indx:indx+r,indx-r:indx) = KK_part(indx:indx+r,indx-r:indx) + [K3;K43]*[K21 K2];

    % KK_part(1:4,4:7) = KK_part(1:4,4:7) + [K12;K2]*[K3 K34];
    % KK_part(4:7,1:4) = KK_part(4:7,1:4) + [K3;K43]*[K21 K2];
    clear K12 K21 K2 K3 K34 K43
end
KK_part = sparse(KK_part);
KK_sq = sparse(KK_sq);
end

