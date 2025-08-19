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

Nel_L = 2*100;              % Number of elements in the mesh grid
h = 2*L./Nel_L;             % Mesh size

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

error_vector = zeros(length(m_L),length(dT_L));
UNP_ALL = zeros(length(dT_L),length(m_L),length(U0));
CFL = zeros(1,length(dT_L));
ll = zeros(length(Nt_L),length(m_L),size(Me_L,3));
ll_X = zeros(length(Nt_L),length(m_L));

for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/(h);      % c*dT/(h/10);
    error_X = zeros(1,length(m_L));
    UNP_ALL_X = zeros(length(m_L),length(U0));

    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        NN_glob = zeros(size(K));
        KKe_L = zeros(size(Ke_L));
        N1g = (round(size(Me_L,3)*0.3)+1);
        % the following loop calculates the matrices for leapfrog per element
        % except for the smaller element that is done up to the order m of 
        % the modified equation in the loop below
        for lll=1:size(Me_L,3)
            if lll==N1g
                continue
            else
                NNe = (((c*dT)^2)/2)*inv(Me_L(:,:,lll))*Ke_L(:,:,lll);
                KKe_L(:,:,lll) = NNe;
            end
            clear NNe
        end
        for oo=1:m
            NNe_g = inv(Me_L(:,:,N1g))*Ke_L(:,:,N1g);
            KKe_L(:,:,N1g) = KKe_L(:,:,N1g) + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NNe_g)^(oo));
            clear NNe_g
        end
% the following 2 loops combine and give us a global matrix used for the scheme
        for il=1:(r):length(K)-r
            NN_glob(il:il+r,il:il+r) = NN_glob(il:il+r,il:il+r) + KKe_L(:,:,(il+r-1)/r);
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
        UNP_ALL_X(jm,:) = Unp;
    end
    error_vector(:,jt) = error_X;
    UNP_ALL(jt,:,:) = UNP_ALL_X;
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
ylim([1 m_L(end)+0.99])
%% max stabl CFL curve per order
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
xlim([1 max(m_L)])
ylim([0 max(cfl_f)*1.5])
xlabel("order"), ylabel("CFL");
