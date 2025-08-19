clear variables
% close all
clc
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tic
% A: FD 4th order + Leap-frog
% B: FD 4th order + Modified Equation (2 derivation orders)
% C: SEM 4th order+ Leap-frog
% D: SEM 4th order+ Modified Equation (2 derivation orders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % PROOFS % % % % % %
% A: 99/355
% B: (6.32) 101/355 - (2.11) 20/43 + (B.1) 42/43 [6]
% C: Table (11.2) 197/355
% D: for r=2 : (11.65) 198/355 - for r=3 : (11.66) 198/355
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 1; L = 1; TM = 0.6;
Lambda = L/10;
or= 2;                                  % order approx for FD
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:10;                                  % order of taylor approx

Nt_L = [31:2:99 101:10:291 301:20:1001 1101:50:2001];
    % Nt_L = [11:2:99 101:20:201 251 301:100:1001 1501:500:5001];
dT_L = TM./Nt_L;


    % first for SEM, second for FD high order approx
Nel_L = 2*100; h = 2*L./Nel_L;            % use 200 for SEM
[M,K,X] = SEM_MK(r,h,L);
        % Nel_L = 1000; h = 2*L./Nel_L;       % use 1000 for FD
        % [M,K,X] = FD_MK(or,h,L);

M = sparse(M);
K = sparse(K);
U0 = exp(-(X/Lambda).^2);               % Initial wave

                    tic
error_vector = zeros(length(m_L),length(dT_L));
CFL = zeros(1,length(dT_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/h;

    error_X = zeros(1,length(Nel_L));
    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        NN = inv(M)*K; %#ok<*MINV>
        KK = sparse(zeros(size(K)));
        for oo=1:m
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
        
        for kt=2:N_t
            Unp = 2*Un-Unm+(2*KK)*Un;
            Unm = Un;
            Un = Unp;
        end
        X1 = -L:(h/100):L;
        Unp1 = interp1(X,Unp,X1);
        theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        error_X(jm) = error_val;
    end
    error_vector(:,jt) = error_X;
end
                    toc
%% plotting: 3D

error_vector(error_vector>1) = 1;
error_vector(isnan(error_vector)) = 1;

figure(1)
surf(dT_L,[m_L 11],[error_vector ;error_vector(10,:)]), shading flat;
xlabel("\Delta t"), ylabel("order m"), zlabel("error"), zscale log, xscale log,% yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
% view(2),
view(90,-90),
yticks([m_L 11])
ylim([1 10.99])
% xlim([min(dT_L) max(dT_L)])


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
figure(21),hold on
plot(XBound,c*dt_f/h,"-*",LineWidth=1.4,MarkerSize=10); grid on, hold on;
% plot(XBound(1:1:3),c*dt_f(1:1:3)/h,"o",LineWidth=2,MarkerSize=14);
% ylim([0 floor(round(max(10*c*dt_f/h))/10 + 1)]), 
ylim([0 round(max(c*dt_f/h))+1])
xlabel("order"), ylabel("CFL");
% fontsize(figure(21), 16, "points")


% figure(123),
% pcolor(dT_L,m_L,error_vector), hold on;
% xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error");
% shading flat, xscale log,yscale log;
% colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);

% % % figure(),plot(X,U0,X,Unp)
% % % phase = -unwrap(angle(fft(Unp',10000))-angle(fft(U0,10000)));
% % % figure(),plot(1:1:10000,phase)
