clear variables
% close all
clc

c = 1; L = 1; TM = 0.6;%*(5/4);
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:10;                           % order of taylor approx

% Nel_L = [20:10:100 120:20:200 250:50:2000];

% Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];     % 2251:250:9001
% Nt_L = [51:5:201 201:10:601]; % t confirm r=3 normal mesh
% Nt_L = [1601:10:2001 2201:200:5001];   % for mesh with tiny element
    % Nt_L = Nt_L([12:16 51:55]);
        Nt_L = [51:5:201 201:10:2001 2201:100:5001];
        % Nt_L = [2601 1501 1901 1121 931];

% Nt_L = 1350;
Nel_L = 2*100;
h_L = 2*L./Nel_L;

% Nt_L = 36*(5/4)*30; % 36*(5/4)*120;
dT_L = TM./Nt_L;

h = h_L;
[M,K,X,hh,N1,Me_S,Ke_S] = SEM_MK_hetero(r,h,L);

                                        MMg = M + zeros(size(M));
                                        KKg = K + zeros(size(M));
                                        dt_critic = 2/sqrt(max(eig(KKg,-MMg)));
                                        disp(dt_critic/(h/10))

M = sparse(M);
K = sparse(K);
U0 = exp(-(X/Lambda).^2);               % Initial wave

N0 = length(X);
ZZ = cumsum(hh)-1;

NN = inv(M)*K; %#ok<*MINV>
NNe = inv(Me_S)*Ke_S;

NN_4_pow1 = NN(N1:N1+r,:);
                    tic

Unp_L = zeros(length(U0),length(m_L),length(dT_L));
error_vector = zeros(length(m_L),length(dT_L));
ll = zeros(length(m_L),length(dT_L));
CFL = zeros(1,length(dT_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/(h);      % c*dT/(h/10);
    error_X = zeros(1,length(Nel_L));
    ll_X = zeros(1,length(Nel_L));
    Unp_L_M = zeros(length(U0),length(m_L));
    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        KK = zeros(size(K));
        KKe = zeros(size(Ke_S));
                for oo=2:m        % for high order Modified equation overe all mesh
                    KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
                end
        % NN1 = zeros(size(KK));
        % NN1(N1:N1+r,N1:N1+r) = inv(Me_S)*Ke_S;
        % for oo=2:m
        %     KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN1)^(oo));
        %     KKe = KKe +( (c*dT)^(2*oo)/factorial(2*oo) )*((NNe)^(oo));
        % end
        KKe_eigs =  ((c*dT)^2)*NNe + 2*KKe;
        ll_X(jm) = 2*sqrt(3)/max(eigs(KKe_eigs));
                    % KK_temp = zeros(size(K));
                    % KK_temp(N1:N1+r,:) = NN_4_pow1; % NN(N1:N1+r,:);
                    % KK_temp = sparse(KK_temp);
                    % for oo=2:m
                    %     KK_temp = KK_temp*NN;
                    %     KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*(KK_temp);
                    % end
        KK = sparse(KK);
        for kt=2:N_t
            Unp = 2*Un-Unm+( ((c*dT)^2)*NN + 2*KK )*Un;
            Unm = Un;
            Un = Unp;
            
            % if mod(kt,10)==0
            %     plot(X,Unp), % hold on;
            %     % plot(ZZ,zeros(N0,1),"--k"),
            %     % plot([ZZ(299) ZZ(300)],[0 0],"*-r")
            %     % hold off
            %     ylim([-0.1 1.1]);
            %     xlim([-L L])
            %     xlabel("Position (m)");
            %     ylabel("amplitude");
            %     grid on;
            %     % fontsize(gcf, 14, "points")
            %     drawnow;
            % end

        end
        Unp_L_M(:,jm) = Unp;
        % X1 = min(X):(h/100):max(X);
        % Unp1 = interp1(X,Unp,X1);
        % theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        % error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
        error_val = trapz(X,(theo_fn-Unp').^2)/trapz(X,theo_fn.^2);
        error_X(jm) = error_val;
    end
    error_vector(:,jt) = error_X;
    ll(:,jt) = ll_X;
    Unp_L(:,:,jt) = Unp_L_M;
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

%% old one
% XBound = zeros(length(m_L),1); 
% YBound = zeros(length(m_L),1);
% for kk=1:length(m_L)
%     idd = find(error_vector(kk,:)<1);
%     if ~isempty(idd)
%         YBound(kk) = idd(1);
%         XBound(kk) = kk;
%     end
% end
% XBound = nonzeros(XBound); YBound = nonzeros(YBound);
% 
% dt_f = dT_L(YBound);
% cfl_f = c*dt_f/(h/10);
% figure(21),hold on
% plot(XBound,cfl_f,"-*",LineWidth=1.4,MarkerSize=10); grid on, hold on;
% xlim([0 max(m_L)+1])
% ylim([0 max(cfl_f)*1.5])
% xlabel("order"), ylabel("CFL");
% fontsize(figure(21), 16, "points")

%% new one

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
xlim([0 max(m_L)+1])
ylim([0 max(cfl_f)*1.5])
xlabel("order"), ylabel("CFL");
% fontsize(figure(21), 16, "points")

