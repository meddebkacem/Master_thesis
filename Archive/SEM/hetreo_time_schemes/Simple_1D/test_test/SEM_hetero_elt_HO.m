clear variables
close all
clc

c = 1; L = 1; TM = 0.6;%*(5/4);
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:10;                           % order of taylor approx

% Nel_L = [20:10:100 120:20:200 250:50:2000];

% Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];     % 2251:250:9001
% Nt_L = [51:5:201 201:10:601]; % t confirm r=3 normal mesh
% Nt_L = [1601:10:2001 2201:200:5001];   % for mesh with tiny element
        Nt_L = [51:5:201 201:10:2001 2201:100:5001];
                % Nt_L = [1401:10:2001 2101:50:3601];
                % Nt_L = [2701:10:3601];
                % Nt_L = [2901 2971 2981 3381];
                                    % Nt_L = Nt_L([76:93 117:126]);

Nel_L = 2*100;
h_L = 2*L./Nel_L;

dT_L = TM./Nt_L;

h = h_L;
[M,K,X,hh,Me_L1,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L);

M = sparse(M); K = sparse(K);

                                MMg = M + zeros(size(M));
                                KKg = K + zeros(size(M));
                                dt_critic = 2/sqrt(max(eig(KKg,-MMg)));
                                disp(dt_critic/(h/10))

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
                % ll(jt,jm,lll) = sqrt(12/max(eigs(KKe_L(:,:,lll) )));
                % disp(sqrt(12/max(eigs(KKe_L(:,:,lll)))))
            clear NNe_g
        end
                    % for oo=1:m
                    %     for llk=1:size(Me_L,3)
                    %         NNe1 = inv(Me_L(:,:,llk))*Ke_L(:,:,llk);
                    %         KKe_L(:,:,llk) = KKe_L(:,:,llk) + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NNe1)^(oo));
                    %     end
                    % end
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
        UNP_ALL_X(jm,:) = Unp;
    end
    error_vector(:,jt) = error_X;
    UNP_ALL(jt,:,:) = UNP_ALL_X;
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

% cfreq = 5;
% N_time = 300;
% dtt = 0.6/N_time;
% ds = 0.005;
% xM = 2;
% yM = 3;
% xl = -xM:ds:xM;
% yl = -yM:ds:yM;
% [x,y] = meshgrid(xl,yl);
% xs=0; ys=0;
% 
% a=8*cfreq/sqrt(pi);
% t=((1:N_time)/(1/dtt)-4/a); 
% figure,
% % % for ik=1:N_time
% % %     d = (x-xs).^2+(y-ys).^2;
% % %     w=-(exp(-a^2*(d)/2).*(a^2*(d)-1));
% % %     surf(x,y,w), shading flat;
% % %     pause(0.02)
% % %     drawnow
% % % end

% tut = zeros(size(KK));
% tut(N1:N1+r,:) = KK(N1:N1+r,:);
% tut(:,N1:N1+r) = KK(:,N1:N1+r);
% tut1 = zeros(size(KK));
% tut1(N1:N1+r,N1:N1+r) = KK(N1:N1+r,N1:N1+r);
% figure,spy(tut);hold on, spy(tut1,"ro")

