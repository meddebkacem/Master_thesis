clear variables
% close all
% clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% * To confirm the power assembly works correctly and similarly to the 
%%%%% case of doing the power over the global matrix
%%%%% * Compare it to the results of "SEM_heteor_time_HO.m", once valid you
%%%%% can also use/show the results of "SEM_heteor_elt_HO.m" that makes the
%%%%% power matrix only on the small element, and the obtained results of 
%%%%% "SEM_heteor_elt_all_HO.m" where we do it on a small surrounding of
%%%%% the small element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 1; L = 1; TM = 0.6;%*(5/4);
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:5 ;                           % order of taylor approx

% Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
        Nt_L = [51:5:201 201:10:2001 2201:100:5001];
        % Nt_L = Nt_L([65:68 76:93 117:126]);
                % Nt_L = [1401:10:2001 2101:50:3601];
                % Nt_L = [2701:10:3601];
                % Nt_L = [2901 2971 2981 3381];

N__el = 2*100;
h = 2*L./N__el;

dT_L = TM./Nt_L;

[M,K,X,hh,Me_L1,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L);

M = sparse(M); K = sparse(K);

% NN = inv(M)*K;

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

KKe_L = zeros(size(Ke_L));
for lll=1:size(Me_L,3)
    NNe_g = inv(Me_L(:,:,lll))*Ke_L(:,:,lll);
    KKe_L(:,:,lll) = KKe_L(:,:,lll) + NNe_g;
    clear NNe_g
end
KK_glob = zeros(size(K));
for il=1:(r):length(K)-r
    KK_glob(il:il+r,il:il+r) = KK_glob(il:il+r,il:il+r) + KKe_L(:,:,(il+r-1)/r);
end

[KK_part,KK_sq] = square_parts_KM(Me_L,Ke_L,r,N__el);
KK_pow_2 = KK_part+KK_sq;

                    tic
error_vector = zeros(length(m_L),length(dT_L));
UNP_ALL = zeros(length(dT_L),length(m_L),length(U0));
CFL = zeros(1,length(dT_L));
ll = zeros(length(Nt_L),length(m_L),size(Me_L,3));
ll_X = zeros(length(Nt_L),length(m_L));
% figure,
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

        switch m
            case 1
                NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob;
            case 2
                NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob ...
                    + ((c*dT)^(4)/factorial(4))*KK_pow_2;
            case 3
                NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob ...
                    + ((c*dT)^(4)/factorial(4))*KK_pow_2 ...
                    + ((c*dT)^(6)/factorial(6))*KK_pow_2*KK_glob;
            case 4
                NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob ...
                    + ((c*dT)^(4)/factorial(4))*KK_pow_2 ...
                    + ((c*dT)^(6)/factorial(6))*KK_pow_2*KK_glob ...
                    + ((c*dT)^(8)/factorial(8))*KK_pow_2*KK_pow_2;
            case 5
                NN_glob = ((c*dT)^(2)/factorial(2))*KK_glob ...
                    + ((c*dT)^(4)/factorial(4))*KK_pow_2 ...
                    + ((c*dT)^(6)/factorial(6))*KK_pow_2*KK_glob ...
                    + ((c*dT)^(8)/factorial(8))*KK_pow_2*KK_pow_2 ...
                    + ((c*dT)^(10)/factorial(10))*KK_pow_2*KK_pow_2*KK_glob;
        end



        NN_glob = 2*sparse(NN_glob);
% local application of power 2

        
            % figure,spy(KK_pow_2,"ro")
            % figure,spy(KK_glob^2)
            % figure,spy(KK_pow_2-KK_glob^2)
            % AAA = (KK_pow_2-KK_glob^2);
            % AAA(abs(AAA)<10^(-5))=0;

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
fontsize(figure(21), 16, "points")

%%

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

