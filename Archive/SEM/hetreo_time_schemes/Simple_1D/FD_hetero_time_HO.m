clear variables
% close all
clc
        tic

c = 1; L = 1; TM = 0.6*(4/3);
Lambda = L/10;
or= 2;                                  % order approx for FD
m_L = 20;%1:1:10;                       % order of taylor approx

% Nel_L = [20:10:100 120:20:200 250:50:2000];
% Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
Nel_L = 100;
h_L = 2*L./Nel_L;

Nt_L = 36*(4/3)*30;
dT_L = TM./Nt_L;

h = h_L;
[M,K,X,hh,N1] = FD_MK_hetero(or,h,L);        % verified
M = sparse(M);
K = sparse(K);
U0 = exp(-(X/Lambda).^2);               % Initial wave

N0 = length(X);

NN = inv(M)*K; %#ok<*MINV>

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
        KK = zeros(size(K));
        % for oo=2:m        % for high order Modified equation overe all mesh
        %     KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        % end
                KK_1 = NN(N1,:);
                KK_temp = NN(N1,:);
                for oo=2:m
                    for ii=N1-or*oo:1:N1+or*oo
                        KK_shift = circshift(KK_1,ii-N1)';
                        if ii<=0
                            ii = N0+ii;
                        end
                        KK_temp(ii) = KK_temp(ii) + ( (c*dT)^(2*oo)/factorial(2*oo))*(KK_temp*KK_shift);
                    end
                end

        NN(N1,:) = NN(N1,:)+ KK_temp;

        figure(123)
        for kt=2:N_t
            Unp = 2*Un-Unm+( ((c*dT)^2)*NN + 2*KK )*Un;
            Unm = Un;
            Un = Unp;
            
            if mod(kt,19)==1
                plot(X,Unp),
                ylim([-0.1 1.1]);
                xlim([-L L])
                xlabel("Position (m)");
                ylabel("amplitude");
                grid on;
                % fontsize(gcf, 14, "points")
                drawnow;
            end

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

