clear variables
% close all
clc

c = 1; L = 1; TM = 0.6;
Lambda = L/10;
or= 2;                                  % order approx for FD
r = 2;                                  % order of interp r-1 for SEM
m = 10;                                 % order of taylor approx

N_elts = 200;
h = 2*L/N_elts;
N_t = [11:2:101 151 201:100:1001];
dTList = TM./N_t;

    % first for SEM, second for FD high order approx
[M,K,X] = SEM_MK(r,h,L);
% [M,K,X] = FD_MK(or,h,L);

M = sparse(M);
K = sparse(K);
U0 = exp(-(X/Lambda).^2);               % Initial wave

error_vector = zeros(m,length(dTList));
for jt=1:length(dTList)
    for mm=1:m
        dT = dTList(jt);
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
    
        NN = inv(M)*K; %#ok<*MINV>
        KK = sparse(zeros(size(K)));
        for oo=2:mm
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
        
        for kt=2:N_t(jt)
            Unp = 2*Un-Unm+( ((c*dT)^2)*NN + 2*KK)*Un;      % FD
            Unm = Un;
            Un = Unp;
            % if mod(kt,2) == 0
            %     plot(X,Unp); % hold on
            %     % plot(X,theo_fn,"--"),hold off;
            %     ylim([-0.1 1.1]);
            %     xlabel("Position (m)");
            %     ylabel("amplitude");
            %     grid on;
            %     drawnow;
            % end
        end
        X1 = -L:(h/100):L;
        Unp1 = interp1(X,Unp,X1);
        theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        error_vector(mm,jt) = error_val;
    end
end
%% plotting: 2D 2nd & 4th order

error_vector(error_vector>1) = 1;

figure(),
semilogy(dTList,error_vector(1,:),"*-");hold on;
semilogy(dTList,error_vector(2,:),"*-");hold off;
legend("2^{nd} order","4^{nd} order")
xlabel("dt"),ylabel("error");

figure()
surf(dTList,(1:1:m),error_vector), shading flat;
xlabel("dt"), ylabel("order"), zlabel("error"), zscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);

%% plotting for multiple orders

max_t = zeros(1,m);
% error_vector = zeros(m,length(dTList));
for jk=1:m
    M_err = find(error_vector(jk,:)<0.0025); % <0.02
    max_t(jk) = max(dTList(M_err));
end

figure()
plot((2:2:2*m),max_t/h,"-*")
xlim([0 2*m+1])
ylim([0 max(max_t/h)+.3])
xlabel("order"),ylabel("max dt")
