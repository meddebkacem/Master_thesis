clear variables
% close all
clc

c = 1; L = 1; TM = 0.6;
Lambda = L/10;
r = 1;                                  % order of interp r-1
m = 1;                                  % order of taylor approx

% % hList = 4*divisors(1000)/10000;
dT =  0.001; % 0.001                   % Δt fixed
N_elts = [6:5:51 101:100:1001 1501:500:2501];
hList = 2*L./N_elts;
alp = dT./hList;

error_vector = zeros(1,length(hList));
parfor jx=1:length(hList)
    h = hList(jx);
    [M,K,X] = SEM_MK(r,h,L);
    M = sparse(M);
    K = sparse(K);

    U0 = exp(-(X/Lambda).^2);               % Initial wave
    Unm = U0';
    Un = U0';
    % Unp = 0*U0';

    NN = inv(M)*K; %#ok<*MINV>
    theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
    TF = ((TM/dT)+1);
    for kt=2:TF
        Unp = 2*Un-Unm-(((c*dT)^2)*NN)*Un;
        Unm = Un;
        Un = Unp;
        if mod(kt,2) == 0
            plot(X,Unp);hold on
            plot(X,theo_fn,"--"),hold off;
            ylim([-0.1 1.1]);
            xlabel("Position (m)");
            ylabel("amplitude");
            grid on;
            drawnow;
        end
    end
    error_val = trapz(X,(theo_fn'-Unp).^2)/trapz(X,theo_fn.^2);
    error_vector(jx) = error_val;
end
%%
figure(),
% figure(1), hold on;
loglog(hList,error_vector,"*-");            % (1:end-2): truncate for r = 1
xlabel("\Delta x"),ylabel("error");
