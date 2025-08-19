clear variables
% close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = 1; N_t = [31:10:101 201:100:1001 2001:1000:10001];
% r = 2; N_t = [1001 2001:2000:10001 20001:20000:100001 200001:200000:1000001];
% r = 3; N_t = [100001 200001:200000:1000001 2000001:2000000:10000001];
c = 1; L = 1; TM = 0.6;
Lambda = L/10;
r = 3;                                  % order of interp r-1
m = 1;                                  % order of taylor approx

N_elts = 1000;
h = 2*L/N_elts;
N_t = [100001 200001:200000:1000001 2000001:2000000:10000001];
dTList = TM./N_t;

[M,K,X] = SEM_MK(r,h,L);
M = sparse(M);
K = sparse(K);
U0 = exp(-(X/Lambda).^2);               % Initial wave
                                tic
error_vector = zeros(1,length(dTList));
parfor jt=1:length(dTList)
    dT = dTList(jt);
    Unm = U0';
    Un = U0';
    % Unp = 0*U0';

    NN = inv(M)*K; %#ok<*MINV>
    theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);
    % TF = (TM/dT);
    for kt=2:N_t(jt)
        Unp = 2*Un-Unm-(((c*dT)^2)*NN)*Un;
        Unm = Un;
        Un = Unp;
        % if mod(kt,2) == 0
        %     plot(X,Unp);hold on
        %     plot(X,theo_fn,"--"),hold off;
        %     ylim([-0.1 1.1]);
        %     xlabel("Position (m)");
        %     ylabel("amplitude");
        %     grid on;
        %     drawnow;
        % end
    end
    error_val = trapz(X,(theo_fn'-Unp).^2)/trapz(X,theo_fn.^2);
    error_vector(jt) = error_val;
end
                                toc
%%
% figure(),
% % figure(1), hold on;
% loglog(c.*dTList/h,error_vector,"*-");   % (1:end-2): maybe need truncate 4 r=1
% xlabel("CFL \alpha"),ylabel("error");

figure(),
% figure(1), hold on;
loglog(N_t,error_vector,"*-");   % (1:end-2): maybe need truncate 4 r=1
xlabel("nb elts"),ylabel("error");


%%

% % % clear
% % % close
% % % clc
% % % 
% % % c = 1; L = 1; TM = 0.6;
% % % 
% % % N_elts = 1000; % N_elts = 101;
% % % h = 2*L/N_elts;
% % % 
% % % 
% % % N_t = [301:50:1001 ... 
% % %     2001:1000:10001 ...
% % %     20001:20000:100001 ...
% % %     200001:200000:1000001 ...
% % %     2000001:2000000:10000001 ...
% % %     20000001:20000000:100000000];
% % % dTList = TM./N_t;
% % % 
% % % cfl = c*dTList./h;