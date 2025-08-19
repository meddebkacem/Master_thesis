clear variables
close all
clc

% Initial parameters
c = 1; L = 1; TM = 0.6;
Lambda = L/10;
or= 2;                      % order approx for FD
r = 3;                      % order of interp r-1 for SEM
m = 2;                      % order of modified equation (taylor approx)

Nelt_L = [10:10:100 120:20:200 250:50:2000];                        % List of number of total elements
h_L = 2*L./Nelt_L;                                                  % List of mesh size
Nt_L = [11:2:99 101:20:201 251 301:100:1001 1501:500:5001];         % List of number of iterations
dT_L = TM./Nt_L;                                                    % List of time steps

error_vector = zeros(length(h_L),length(dT_L));
CFL = zeros(length(h_L),length(dT_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    error_X = zeros(1,length(Nelt_L));
    CFL_X = zeros(1,length(Nelt_L));
    for jx=1:length(Nelt_L)
        h = h_L(jx);
        CFL_X(jx) = c*dT/h;

% Select the space scheme to test by commenting one of the following 2 lines and keeping the other.
% The SEM_MK for spectral elements and FD_MK for finite difference (high order)
% the output is: M & K matrices and the nodes vector X
        % [M,K,X] = SEM_MK(r,h,L);
        [M,K,X] = FD_MK(or,h,L);
        
        M = sparse(M);
        K = sparse(K);
        U0 = exp(-(X/Lambda).^2);               % Initial wave
        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        NN = inv(M)*K; %#ok<*MINV>
        KK = sparse(zeros(size(K)));
% the loop calculates the high order of the modified equation to be added
% to the initial term of the leapfrog time schemes
        for oo=2:m
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
        
        for kt=2:N_t
            Unp = 2*Un-Unm+(((c*dT)^2)*NN + 2*KK)*Un;
            Unm = Un;
            Un = Unp;
        end
        X1 = -L:(h/100):L;
        Unp1 = interp1(X,Unp,X1);
        theo_fn = 0.5*exp(-((X1-c*TM)/Lambda).^2) + 0.5*exp(-((X1+c*TM)/Lambda).^2);
        error_val = trapz(X1,(theo_fn-Unp1).^2)/trapz(X1,theo_fn.^2);
        error_X(jx) = error_val;
    end
    error_vector(:,jt) = error_X;
    CFL(:,jt) = CFL_X;
end
%% plotting: 3D
error_vector(error_vector>1) = 1;
error_vector(isnan(error_vector)) = 1;

figure()
surf(dT_L,h_L,error_vector), shading flat;
xlabel("\Delta t"), ylabel("h"), zlabel("error"), zscale log;         xscale log, yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
view(2),

%% plot 2D error dt
figure()
loglog(dT_L,error_vector(51,:),"*-"),hold on;
loglog(dT_L,error_vector(23,:),"*-");
loglog(dT_L,error_vector(11,:),"*-");
loglog(dT_L,error_vector(06,:),"*-");
xlim([10^(-4) 10^(-1)]),ylim([10^(-8) 10^(0)]);
xlabel("\Delta t"),ylabel("error"),grid on;

% xlines = 100*[dT_L(60) dT_L(64)];
% ylines = 100*[error_vector(51,60) error_vector(51,64)];
    xlines = 10*[dT_L(62) dT_L(68)];
    ylines = 10*[error_vector(23,62) error_vector(23,68)];
    % ylines = 100*[error_vector(50,60) error_vector(50,64)];
loglog([xlines(1) xlines(2)],[ylines(2) ylines(2)],"k",linewidth=1.5)
loglog([xlines(1) xlines(1)],[ylines(2) ylines(1)],"k",linewidth=1.5)
loglog(xlines,ylines,"k",linewidth=1.5)
slp1 = (log10(ylines(1))-log10(ylines(2)))/(log10(xlines(1))-log10(xlines(2)));
text(0.5*(xlines(1)+xlines(2)),ylines(2)/2,"1",FontSize=14)
text(round(1000*xlines(1)+1)/1000,(ylines(1)+3*ylines(2))/4,num2str(round(slp1)),FontSize=14)

legend(append("h=",num2str(h_L(51))), append("h=",num2str(h_L(23))),...
    append("h=",num2str(h_L(11))),append("h=",num2str(h_L(06))),"","","",Location="southeast");
%% plot 2D error dx
figure()
loglog(h_L,error_vector(:,68),"*-"),hold on;
loglog(h_L,error_vector(:,52),"*-");
loglog(h_L,error_vector(:,25),"*-");
loglog(h_L,error_vector(:,05),"*-");
xlim([10^(-3) 2*10^(-1)]),ylim([10^(-8) 10^(0)]);
xlabel("\Delta x"),ylabel("error"), grid on;

% % for FD
% xl = 2*[h_L(10) h_L(15)];
% yl = 2*[0.5*error_vector(10,68) error_vector(15,68)];
%     xl1 = 2*[h_L(1) h_L(2)];
%     yl1 = 2*[error_vector(1,68) error_vector(2,68)];

% for SEM
xl = 2*[h_L(4) h_L(8)];
yl = 2*[error_vector(4,68) error_vector(8,68)];
      % xl1 = xl; yl1 = yl;
    xl1 = 2*[h_L(1) h_L(5)];
    yl1 = 2*[error_vector(1,68) error_vector(5,68)];

loglog([xl(1) xl(2)],[yl(2) yl(2)],"k",linewidth=1.5)
loglog([xl(1) xl(1)],[yl(2) yl(1)],"k",linewidth=1.5)
loglog([xl(1) xl(2)],yl,"k",linewidth=1.5)
slp2 = (log10(yl1(1))-log10(yl1(2)))/(log10(xl1(1))-log10(xl1(2)));
text(0.5*(xl(1)+xl(2)),yl(2)/2,"1",FontSize=14)
text(round(100*xl(1)+1)/100,(yl(1)+3*yl(2))/4,num2str(round(slp2)),FontSize=14)

legend(append("\Delta t=",num2str(dT_L(68))), append("\Delta t=",num2str(dT_L(52))),...
    append("\Delta t=",num2str(dT_L(25))),append("\Delta t=",num2str(dT_L(05))),"","","",Location="southeast");

%%
EV = error_vector;
EV(EV<1) = 0;
EV = -EV + 1;
FFF = EV.*CFL;
VV = zeros(length(h_L),1);
for kkk=1:length(h_L)
    UU = find(FFF(kkk,:)>0);
    VV(kkk) = UU(1);
end
mean(dT_L(VV)./h_L)

%% Average CFL from result
XBound = zeros(length(dT_L),1); 
YBound = zeros(length(dT_L),1);
for kk=1:length(dT_L)
    idd = find(error_vector(:,kk)<1);
    if ~isempty(idd)
        YBound(kk) = idd(end);
        XBound(kk) = kk;
    end
end
XBound = nonzeros(XBound); YBound = nonzeros(YBound);

dt_f = dT_L(XBound); dx_f = h_L(YBound);
alp = mean(dt_f(2:end)-dt_f(1:end-1))/mean(dx_f(2:end)-dx_f(1:end-1));
disp(alp)