clear variables
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
m = 2;                                  % order of taylor approx

Nel_L = [20:10:100 120:20:200 250:50:2000];
% Nel_L = 400;
h_L = 2*L./Nel_L;
Nt_L = [11:2:99 101:20:201 251 301:100:1001 1501:500:5001];
% Nt_L = 600;
dT_L = TM./Nt_L;
                    tic
error_vector = zeros(length(h_L),length(dT_L));
CFL = zeros(length(h_L),length(dT_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);

    error_X = zeros(1,length(Nel_L));
    CFL_X = zeros(1,length(Nel_L));
    for jx=1:length(Nel_L)
        h = h_L(jx);
        CFL_X(jx) = c*dT/h;

            % first for SEM, second for FD high order approx
        [M,K,X] = SEM_MK(r,h,L);
        % [M,K,X] = FD_MK(or,h,L);
        M = sparse(M);
        K = sparse(K);
        U0 = exp(-(X/Lambda).^2);               % Initial wave

        Unm = U0';
        Un = U0';
        Unp = 0*U0';
        NN = inv(M)*K; %#ok<*MINV>
        KK = sparse(zeros(size(K)));
        for oo=1:m
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
        
        for kt=2:N_t
            Unp = 2*Un-Unm+(2*KK)*Un;      % FD
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
        error_X(jx) = error_val;
    end
    error_vector(:,jt) = error_X;
    CFL(:,jt) = CFL_X;
end
                    toc
%% plotting: 3D

error_vector(error_vector>1) = 1;
error_vector(isnan(error_vector)) = 1;

figure(1)
surf(dT_L,h_L,error_vector), shading flat;
xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error"), zscale log;         xscale log, yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
view(2),%title("FD+LP");

% figure()
% surf(Nt_L,Nel_L,error_vector), shading flat;
% xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error"), zscale log;         xscale log, yscale log;
% colorbar, set(gca,'ColorScale','log'), clim([10^(-7) 10^(0)]);

%%
figure(2)
loglog(dT_L,error_vector(50,:),"*-"),hold on;
loglog(dT_L,error_vector(22,:),"*-");
loglog(dT_L,error_vector(10,:),"*-");
loglog(dT_L,error_vector(05,:),"*-");
xlim([10^(-4) 10^(-1)]),ylim([10^(-8) 10^(0)]);
xlabel("\Delta t"),ylabel("error"),grid on;

xlines = 100*[dT_L(60) dT_L(64)];
ylines = 100*[error_vector(22,60) error_vector(22,64)];
loglog([xlines(1) xlines(2)],[ylines(2) ylines(2)],"k",linewidth=1.5)
loglog([xlines(1) xlines(1)],[ylines(2) ylines(1)],"k",linewidth=1.5)
loglog(xlines,ylines,"k",linewidth=1.5)
slp1 = (log10(ylines(1))-log10(ylines(2)))/(log10(xlines(1))-log10(xlines(2)));
text(0.5*(xlines(1)+xlines(2)),ylines(2)/2,"1",FontSize=14)
text(round(100*xlines(1)+1)/100,(ylines(1)+3*ylines(2))/4,num2str(round(slp1)),FontSize=14)


legend(append("h=",num2str(h_L(50))), append("h=",num2str(h_L(22))),...
    append("h=",num2str(h_L(10))),append("h=",num2str(h_L(05))),"","","",Location="southeast");
fontsize(figure(2), 16, "points")
%%
figure(3)
loglog(h_L,error_vector(:,68),"*-"),hold on;
loglog(h_L,error_vector(:,52),"*-");
loglog(h_L,error_vector(:,25),"*-");
loglog(h_L,error_vector(:,05),"*-");
xlim([10^(-3) 10^(-1)]),ylim([10^(-8) 10^(0)]);                 xlim([4*10^(-4) 10^(-1)]);
xlabel("\Delta x"),ylabel("error"), grid on;

% xl = 2*[h_L(7) h_L(14)];
% yl = 2*[error_vector(7,68) error_vector(14,68)];

% for SEM
xl = 2*[h_L(4) h_L(8)];
yl = 2*[error_vector(4,68) error_vector(8,68)];

loglog([xl(1) xl(2)],[yl(2) yl(2)],"k",linewidth=1.5)
loglog([xl(1) xl(1)],[yl(2) yl(1)],"k",linewidth=1.5)
loglog([xl(1) xl(2)],yl,"k",linewidth=1.5)
slp2 = (log10(yl(1))-log10(yl(2)))/(log10(xl(1))-log10(xl(2)));
text(0.5*(xl(1)+xl(2)),yl(2)/2,"1",FontSize=14)
text(round(100*xl(1)+1)/100,(yl(1)+3*yl(2))/4,num2str((round(slp2*10)/10)),FontSize=14)

legend(append("\Delta t=",num2str(dT_L(68))), append("\Delta t=",num2str(dT_L(52))),...
    append("\Delta t=",num2str(dT_L(25))),append("\Delta t=",num2str(dT_L(05))),"","","",Location="west");
fontsize(figure(3), 16, "points")
%%
% error_vector = zeros(length(h_L),length(dT_L));
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


for ik=1:length(XBound)-11
    figure(1), hold on,
    plot3(dT_L(XBound(ik)),h_L(YBound(ik)),error_vector(YBound(ik),XBound(ik)),"*k");
    zscale log, xscale log, yscale log;
end

dt_f = dT_L(XBound);
dx_f = h_L(YBound);
alp = mean(dt_f(2:end)-dt_f(1:end-1))/mean(dx_f(2:end)-dx_f(1:end-1));
disp(alp)

figure(123),
pcolor(dT_L,h_L,error_vector), hold on;
xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error");
shading flat, xscale log,yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
plot(alp*dx_f,dx_f)
%%
% make lists with unique values
    % % xxx_dx = dx_f'; yyy_dt = dt_f';
    % % for lm=length(dx_f):-1:1
    % %     if sum(dx_f == dx_f(lm))>1
    % %         xxx_dx(lm) = [];
    % %         yyy_dt(lm) = [];
    % %     end
    % % end
    % % bb = xxx_dx\yyy_dt;
    % % disp(bb)
% figure(123), hold on;
% plot(bb*dx_f,dx_f)

% % figure(111)
% % surf(dT_L,h_L,CFL), shading flat;
% % xlabel("\Delta t"), ylabel("\Delta x"), zlabel("CFL"), 
% % zscale log, xscale log, yscale log;
% % colorbar, set(gca,'ColorScale','log')
% % XBound = zeros(length(dT_L),1); YBound = zeros(length(dT_L),1);
% % for kk=1:length(dT_L)
% %     idd = find(CFL(:,kk)>=0.99);
% %     if ~isempty(idd)
% %         YBound(kk) = idd(1);
% %         XBound(kk) = kk;
% %     end
% % end
% % XBound = nonzeros(XBound); YBound = nonzeros(YBound);
% % for ik=1:length(XBound)
% %     figure(111), hold on,
% %     plot3(dT_L(XBound(ik)),h_L(YBound(ik)),CFL(YBound(ik),XBound(ik)),"-*k");
% %     zscale log, xscale log, yscale log;
% % end
%%
unc = unique(dx_f);
dx_ff = zeros(length(unc),1);
dt_ff = zeros(length(unc),1);
for ik=length(unc):-1:1
    dx_ff(ik) = unc(ik);
    dt_ff(ik) = max(dt_f( (dx_f==unc(ik)) ) );
end
alp_hope = mean(dt_ff./dx_ff);
disp(alp_hope)

figure(456),
pcolor(dT_L,h_L,error_vector), hold on;
xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error");
shading flat, xscale log,yscale log;
colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
plot(alp_hope*dx_ff,dx_ff)
plot(dt_ff,dx_ff,"*k")

% [ZZ_dt,ZZ_h]=meshgrid(dT_L,h_L);
% ZZ_CFL = ZZ_dt./ZZ_h;
% ZZ_CFL((ZZ_CFL<1))=1;
% ZZ_CFL((ZZ_CFL<=2) & (ZZ_CFL>1.5) & (ZZ_CFL<1) )=1;
% ZZ_CFL((ZZ_CFL<=2) & (ZZ_CFL>1.5))=1;
% ZZ_CFL(ZZ_CFL>2)=0;
% ZZ_CFL((ZZ_CFL<=1.5) & (ZZ_CFL>1))=0;
% figure()
% surf(dT_L,h_L,ZZ_CFL),shading flat, xscale log,yscale log;
% colorbar, set(gca,'ColorScale','log'),clim([10^(-1) 10^(0)]);
% view (2)

% zzz = -ZZ_CFL+1;
% figure()
% surf(dT_L,h_L,error_vector.*zzz), shading flat;
% xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error"), zscale log;         xscale log, yscale log;
% colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
% view(2),
% figure()
% surf(dT_L,h_L,error_vector.*ZZ_CFL), shading flat;
% xlabel("\Delta t"), ylabel("\Delta x"), zlabel("error"), zscale log;         xscale log, yscale log;
% colorbar, set(gca,'ColorScale','log'), clim([10^(-6) 10^(0)]);
% view(2),