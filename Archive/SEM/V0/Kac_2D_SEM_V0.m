clear variables
close all
clc

% load("MK_Full_Room.mat");
% load("MK_trunc_Room.mat");
load("MK_thin_tube_Room.mat");

M = -M;

c = 1; L = 8; TM = 15; % 0.6
Lambda = L/10;

m_L = 1; % 1:1:20;

xNd = Mesh.Nodes(:,1);
yNd = Mesh.Nodes(:,2);

Xsrc = 2.5;
Ysrc = 5.5;
d = (xNd-Xsrc).^2+(yNd-Ysrc).^2;

a = (1/pi/(Lambda^4));
b = (1-0.5*(d)./(Lambda^2));
U0 = a*b.*exp(-d./2./(Lambda^2));

Nt_L = 301;
dT_L = TM./Nt_L;
h = 0.5;

[XX,YY] = meshgrid(unique(xNd),unique(yNd));
% % for full room
% vizz = reshape(U0,[41,61]);
% figure,surf(XX,YY,vizz),shading flat

% % for truncated room

rect1 = find(xNd<=5)';

% x_L = length(unique(xNd));
y_L = length(unique(yNd(rect1)));

rect1resh = [y_L,length(rect1)/y_L];

rect12 = find((xNd==5)&(yNd>=4)&(yNd<=5))';
rect23 = find((xNd==8)&(yNd>=4)&(yNd<=5))';
rect2 = find((xNd>5)&(xNd<8))';
rect2resh = [5,2+length(rect2)/5];

rect3 = find(xNd>=8)';
rect3resh = [y_L,length(rect3)/y_L];

xrect1 = unique(xNd(rect1));
xrect2 = reshape(xNd([rect12 rect2 rect23]),[5 13]);
xrect3 = unique(xNd(rect3));

yrect1 = unique(yNd(rect1));
yrect2 = reshape(yNd([rect12 rect2 rect23]),[5 13]);
yrect3 = unique(yNd(rect3));


CFL = zeros(1,length(dT_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/h;
                            % Mo(N_t) = struct("cdata",[],"colormap",[]);
    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0;
        Un = U0;
        NN = inv(M)*K; %#ok<*MINV>
        KK = sparse(zeros(size(K)));
        for oo=1:m
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
                            % f = figure;
                            % f.Visible = "off";
                            % Unp_vis = reshape(U0,[41,61]);
                            % surf(XX,YY,Unp_vis);
                            % shading flat;
                            % ylim([-0.1 10.1]);
                            % xlim([-0.1 15.1]);
                            % zlim([-0.5 1.5])
                            % grid on;
                            % colorbar
                            % clim([-0.3 1]);
                            % % fontsize(gcf, 14, "points")
                            % drawnow;
                            % Mo(1) = getframe(gcf);
        for kt=2:N_t
            Unp = 2*Un-Unm+(2*KK)*Un;
            Unm = Un;
            Un = Unp;
            if mod(kt,1)==0
        % Unp_vis = reshape(Unp,[41,61]);   % for full room
        % surf(XX,YY,Unp_vis);              % for full room
                              
                Unp1 = reshape(Unp(rect1),rect1resh);
                Unp2 = reshape(Unp([rect12 rect2 rect23]),rect2resh);
                Unp3 = reshape(Unp(rect3),rect3resh);

                surf(xrect1,yrect1,Unp1), hold on;
                surf(xrect2,yrect2,Unp2)
                surf(xrect3,yrect3,Unp3), hold off;
                shading flat;
                ylim([-0.1 10.1]);
                xlim([-0.1 15.1]);
                zlim([-0.5 1.5])
                grid on;
                colorbar
                % view(2)
                clim([-0.3 1]);
                % fontsize(gcf, 14, "points")
                drawnow;
                pause(0.05)

                            % Mo(kt) = getframe(gcf);
            end
        end
    end
end

                            % f.Visible = "on";
                            % movie(Mo)
                            % wvidObj = VideoWriter("plot2D","MPEG-4");
                            % open(wvidObj)
                            % writeVideo(wvidObj,Mo)
                            % close(wvidObj)


%% transformation function

% x1 = 5:0.001:6;
% y1low = x1/4+11/4;
% y1high = -x1/4 + 25/4;
% 
% x3 = 7:0.001:8;
% y3low = -x3/4+6;
% y3high = x3/4+3;
% 
% x2 = 6:0.001:7;
% y2low = 0*x+4;
% y2high = 0*x+5;
% yc = 4.5;
% y2low = yc + (yc-y2low)/2;
% y2high = yc - (y2high-yc)/2;
% 
% 
% figure, plot(x1,y1low,x1,y1high)
% hold on
% plot(x3,y3low,x3,y3high)
% plot(x2,y2low,x2,y2high)

%% Boundary condititions

% xx1 = find(xNd==0);
% xx2 = find(xNd==15);
% figure,plot(xNd(xx1), yNd(xx1),"*b", xNd(xx2), yNd(xx2),"*b");
% 
% yy1 = find(yNd==0);
% yy2 = find(yNd==10);
% hold on, plot(xNd(yy1), yNd(yy1),"*b", xNd(yy2), yNd(yy2),"*b");
% 
% yy11 = find((xNd>=5)&(xNd<=8)&(yNd==4));
% yy12 = find((xNd>=5)&(xNd<=8)&(yNd==5));
% hold on, plot(xNd(yy11), yNd(yy11),"*b", xNd(yy12), yNd(yy12),"*b");
% 
% xx111 = find((yNd<=4)&(xNd==5));
% xx112 = find((yNd>=5)&(xNd==5));
% hold on, plot(xNd(xx111), yNd(xx111),"*b", xNd(xx112), yNd(xx112),"*b");
% 
% xx121 = find((yNd<=4)&(xNd==8));
% xx122 = find((yNd>=5)&(xNd==8));
% hold on, plot(xNd(xx121), yNd(xx121),"*b", xNd(xx122), yNd(xx122),"*b");
% % % % % % figure,
% % % % % % plot3(xNd,yNd,rkr,".")
% % % % % % shading flat