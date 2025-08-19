clear variables
close all
clc

load("MK_Room_el_M_el_K.mat");
dd = plotting(Mesh,1,0,0);

% find the connecting element
sep_id = find(dd<=0.3);
for ki=size(sep_id):-1:1
    if (any(Mesh.Nodes(Mesh.Connectivity(sep_id(ki),2:5),1)<=4.2))
% &&any(Mesh.Nodes(Mesh.Connectivity(ki,2:5),1)==5.1))&&(any(Mesh.Nodes(Mesh.Connectivity(ki,2:5),2)==2)
        sep_id(ki) = []; % Mesh.Connectivity(ki,:);
    end
end


M = -M;
M_el = -M_el;

c = 1; L = 8; TM = 15; % 0.6
Lambda = L/10;

m_L = 2;%1:1:20;

xnd = Mesh.Nodes(:,1);
ynd = Mesh.Nodes(:,2);

Xsrc = 2.5;
Ysrc = 3.5;
d = (xnd-Xsrc).^2+(ynd-Ysrc).^2;

a = (1/pi/(Lambda^4));
b = (1-0.5*(d)./(Lambda^2));
U0 = a*b.*exp(-d./2./(Lambda^2));
% triii = delaunay(xnd,ynd);figure,trisurf(triii,xnd,ynd,abs(fftshift(fft(U0))));
Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
        % Nt_L = Nt_L(35:85);
        Nt_L = Nt_L(51);
dT_L = TM./Nt_L;
h = 0.1;

% % for truncated room
rect1 = find(xnd<=5)';
rect2 = find((xnd>=5)&(xnd<=5.1)&(ynd>=2)&(ynd<=2.2))';
rect3 = find(xnd>=5.1)';
xn1 = xnd(rect1); xn2 = xnd(rect2); xn3 = xnd(rect3);
yn1 = ynd(rect1); yn2 = ynd(rect2); yn3 = ynd(rect3);

% prepare the element matrix of the separation
NN_el = zeros(size(M));

for kl=1:length(sep_id)
    indx = Mesh.Connectivity(sep_id(kl),2:end);
    NN_el(indx,indx) = inv(M_el(:,:,sep_id(kl)))*K_el(:,:,sep_id(kl));
end

CFL = zeros(1,length(dT_L));
chk = zeros(length(Nt_L),length(m_L));
for jt=1:length(Nt_L)          % you can add "par"
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/h;
                            % Mo(N_t) = struct("cdata",[],"colormap",[]);
    chk_m = zeros(length(m_L),1);
    for jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0;
        Un = U0;
        NN = inv(M)*K; %#ok<*MINV>
        KK = zeros(size(K));
        for oo=2:m
            KK = KK+ ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN_el)^(oo));
        end
            % TuT = ((c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo)); % 1st elt modified equation
            % TuT1 = zeros(size(M));
            % TuT1 = zeros(size(M));
            % for kl=1:length(sep_id)
            %     indx = Mesh.Connectivity(sep_id(kl),2:end);
            %     TuT1(indx,indx) = TuT(indx,indx) ;
            % end
        KK = sparse(KK);
            % f = figure;
            % f.Visible = "off";
            % tri1 = delaunay(xn1,yn1);
            % tri2 = delaunay(xn2,yn2);
            % tri3 = delaunay(xn3,yn3);
            % trisurf(tri1, xn1, yn1, Unp(rect1)), hold on;
            % trisurf(tri2, xn2, yn2, Unp(rect2));
            % trisurf(tri3, xn3, yn3, Unp(rect3)); hold off;
            % shading interp;
            % xlim([-0.1 9.3]);
            % ylim([-0.1 6.1]);
            % zlim([-1.1 1.1])
            % grid on;
            % view(2);
            % colorbar,
            % clim([-1.0 1.5]);
            % fontsize(gcf, 14, "points")
            % drawnow;
            % Mo(1) = getframe(gcf);
        max_u = zeros(N_t,1);
        min_u = zeros(N_t,1);
        max_u(1) = max(U0);
        min_u(1) = min(U0);
        figure
        for kt=2:N_t
            Unp = 2*Un-Unm+( ((c*dT)^2)*NN + 2*KK)*Un;
            Unm = Un;
            Un = Unp;
            max_u(kt) = max(Unp);
            min_u(kt) = min(Unp);
            
            % if (abs(max(max_u(kt),min_u(kt)))==inf)||isnan(max_u(kt))||isnan(min_u(kt))||(abs(max(max_u(kt),min_u(kt)))>50)
            %     break
            % end

            if mod(kt,1)==0
                tri1 = delaunay(xn1,yn1);
                tri2 = delaunay(xn2,yn2);
                tri3 = delaunay(xn3,yn3);
                trisurf(tri1, xn1, yn1, Unp(rect1)), hold on;
                trisurf(tri2, xn2, yn2, Unp(rect2));
                trisurf(tri3, xn3, yn3, Unp(rect3)); hold off;
                shading interp;
                xlim([-0.1 9.3]);
                ylim([-0.1 6.1]);
                zlim([-1.1 1.1])
                grid on;
                % view(2);
                colorbar,
                clim([-1.0 1.5]);
                % fontsize(gcf, 14, "points")
                drawnow;
                            % Mo(kt) = getframe(gcf);
            end
        end
        if any(max_u>=10)||any(min_u<=-10)
            chk_m(jm) = inf;
        end
        if any(max_u==inf)||any(min_u==-inf)
            chk_m(jm) = inf;
        end
        if any(isnan(max_u))||any(isnan(min_u))
            chk_m(jm) = NaN;
        end
    end
    chk(jt,:) = chk_m;
end
chk(chk==inf)=1;
chk(isnan(chk))=1;
chk = -chk+1;
dt_M = repmat(dT_L,[length(m_L),1])';

figure,
surf(m_L,dT_L,dt_M.*chk)
shading flat, yscale log;

figure,plot(m_L,dT_L(end-sum(chk)+1))

                            % f.Visible = "on";
                            % movie(Mo)
                            % wvidObj = VideoWriter("plot2D","MPEG-4");
                            % open(wvidObj)
                            % writeVideo(wvidObj,Mo)
                            % close(wvidObj)
%%
plt_dt = zeros(length(m_L),1);
for lk=1:length(m_L)
    tmp = chk(:,lk);
    idd = find(tmp==1);
    plt_dt(lk) = dT_L(idd(1));
end
figure,
plot(m_L,c.*plt_dt./h,"o-");
ylim([0 1]);
disp("max delta t:")
disp(2/max(eigs(-M,K,6,'largestabs')))

%% Boundary condititions

function dd = plotting(Mesh,plt_lines,plt_nodes,txt_plt)
dd = zeros(size(Mesh.Connectivity,1),1);
if plt_lines
figure;
hold on;
axis equal;
% Plot elements
for el = 1:size(Mesh.Connectivity,1)
nodeIDs = Mesh.Connectivity(el,2:5);
x_coords = Mesh.Nodes(nodeIDs,1);
y_coords = Mesh.Nodes(nodeIDs,2);
% Close the quadrilateral
x_coords(5) = x_coords(1);
y_coords(5) = y_coords(1);
dist = (x_coords(2:5)-x_coords(1:4)).^2 + (y_coords(2:5)-y_coords(1:4)).^2 ;
dd(el) = min(sqrt(dist));
plot(x_coords, y_coords, 'b-');
if txt_plt
    text(mean(x_coords),mean(y_coords),num2str(el),"FontSize",14)
end
end
end
if plt_nodes
    plot(Mesh.Nodes(:,1), Mesh.Nodes(:,2), 'ro', 'MarkerFaceColor', 'r');
end
end

%% for viz of matrice parts

% tt = zeros(size(K));
% tt(sep_id(2:end),sep_id(2:end))=1;
% K1 = K;
% K1(Mesh.Connectivity(sep_id(1),2:end),Mesh.Connectivity(sep_id(1),2:end))=0;
% spy(K1,"r")
% hold on, spy(tt)