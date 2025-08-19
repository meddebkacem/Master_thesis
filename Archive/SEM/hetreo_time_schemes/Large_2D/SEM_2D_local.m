clear variables
close all
clc

% load("MK_Full_Room1_r=4.mat","Mesh","K","M","K_el","M_el");

% load("MK_Full_Room_V7_2507.mat","Mesh","K","M","K_el","M_el")    % default for best results "MK_Full_Room_V2.mat"

    load("MK_Full_Room_test_II.mat","Mesh","K","M","K_el","M_el")
    % Room_test_I >> 2.16 improvement in dt starting 4th power of MO 
    % using 29 elts over 480 elts total, runtime 1883.31 seconds

    % Room_test_II >> 2.29 improvement in dt starting 4th power of MO 
    % using 20 elts over 418 elts total, runtime 1399 seconds

cnctvty = Mesh.Connectivity(:,2:end);
M_el_C = zeros(size(M_el));

M = - M;
M_el = - M_el;

%% finding connectivity per node
% N_elts = size(cnctvty,1);
% N_per_elt = size(cnctvty,2);
% for lk=1:N_elts
%     M_el_temp = zeros(size(M_el(:,:,1)));
%     ord = cnctvty(lk,:);
%     for jk=1:length(ord)
%         [x_ind,y_ind] = find(cnctvty==ord(jk));
%         sorted_ord = sort(ord);
%         indd = sum(ord(jk)>=sorted_ord);
%         for nn=1:length(x_ind)
%             M_el_temp(indd,indd) = M_el_temp(indd,indd)+M_el(y_ind(nn),y_ind(nn),x_ind(nn));
%         end
%     end
%     M_el_C(:,:,lk) = M_el_temp;
% end
% M_el_C = - M_el_C;
% %% assembly V2.0
% MMM = zeros(length(unique(cnctvty)));
% KKK = zeros(length(unique(cnctvty)));
% for ll=1:N_elts
%     MMM(cnctvty(ll,:),cnctvty(ll,:)) = M_el_C(:,:,ll);
%     KKK(cnctvty(ll,:),cnctvty(ll,:)) = K_el(:,:,ll);
% end
% MMM = sparse(MMM);
% KKK = sparse(KKK);

%%
dd = plotting(Mesh,1,0,0,0);

kk = 0.1:0.05:3;
class_of_elts = zeros(length(kk),1);
for kkl=1:length(kk)
    class_of_elts(kkl) = sum(dd<=kk(kkl));
end
figure,
plot(kk,class_of_elts,"x-k",MarkerSize=10);

% find the connecting/problematic element
% sep_id = find(dd<=0.21);
sep_id = find(dd<=0.41);   % 0.3 V2.mat et 0.21 pour la V3.matc
        % sep_id = [310 311 320 321 340 341 344:355 364 365 372 373 391 392]

NN = inv(M)*K; %#ok<*MINV>

c = 1; L = 8; TM = 15+20*0;
Lambda = L/10;

m_L = 1:1:10;

xnd = Mesh.Nodes(:,1);
ynd = Mesh.Nodes(:,2);

Xsrc = 5;
Ysrc = 7;
d = (xnd-Xsrc).^2+(ynd-Ysrc).^2;

a = (1/pi/(Lambda^4));
b = (1-0.5*(d)./(Lambda^2));
U0 = 2*a*b.*exp(-d./2./(Lambda^2));

Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
        % Nt_L = Nt_L(20:67); % reduce iterations
dT_L = TM./Nt_L;
h = min(dd);

% rect1 = find(xnd<=20)';
% rect2 = find((xnd>=20)&(xnd<=20.8)&(ynd>=8)&(ynd<=8.8))';
% rect3 = find(xnd>=20.8)';
% xn1 = xnd(rect1); xn2 = xnd(rect2); xn3 = xnd(rect3);
% yn1 = ynd(rect1); yn2 = ynd(rect2); yn3 = ynd(rect3);
            % % rect1 = find(xnd<=10)';
            % % rect2 = find((xnd>=10)&(xnd<=10.4)&(ynd>=4)&(ynd<=4.4))';
            % % rect3 = find(xnd>=10.4)';
            % % xn1 = xnd(rect1); xn2 = xnd(rect2); xn3 = xnd(rect3);
            % % yn1 = ynd(rect1); yn2 = ynd(rect2); yn3 = ynd(rect3);
% tri1 = delaunay(xn1,yn1);
% tri2 = delaunay(xn2,yn2);
% tri3 = delaunay(xn3,yn3);
% figure,
% trisurf(tri1, xn1, yn1, U0(rect1)), hold on;
% trisurf(tri2, xn2, yn2, U0(rect2));
% trisurf(tri3, xn3, yn3, U0(rect3)); hold off;
% shading interp;
% xlim([-0.1 37]);
% ylim([-0.1 24.2]);
% zlim([-2 2]);
%             % xlim([-0.1 18.5]);
%             % ylim([-0.1 12.1]);
%             % zlim([-2 2])
% 
% grid on;
% % view(2);
% colorbar,
% % clim([-0.5 1]);
% clim([-0.25 0.25]);

            tic
CFL = zeros(1,length(dT_L));
chk = zeros(length(Nt_L),length(m_L));
for jt=1:length(Nt_L)
    dT = dT_L(jt);
    N_t = Nt_L(jt);
    CFL(jt) = c*dT/h;
    chk_m = zeros(length(m_L),1);
    parfor jm=1:length(m_L)
        m = m_L(jm);
        Unm = U0;
        Un = U0;
        KK = zeros(size(K));
                    NN_Mo = NN(unique(cnctvty(sep_id,:)),unique(cnctvty(sep_id,:)));
                    KK_Mo = zeros(size(NN_Mo));
                    for oo=2:m
                        KK_Mo = KK_Mo + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN_Mo)^(oo));
                    end
                    KK(unique(cnctvty(sep_id,:)),unique(cnctvty(sep_id,:))) = KK_Mo;
        KK = sparse(KK);
        max_u = zeros(N_t,1);
        min_u = zeros(N_t,1);
        max_u(1) = max(U0);
        min_u(1) = min(U0);
        % figure
        for kt=2:N_t
            Unp = 2*Un-Unm+( ((c*dT)^2)*NN + 2*KK )*Un;
            Unm = Un;
            Un = Unp;
            max_u(kt) = max(Unp);
            min_u(kt) = min(Unp);
            if (abs(max(max_u(kt),min_u(kt)))==inf)||isnan(max_u(kt))||isnan(min_u(kt))||(abs(max(max_u(kt),min_u(kt)))>50)
                break
            end

           % if mod(kt,10)==0
           %      tri1 = delaunay(xn1,yn1);
           %      tri2 = delaunay(xn2,yn2);
           %      tri3 = delaunay(xn3,yn3);
           %      trisurf(tri1, xn1, yn1, Unp(rect1)), hold on;
           %      trisurf(tri2, xn2, yn2, Unp(rect2));
           %      trisurf(tri3, xn3, yn3, Unp(rect3)); hold off;
           %      shading interp;
           %      xlim([-0.1 37]);
           %      ylim([-0.1 24.2]);
           %      zlim([-2 2]);
           %              % xlim([-0.1 18.5]);
           %              % ylim([-0.1 12.1]);
           %              % zlim([-2 2])
           %      grid on;
           %      % view(2);
           %      colorbar,
           %      % clim([-0.5 1]);
           %      clim([-0.25 0.25]);
           %      % fontsize(gcf, 14, "points")
           %      drawnow;
           %  end

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
            toc
chk(chk==inf)=1;
chk(isnan(chk))=1;
chk = -chk+1;

plt_dt = zeros(length(m_L),1);
for lk=1:length(m_L)
    tmp = chk(:,lk);
    idd = find(tmp==1);
    plt_dt(lk) = dT_L(idd(1));
end
cfl_f = c.*plt_dt./h;
%% max stable CFL per order
figure(21),
plot(m_L,cfl_f,"-*",LineWidth=1.4,MarkerSize=10); grid on;
ylim([0 0.4]);
xlabel("order"), ylabel("CFL");
% fontsize(figure(21), 16, "points")
%% making plot of the domain 18.4 by 12 (the small domain)
x1 = [0 10 10 10.4 10.4 18.4 18.4 10.4 10.4 10 10 0 0];
y1 = [0 0 4 4 0 0 12 12 4.4 4.4 12 12 0];
figure, plot(x1,y1,"-k"),grid on, axis equal
xlim([-0.1 18.5])
xlim([-0.5 18.9])

%% making plot of the domain 36.8 by 24 (the larger domain)
x2 = [0 20 20 20.8 20.8 36.8 36.8 20.8 20.8 20 20 0 0];
y2 =     [0 0 8 8 0 0 24 24 8.8 8.8 24 24 0];
figure, plot(x2,y2,"-k"),grid on, axis equal
xlim([-1 38])


%% functions used locally
function dd = plotting(Mesh,plt_lines,plt_nodes,txt_plt,plt_close)
dd = zeros(size(Mesh.Connectivity,1),1);
if plt_lines
f = figure;
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
if plt_close
    close(f)
end
end