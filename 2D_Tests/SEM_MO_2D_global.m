clear variables
close all
clc
% loading the mesh data, the M & K matrices
load("MK_Full_Room_V7_2507.mat","Mesh","K","M","K_el","M_el")

% plotting the mesh
dd = plotting(Mesh,1,0,0);

% adjusting the matrices for the defined update equations
M = - M;
M_el = - M_el;
NN = inv(M)*K; %#ok<*MINV>

% Initial parameters
c = 1; L = 8; TM = 15;
Lambda = L/10;

m_L = 1:1:10;       % List of orders of the modified equation

xnd = Mesh.Nodes(:,1);      % defining x coords vector of the nodes 
ynd = Mesh.Nodes(:,2);      % defining y coords vector of the nodes 

% Source position
Xsrc = 2.5;
Ysrc = 3.5;
d = (xnd-Xsrc).^2+(ynd-Ysrc).^2;

% Initial condition
a = (1/pi/(Lambda^4));
b = (1-0.5*(d)./(Lambda^2));
U0 = a*b.*exp(-d./2./(Lambda^2));

% time steps and mesh size
Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
Nt_L = Nt_L(35:85);  
dT_L = TM./Nt_L;
h = min(dd);

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
        KK = sparse(zeros(size(K)));
        % calculating the modified equation matrix to the order m
        for oo=1:m
            KK = KK + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN)^(oo));
        end
        KK = sparse(KK);
        max_u = zeros(N_t,1);
        min_u = zeros(N_t,1);
        max_u(1) = max(U0);
        min_u(1) = min(U0);
        for kt=2:N_t
            Unp = 2*Un-Unm+( 2*KK )*Un;  % ((c*dT)^2)*NN
            Unm = Un;
            Un = Unp;
            max_u(kt) = max(Unp);
            min_u(kt) = min(Unp);
            if (abs(max(max_u(kt),min_u(kt)))==inf)||isnan(max_u(kt))||isnan(min_u(kt))||(abs(max(max_u(kt),min_u(kt)))>50)
                break
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
chk = -chk + 1;

%% max stable CFL per order
plt_dt = zeros(length(m_L),1);
for lk=1:length(m_L)
    tmp = chk(:,lk);
    idd = find(tmp==1);
    plt_dt(lk) = dT_L(idd(1));
end
cfl_f = c.*plt_dt./h;
figure(3),
plot(m_L,cfl_f,"o-",LineWidth=1.4,MarkerSize=10);
ylim([0 1.5*c*max(plt_dt)/h]);
grid on
xlabel("order"), ylabel("CFL");

%% function used locally only

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