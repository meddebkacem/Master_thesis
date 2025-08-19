% clear variables
% close all
% clc

% the average run time is 20 to 30 minutes,depending on the chosen mesh,
% using the built-in parallel calculation of matlab with "parfar"

% uncomment the load line of the chosen mesh to test

% load("MK_Full_Room_V7_2507.mat","Mesh","K","M","K_el","M_el"); max_h_loc = 0.21; 
load("MK_Full_Room_test_I.mat","Mesh","K","M","K_el","M_el"); max_h_loc = 0.41;
% load("MK_Full_Room_test_II.mat","Mesh","K","M","K_el","M_el"); max_h_loc = 0.41;

cnctvty = Mesh.Connectivity(:,2:end);
M_el_C = zeros(size(M_el));

M = - M;
M_el = - M_el;

dd = plotting(Mesh,1,0,0,0);

% find the connecting/problematic elements

sep_id = find(dd<=max_h_loc);

NN = inv(M)*K; %#ok<*MINV>
% Initial parameters
c = 1; L = 8; TM = 35;
Lambda = L/10;

m_L = 1:1:10;               % List of orders of the modified equation

xnd = Mesh.Nodes(:,1);      % defining x coords vector of the nodes 
ynd = Mesh.Nodes(:,2);      % defining y coords vector of the nodes 

% Source position
Xsrc = 5;
Ysrc = 7;
d = (xnd-Xsrc).^2+(ynd-Ysrc).^2;

% Initial condition
a = (1/pi/(Lambda^4));
b = (1-0.5*(d)./(Lambda^2));
U0 = 2*a*b.*exp(-d./2./(Lambda^2));

% time steps and mesh size
Nt_L = [51:2:99 101:10:291 301:20:1001 1101:50:2001];
% Nt_L = Nt_L(20:67);       % uncomment to reduce the number of iterations
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
        KK = zeros(size(K));

        % matrix with the selected element + calculation of modified 
        % equation matrix to order m locally from m = 2
        NN_Mo = NN(unique(cnctvty(sep_id,:)),unique(cnctvty(sep_id,:)));
        KK_Mo = zeros(size(NN_Mo));
        for oo=2:m
            KK_Mo = KK_Mo + ( (c*dT)^(2*oo)/factorial(2*oo) )*((NN_Mo)^(oo));
        end
        % adding the obtained matrix to the global one that
        % represents the first order of calculation (m=1) 
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