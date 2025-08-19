clear variables
close all
clc

c = 1; L = 1; TM = 0.6;
h = 0.02;
Lambda = L/10;
dT = 0.001;
T = 0:dT:TM;
N0=round(2*L/h+1);
r = 3;                              % order of interp r-1
[x_gll,w_gll] = lglnodes(r);        % positions and weights in [0;1]
n_gll = length(x_gll);              % length of the corresponding vector


X = linspace(-L,L,(N0-1)*(n_gll-2)+N0);     % position (interp pts included)
M = zeros(1,length(X));                     % Mass matrix diagonal

for ij=1:(n_gll-1):length(X)-n_gll+1
    X(ij:ij+n_gll-1)=x_gll*(X(ij+n_gll-1)-X(ij))+X(ij);     % positions for plotting purpose
    M(ij:ij+n_gll-1)=M(ij:ij+n_gll-1) + (X(ij+n_gll-1)-X(ij))*w_gll';
end
M = sparse(diag(M));    % making sparse mass matrix
D = zeros(r+1,r+1);     % elementary derivatives with gauss lobatto interp
for i=1:r+1
    for j=1:r+1
        D(i,j) =alternative_dl(i,x_gll,x_gll(j));
    end
end
Ke = zeros(r+1);        % elementary stiffness matrix
for ii=1:r+1
    for jj=1:r+1
        for kk=1:r+1
            Ke(ii,jj)=Ke(ii,jj)+w_gll(kk)*D(ii,kk)*D(jj,kk);
        end
    end
end
Me=diag(w_gll);


K = zeros(length(X));   % global stiffness matrix (matching the Mass matrix)

for il=1:(n_gll-1):length(X)-n_gll+1
    K(il:il+n_gll-1,il:il+n_gll-1) = K(il:il+n_gll-1,il:il+n_gll-1) + Ke;
end
K = K./h;


U0 = exp(-(X/Lambda).^2);               % Initial wave
V0 = zeros(length(X),1);                % Initial velocity
UU = zeros(length(U0),length(T));       % saving all time instances
UU(:,1) = U0;
Vn_1 = V0 - (dT/2)*(-(c^2)*M\K)*U0';    % Initializing V_{-1/2}


for kt=2:length(T)
    An = (-(c^2)*M\K)*UU(:,kt-1);
    Vn_2 = Vn_1+dT*An;
    UU(:,kt) = UU(:,kt-1)+dT*Vn_2;
    Vn_1 = Vn_2;
end

% final theoretical solution
theo_fn = 0.5*exp(-((X-c*TM)/Lambda).^2) + 0.5*exp(-((X+c*TM)/Lambda).^2);

% % % Plots
figure(11)
for tt=1:20:length(T) 
    plot(X,theo_fn,"--"),hold on
    plot(X,UU(:,tt)),hold off;
    ylim([-1.2 1.2]);
    xlabel("Position (m)"), ylabel("amplitude")
    pause(0.05);
    drawnow
end

figure()
plot(X,theo_fn,"-",X,UU(:,end),"-.",LineWidth=2), grid on;
xlabel("Position (m)"), ylabel("amplitude")
legend("theory","calc"),ylim([-.5 1.2]);
title("Wave propagation");

% % % function: interpolation derivatives for gll
function y = alternative_dl(j,x,z)
y = 0;
n = length(x);
for l=1:n
    if not(l==j)
        k = 1/(x(j)-x(l));
        for m=1:n
            if not(m==j) && not(m==l)
                k = k*(z-x(m))/(x(j)-x(m));
            end
        end
        y = y + k;
    end
end
end