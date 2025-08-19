clear variables
close all
clc

r = 3; % order of interp r+1
[x_gll,w_gll] = lglnodes(r); % positions and weights in [0;1]
n_gll = length(x_gll); % length of the corresponding vector

X = 0:0.001:1;

Y = zeros(r+1,length(X));
dY = zeros(r+1,length(X));

for ii=1:(r+1)
    Y(ii,:)=phi_hat(X,ii,r,x_gll);
end

figure(1),
plot(X,Y);grid on;hold on;
plot(x_gll,zeros(n_gll),"o")
legend("\Phi_1","\Phi_2","\Phi_3","\Phi_4")
title(append("r=",num2str(r)))

D = zeros(r+1,r+1);
for i=1:r+1
    for j=1:r+1
        D(i,j) =alternative_dl(i,x_gll,x_gll(j));
    end
end
figure(2), plot(x_gll,D',"-*")
legend("\partial\Phi_1","\partial\Phi_2","\partial\Phi_3","\partial\Phi_4")
disp("legned only for first 4 functions")

Ke = zeros(r+1);
for ii=1:r+1
    for jj=1:r+1
        for kk=1:r+1
            Ke(ii,jj)=Ke(ii,jj)+w_gll(kk)*D(ii,kk)*D(jj,kk);
        end
    end
end

% % fucnctions used

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


function [PP]=phi_hat(X,jj,r,x_gll) % if input=x_gll, maybe correct for phi_{p,j}

lp = [1:jj-1 jj+1:r+1];
P1 = 1; P2 = 1;
for pl=1:length(lp)
    P1 = P1.*(X-x_gll(lp(pl)));
end
for pk=1:length(lp)
    P2 = P2*(x_gll(jj)-x_gll(lp(pk)));
end
PP = P1./P2;
end

% % for derivatives
% function y=dl(i,x,z)
% n = length(x);
% y = 0;
% for m=1:n
%     if not(m==i)
%         y = y + 1/(x(m)-x(i));
%     end
% end
% % size(y)
% y = y*l(i,x,z);
% end
% function y=l(i,x,z)
% n = length(x);
% % computes h_i(z)
% y = 1;
% for m=1:n
%     if not(m==i)
%         y = y.*(z-x(m))./(x(i)-x(m));
%     end
% end
% end
% interp_nodes = @(X,p,x) (X(p+1)-X(p)).*x+X(p);  % F_p (adjusting the interval)
% nodes_weight = @(X,p,w_gll) (X(p+1)-X(p)).*w_gll;     % weights adjust with F_p
% chi_p = @(x,X,p) (X(p+1)>=x)*(x>=X(p));
% Lambda_p = @(x,X,p,r,x_gll) phi_hat(X,p-1,r+1,r,x_gll)*chi(x,X,p) + phi_hat(X,p+1,1,r,x_gll)*chi(x,X,p+1);
% Lambda_pj = @(x,X,p,jj,r,x_gll) phi_hat(X,p,jj,r+1,x_gll)*chi(x,X,p);