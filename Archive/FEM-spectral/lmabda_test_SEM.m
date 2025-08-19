clear variables
close all
clc

r = 3;

[x_gll,w_gll] = lglnodes(r);        % positions and weights in [0;1]
Me = diag(w_gll);

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
N = inv(Me)*Ke; %#ok<*MINV> PDF-186 eq.11.16

LL = max(eig(N));

cfl_P2_P3 = (2*sqrt(3)/sqrt(LL));
disp(append("r=",num2str(r),", cfl P₂_2 P₃_3: "))
disp(cfl_P2_P3)
disp(append("r=",num2str(r),", cfl P₄_4"))
cfl_P4 = sqrt(2*(5^(3/2)-25^(1/3)+5))/sqrt(LL);
disp(cfl_P4)


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