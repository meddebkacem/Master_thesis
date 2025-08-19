function [M,K,X] = SEM_MK(r,h,L)

[x_gll,w_gll] = lglnodes(r);
n_gll = length(x_gll);
Me = diag(w_gll);
D = zeros(r+1,r+1);     % elementary derivatives with gll
for iii=1:r+1
    for jjj=1:r+1
        D(iii,jjj) =alternative_dl(iii,x_gll,x_gll(jjj));
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
N0 = round((2*L*r/h)+1);
X = linspace(-L,L,N0);
M = zeros(N0); K = zeros(N0);
for il=1:(n_gll-1):length(X)-n_gll+1
    X(il:il+n_gll-1)=x_gll*(X(il+n_gll-1)-X(il))+X(il);
    M(il:il+n_gll-1,il:il+n_gll-1) = M(il:il+n_gll-1,il:il+n_gll-1) + h*Me;
    K(il:il+n_gll-1,il:il+n_gll-1) = K(il:il+n_gll-1,il:il+n_gll-1) + Ke;
end
K = K./h;
M = -M;

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
end