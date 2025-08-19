function [M,K,X,hh,N1,Me_S,Ke_S] = SEM_MK_hetero(r,h,L)

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

hh = h.*ones(N0,1);
N1 = round(0.3*N0)+1;
hh(N1:N1+r-1) = h/10;


X = linspace(-L,L,N0);
M = zeros(N0); K = zeros(N0);
for il=1:(r):length(X)-n_gll+1
    % X(il:il+r)=x_gll*(X(il+r)-X(il))+X(il);
            X(il:il+r)=x_gll.*hh(il)+X(il);
    M(il:il+r,il:il+r) = M(il:il+r,il:il+r) + hh(il)*Me;
    K(il:il+r,il:il+r) = K(il:il+r,il:il+r) + Ke/hh(il);
    if il == N1
        Me_S = hh(il)*Me;
        Ke_S = Ke/hh(il);
    end
end
M = -M;
Me_S = -Me_S;

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