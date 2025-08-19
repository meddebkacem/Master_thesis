function [M,K,X,hh,N1] = FD_MK_hetero(or,h,L)

N0 = round((2*L/h)+1);
X = -L:h:L;
M = eye(N0);
K = zeros(N0);

hh = h.*ones(N0,1);
N1 = round(0.3*N0);
hh(round(0.3*N0)) = h/10;

alpha_pm = @(m,p) ((-1)^(p-1))*(2*factorial(m)^2)/(factorial(m-p))/(factorial(m+p));

for p=1:or
    Kph = (-2*diag(ones(N0,1))+diag(ones(N0-p,1),p)+diag(ones(N0-p,1),-p));
    Kph(1,N0-p+1) = 1;
    Kph(N0,p) = 1;
    Kph = sparse( Kph./ ((p.*hh).^2) );
    K = K + alpha_pm(or,p)*Kph;
end

end