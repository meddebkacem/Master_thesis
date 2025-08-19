clear variables
close all
clc

c = 1; L = 1; TM = 0.6;%*(5/4);
Lambda = L/10;
r = 3;                                  % order of interp r-1 for SEM
m_L = 1:1:2 ;                           % order of taylor approx

Nt_L = [51:5:201 201:10:2001 2201:100:5001];
                                        % Nt_L = Nt_L([65:68 76:93 117:126]);
Nel_L = 2*100;
h = 2*L./Nel_L;
dT_L = TM./Nt_L;

[M,K,X,hh,Me_L1,Ke_L,N1] = SEM_MK_hetero_elt(r,h,L);
AAA = inv(M)*K;

M = sparse(M); K = sparse(K);

Me_L = zeros(size(Me_L1));

for kkl=2:size(Me_L1,3)-1
    Me_L(:,:,kkl) = Me_L1(:,:,kkl) ;
    Me_L(1,1,kkl) = Me_L1(1,1,kkl) + Me_L1(end,end,kkl-1);
    Me_L(end,end,kkl) = Me_L1(end,end,kkl) + Me_L1(1,1,kkl+1);
end
Me_L(:,:,1) = Me_L1(:,:,1);
Me_L(end,end,1) = Me_L1(end,end,1) + Me_L1(1,1,2);

Me_L(:,:,end) = Me_L1(:,:,end);
Me_L(1,1,end) = Me_L1(1,1,end) + Me_L1(end,end,end-1);


[K_part_sq,K_sq_bloc] = square_parts_KM(Me_L,Ke_L,r,Nel_L);
KK_test = (K_part_sq + K_sq_bloc) - AAA^2;
%% testing
% 
% N = 3*r+1;
% KK_assemb_glob = zeros(N);
% KK_assemb_sq = zeros(N);
% KK_part = zeros(N);
% 
% A = inv(Me_L(:,:,1))*Ke_L(:,:,1);
% B = inv(Me_L(:,:,2))*Ke_L(:,:,2);
% C = inv(Me_L(:,:,end))*Ke_L(:,:,end);
% 
% KK_assemb_glob(1:4,1:4) = A;
% KK_assemb_glob(4:7,4:7) = KK_assemb_glob(4:7,4:7) + B;
% KK_assemb_glob(7:10,7:10) = KK_assemb_glob(7:10,7:10) + C;
% 
% KK_assemb_sq(1:4,1:4) = A^2;
% KK_assemb_sq(4:7,4:7) = KK_assemb_sq(4:7,4:7) + B^2;
% KK_assemb_sq(7:10,7:10) = KK_assemb_sq(7:10,7:10) + C^2;
% 
% K12= A(1:end-1,4);
% K21= A(4,1:end-1);
% K2 = A(end,end);
% 
% K3 = B(1,1);
% K34= B(1,2:end);
% K43= B(2:end,1);
% 
% KK_part(1:4,4:7) = KK_part(1:4,4:7) + [K12;K2]*[K3 K34];
% KK_part(4:7,1:4) = KK_part(4:7,1:4) + [K3;K43]*[K21 K2];
% 
% clear K12 K21 K2 K3 K34 K43
% 
% K12= B(1:end-1,4);
% K21= B(4,1:end-1);
% K2 = B(end,end);
% K3 = C(1,1);
% K34= C(1,2:end);
% K43= C(2:end,1);
% KK_part(4:7,7:10) = KK_part(4:7,7:10) + [K12;K2]*[K3 K34];
% KK_part(7:10,4:7) = KK_part(7:10,4:7) + [K3;K43]*[K21 K2];
% 
% clear K12 K21 K2 K3 K34 K43 temp1 temp2
% 
% 
% KK_test = KK_assemb_sq + KK_part;
% KK_square_ref = KK_assemb_glob^2;
% Diff = KK_square_ref - KK_test;
% Diff(abs(Diff)<10^(-4))=0;
%%
function [Pw_N] = Nth_power(KK_pow_2,KK_pow_1,n_pow)

switch mod(n_pow,2)
    case 0
        M = KK_pow_2^(floor(n_pow/2));
    case 1
        M = KK_pow_2^(floor(n_pow/2))*KK_pow_1;
end
Pw_N = M;
end



%%


function [KK_part,KK_sq] = square_parts_KM(Me_L,Ke_L,r,N_elts)

N = N_elts*r+1;
KK_part = zeros(N);
KK_sq = zeros(N);

A = inv(Me_L(:,:,1))*Ke_L(:,:,1);
KK_sq(1:r+1,1:r+1) = A^2;

for jj=2:N_elts
    indx = (jj-1)*r+1;
    A = inv(Me_L(:,:,jj-1))*Ke_L(:,:,jj-1);
    B = inv(Me_L(:,:,jj))*Ke_L(:,:,jj);
    KK_sq(indx:indx+r,indx:indx+r) = KK_sq(indx:indx+r,indx:indx+r) + B^2;

    K12= A(1:end-1,4);
    K21= A(4,1:end-1);
    K2 = A(end,end);
    
    K3 = B(1,1);
    K34= B(1,2:end);
    K43= B(2:end,1);
    
    KK_part(indx-r:indx,indx:indx+r) = KK_part(indx-r:indx,indx:indx+r) + [K12;K2]*[K3 K34];
    KK_part(indx:indx+r,indx-r:indx) = KK_part(indx:indx+r,indx-r:indx) + [K3;K43]*[K21 K2];

    % KK_part(1:4,4:7) = KK_part(1:4,4:7) + [K12;K2]*[K3 K34];
    % KK_part(4:7,1:4) = KK_part(4:7,1:4) + [K3;K43]*[K21 K2];
    clear K12 K21 K2 K3 K34 K43
end
KK_part = sparse(KK_part);
KK_sq = sparse(KK_sq);

end