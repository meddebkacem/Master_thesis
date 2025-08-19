M_el_C = zeros(size(M_el));

N_elts = size(cnctvty,1);
N_per_elt = size(cnctvty,2);
elts_cnct = zeros(5,N_per_elt,N_elts);
cnct_ind = zeros(5,N_per_elt,N_elts);

for lk=1:N_elts
    M_el_temp = zeros(size(M_el(:,:,1)));
    ord = cnctvty(lk,:);
    for jk=1:length(ord)
        x_ind = find(cnctvty==ord(jk));
        y_ind = fix(x_ind/N_elts)+1;
        y_ind(y_ind>N_per_elt) = N_per_elt;
        for jj=1:length(x_ind)
            if x_ind(jj)==N_elts
                continue
            elseif mod(x_ind(jj),N_elts)==0
                x_ind(jj) = x_ind(jj)/N_elts;
            else
                x_ind(jj) = mod(x_ind(jj),N_elts);
            end
        end
        % brdge1 = zeros(5,1);
        % brdge1(1:length(x_ind)) = x_ind;
        % brdge2 = zeros(5,1);
        % brdge2(1:length(y_ind)) = y_ind;
        % elts_cnct(:,jk,lk) = brdge1;
        % cnct_ind(:,jk,lk) = brdge2;
        for nn=1:length(x_ind)
            M_el_temp(jk,jk) = M_el_temp(jk,jk)+M_el(y_ind(nn),y_ind(nn),x_ind(nn));
        end
    end
    M_el_C(:,:,lk) = M_el_temp;
end

%%

for lkl=1:N_elts
    M_el_temp = zeros(size(M_el(:,:,1)));
    for kk=1:length(M_el(:,:,1))
        x_ind = nonzeros(elts_cnct(:,kk,lkl));
        y_ind = nonzeros(cnct_ind(:,kk,lkl));
        for nn=1:length(x_ind)
            M_el_temp(kk,kk) = M_el(y_ind(nn),y_ind(nn),x_ind(nn));
        end
    end

end
































% for lk=1:N_elts
%     temp = elts_cnct(:,:,lk);
%     temp(temp==lk) = 0;
%     elts_cnct(:,:,lk) = temp;
%     % elts_cnct(elts_cnct(:,:,lk)==lk) = 0;
% end