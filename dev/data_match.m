l_t=phi_logs(:,2);
d_t=data.time_create_write(:,1);
phi_log_matched = zeros(length(d_t),1).*nan;
phi_log_check = zeros(length(d_t),1);
for ii = 1: length(l_t)
phi_c = phi_logs(ii,3);
l_c=l_t(ii);
t_mask=l_c+17<d_t & l_c+23>d_t;
d_indx=find(t_mask);
phi_log_matched(d_indx,1) = phi_c;
phi_log_check(d_indx,1) = 1;

end