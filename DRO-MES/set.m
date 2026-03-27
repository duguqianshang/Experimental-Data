function[wc_xj,wc_sj,sj_cz,xj_cz,wc_mean] = support_set(wc,K);
wc_mean=mean(wc);
wc_fc=var(wc);
wc_fc_2=repmat(wc_fc_2, K, 1);
theta=(wc-wc_mean).*wc_fc_2;
abs_theta=abs(theta); 
[sorted_values, ~] = sort(abs_theta, 'descend');
second_largest_values = sorted_values(500, :);
L_xj=-second_largest_values;
L_sj=second_largest_values;
wc_xj=wc_fc_3.*L_xj+wc_mean;
wc_sj=wc_fc_3.*L_sj+wc_mean;
sj_cz=wc_sj-wc;
sj_cz=abs(sj_cz);
xj_cz=wc_xj-wc;
xj_cz=abs(xj_cz);
end