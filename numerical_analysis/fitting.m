d = 2
close all
filename = strcat('../temperature_scan_data/',num2str(d),'D_temp_data_40x40.csv')
%filename = strcat('../temperature_scan_data/',num2str(d),'D_temp_data.csv')

r_offset = 1;
M = csvread(filename, r_offset);
temp = 1./M(:,end);

%get critical temperature
T_crit = 2.34;
[c, index] = min(abs(temp-T_crit));

temp_less_than_crit = temp(index:end);
mag_l = M(index:end,2);
chi_l = M(index:end,4);
cv_l = M(index:end,5);

mag = M(:,2);
chi = M(:,4);
cv = M(:,5);

plot(temp, mag, 'linewidth', 2)

figure()
plot(temp, chi, '.-', 'linewidth', 2, 'markersize', 30)
title(strcat('susceptibility for d=', num2str(d)))
xlabel('temperature (k_bT)')

figure()
plot(temp, cv, '.-', 'linewidth', 2, 'markersize', 30)
title(strcat('heat capacity for d=', num2str(d)))
xlabel('temperature (k_bT)')

figure()
subplot(121)
plot(temp, chi, '.-', 'linewidth', 2, 'markersize', 30)
title(strcat('susceptibility for d=', num2str(d)))
xlabel('temperature (k_bT)')
set(gca, 'fontsize', 16)
subplot(122)
plot(temp, cv, '.-', 'linewidth', 2, 'markersize', 30)
title(strcat('heat capacity for d=', num2str(d)))
xlabel('temperature (k_bT)')
set(gca, 'fontsize', 16)

cftool

