%transfer vs amp scratch

amp = [20 23 27 30 25 18 24];
transfer_p =[85 88 86 77 89 77 87];
transfer_unc = [1 1 3 3 1 3 2];
figure(8)
errorbar(amp,transfer_p,transfer_unc,'kx')