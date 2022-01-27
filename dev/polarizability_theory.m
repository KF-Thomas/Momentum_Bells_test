function pol = polarizability_theory(omega)
%theorectical dynamic polarizability of the helium tom at arbitrary
%frequency

hebec_constants
transition_data = 'C:\Users\BEC Machine\cloudstor\PROJECTS\Momentum_Bells_test\docs\he_transitions_nist.csv';
imported_data = readtable(transition_data);

fik = imported_data{:,4};
f_trans = const.c./imported_data{:,1}.*1e9;
Edel = imported_data{:,6}-imported_data{:,5};
SI_conversion = 1./1.64877727436e-41;% convert to atomic units
% pol = nansum(fik./(Edel.^2.*const.electron.^2-const.h.^2.*omega.^2)).*SI_conversion;
pol = 1./(2*pi).^2.*const.electron.^2./const.me.*nansum(fik./(f_trans.^2-omega.^2)).*SI_conversion;


end
