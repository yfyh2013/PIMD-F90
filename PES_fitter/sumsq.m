function err = sumsq(beta)

global x;
global menergy_c;

fitted_energies = fitFunc_pot_nasa(beta,x);

num = length(menergy_c);

%percent error
err = 100*sqrt(sum(  ((fitted_energies-menergy_c')./(menergy_c+.00001)').^2)/num );
disp(err)

end%function