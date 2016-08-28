function f = fitFunc_pot_nasa(b,xx)
global x0
nm = size(xx,1);
f = zeros(nm,1);
f0 = pot_nasa_fit(b,x0);

for im = 1:nm
    f(im) = pot_nasa_fit(b,xx(im,:)) - f0;
end
end%function