function x = xPack(r)
na = size(r,2);
i = 0;
x = zeros(1,3*na);
for ia = 1:na
    for ix = 1:3
        i = i+1;
        x(i) = r(ix,ia);
    end
end
end

