function r = xUnpack(x)
na = numel(x)/3;
r = zeros(3,na);
i = 0;
for ia = 1:na
    for ix = 1:3
        i = i+1;
        r(ix,ia)=x(i);
    end
end
end

