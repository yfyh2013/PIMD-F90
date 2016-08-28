function x = bPack(c5zA,roh,alphaoh,deohA,phh1A,phh2)
    global n_c5zA
    x(1) = deohA;
    x(2) = phh1A;
    x(3) = phh2;
    x(4) = roh;
    x(5) = alphaoh;
    x(6:(5+n_c5zA)) = c5zA(1:n_c5zA);
    %x(6:80) = c5zA(1:75);
    %x(6:45) = c5zA(1:40);
end

