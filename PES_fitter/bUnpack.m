function [c5zA,roh,alphaoh,deohA,phh1A,phh2] = bUnpack(x)
    nc=245;
    deohA = x(1);
    phh1A = x(2);
    phh2 = x(3);
    roh = x(4);
    alphaoh = x(5);
    c5zA = x(6:5+nc);
    %c5zA = x(6:80);
    %c5zA = x(6:45);
end

