%Writes out parameters in a .f90 (Fortran 90) format for transfer into
%SIESTA m_water.f90 code 
%written by Michelle Fritz 
function writeParameters(c5zA_fit,roh,alphaoh,deohA,phh1A,phh2,filename)
file = strcat(filename(1:size(filename,2)-4),'_formated.dat');

fid = fopen(file,'w');
nc = numel(c5zA_fit);
reoh = 0.958649;
thetae = 104.3475;
b1 = 2.0;

fprintf(fid, '%s%i%s \n', 'data c5zA(1:',nc,')/  &');
for ic = 1:nc
    fprintf(fid, '%21.13E,',c5zA_fit(ic));
    if (ic == nc)
        fprintf(fid, '%s \n','/');
    elseif (mod(ic,3) == 0)
        fprintf(fid, '& \n');
    end
end
fprintf(fid, '\n \n %s', ...
    'data reoh,thetae,b1,roh,alphaoh,deohA,phh1A,phh2/');
fprintf(fid, '%10.6E, %s \n',reoh,'&');
fprintf(fid, '%10.6E,%5.1E,%21.13E,%21.13E, & \n',thetae,b1,roh,alphaoh);
fprintf(fid, '%21.13E,%21.13E,%21.13E/',deohA,phh1A,phh2);


%idx = zeros(nc,3);
[idx(:,1),idx(:,2),idx(:,3)] = idx_test;
for ix = 1:3
fprintf(fid, '\n \n %s%i%s \n', 'data idx(:,',ix,')/ &');
for ic = 1:nc
    fprintf(fid, '%i,',idx(ic,ix));
    if (ic == nc)
        fprintf(fid, '%s \n','/');
    elseif (mod(ic,20) == 0)
        fprintf(fid, '& \n');
    end
end
end

fclose(fid);
end %function