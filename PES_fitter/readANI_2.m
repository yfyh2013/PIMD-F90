function [ra,energy] = readANI_2(file1,file2)
currentfolder = pwd;
pathstr = [currentfolder,'/'];
%pathstr = [currentfolder,'\results\'];
if ~isempty(file1)
file = strcat(pathstr,file1)
fid = fopen(file,'r');
end
if exist('file2')
file = strcat(pathstr,file2)
fid2 = fopen(file);
end

if ~isempty(file1)
ii = 0;
while ~feof(fid)
    ii = ii + 1;
    ip = ii;
    na = fscanf(fid,'%i', 1); % number of atoms
    junk = fscanf(fid,'%g', 1); 
    if isempty(na)
        break;
    end
        for ia = 1:na
            atom = fscanf(fid,'%s',1);
            a(ia) = {atom};
            ra(1:3,ip,ia) =  fscanf(fid,'%g %g %g',[1 3]); % atom position       
        end
end
else
    ra =[];
end
if exist('file2')
ie = 0;
while ~feof(fid2)
    ie = ie + 1;
    e = fscanf(fid2,'%g',1); % monomer energy  
    if isempty(e)
        break;
    end
     ip = ie;
     energy(ip) = e; % monomer energy
end
else
    energy = [];
end

end % function