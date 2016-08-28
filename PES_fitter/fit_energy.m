% PES fitting code by Michelle Fritz
% commented and edited slightly by D.C. Elton, 2016
clear all;
doFit = 0;
makePlots = 1;

global menergy_c x c5zA_fit roh alphaoh deohA phh1A phh2 x0 n_c5zA % declare global variables 
mnasa_mod; % declare global variables   
nc = size(c5zA_fit,2);
n_c5zA = nc;

%output filename
filename = 'test_new_parameters_custom_PBE_2916pts_robust.dat';%'new_parameters_custom_BH_stride4.dat'; %'new_parameters_custom_PBE_stride4.dat';

if (0) %% to use old parameters
    filename0 = filename ;
    beta0 = importdata(filename0);
else % start from scratch
    beta0 = bPack(c5zA_fit,roh,alphaoh,deohA,phh1A,phh2);
end

%--- read in data ----
%[ra1, menergy_c_all] = readANI_2('out_PBE_DZP.ANI','out_PBE_DZP.energies');
[ra1, menergy_c_all]  = readANI_2('larger_PBE_data/out_PBE_DZP.ANI','larger_PBE_data/out_PBE_DZP.energies');
%[ra1, menergy_c_all] = readANI_2('out_BH.ANI','out_BH.energies');

%--- downsample ------
nm = numel(menergy_c_all)
stride = 5;
for i=1:floor(nm/stride)
    ra(:,i,:) = ra1(:,i*stride,:);
    menergy_c(i) = menergy_c_all(i*stride);
end

menergy_c = 0.0433634*menergy_c;
nm = numel(menergy_c)

%--- create coords vector for equilibrium config ---
x0 = xPack(reshape(ra(1:3,1,1:3),3,3));

%--- calculate monomer energies from NASA PES ------
for im = 1:nm
    x(im,:) = xPack(reshape(ra(1:3,im,1:3),3,3)); %Angstroms
    menergy(im) = fitFunc_pot_nasa(beta0,x(im,:));
end

%--- monomer energies for equilibrium --------------
menergy0 = menergy(1);
menergy_c0 = menergy_c(1);

%--- subtract equilibrium values, calculate error---
error = 0;
for im = 1:nm
    menergy(im) = menergy(im) - menergy0;
    menergy_c(im) = menergy_c(im) - menergy_c0;
    error = error + (menergy_c(im) - menergy(im))^2;
end
error = sqrt(error/nm);
sse = sum((menergy - menergy_c).^2);

%--- sum squared error -----------------------------
fprintf('%s %f \n \n','SSE is ', sse);
fprintf('%s %f \n \n','RMS error in is ', error*8065.73/nm);

%menergy = menergy*8065.73; % eV -> 1/cm

%% ---------------- fitting ------------------------------------
if (doFit)

   options = statset('Display','iter','MaxIter',1000,'TolFun',1e-2,'FunValCheck','on','robust','on','RobustWgtFun','logistic');
   %options = optimset('Display','iter','MaxIter',5000,'FunValCheck','on','robust','on');
   %options.Algorithm = 'levenberg-marquardt';
   %options.RobustWgtFun = 'bisquare';

   %--- original method --- 
   b = nlinfit(x,menergy_c',@fitFunc_pot_nasa,beta0,options); 
   
   %--- minimization of residuals method 1 -------
   %b = fminsearch(@sumsq,beta0,options);
   
   %--- minimization of residuals method 2 -------
   % fp = @(b)fitFunc_pot_nasa(b,x);
   % b = lsqnonlin(fp,beta0,-1000000000,1000000000,options); 
   
   dlmwrite(filename,b,'precision','%.16f')
   
end%if (doFit)
%---------------- end fitting ----------------------------------------

%% ---------------- plotting -----------------------------------------
%--- calculate error ---
error = 0;
b = importdata(filename);
[c5zA_fit,roh,alphaoh,deohA,phh1A,phh2] = bUnpack(b);
for im = 1:nm
  %  b = importdata(filename)
    %size(xPack(reshape(ra(1:3,id,1:3),3,3))
    %menergy2(im) = pot_nasa_new(x(im,:),filename) - menergy2_0;
    %menergy2(im) = pot_nasa_fit(b,x(im,:)) - menergy2_0;
    menergy2(im) = fitFunc_pot_nasa(b,x(im,:)); 
    %menergy2(im) = pot_nasa_fit(beta0,x(im,:)) - pot_nasa_fit(beta0,x0);
    e_0 = menergy_c(im);
    e_f = menergy2(im);
    error = error + (e_f - e_0)^2;
    e_f2 = fitFunc_pot_nasa(b,x(im,:)); 
    x(im,:) = xPack(reshape(ra(1:3,im,1:3),3,3)); %Angstroms
end

error = sqrt(error/nm);
sse = sum((menergy2-menergy_c).^2); % sum squared error
fprintf('SSE is %f \n', sse);
fprintf('RMS error is %f \n',error);

%---- plots
if (makePlots)
data(1:nm,1) = ones(nm,1);
data(1:nm,2) = menergy_c;
b2 = regress(menergy2',data);
%---- basic regression plot 
h1=figure(2);
e0 = menergy_c;
ecalc = menergy2;
size(e0)
size(ecalc)
plot(e0,ecalc,'*',[min(e0) max(e0)],...
    [b2(1)+min(e0)*b2(2), b2(1)+b2(2)*max(e0)], ...
    [min(e0) max(e0)],[min(e0) max(e0)]);
legend('Calculated','Linear fit','Reference','location','SouthEast')
%plot(dev2*1.e4,error*1.e3,'o');
xlabel('Monomer energy calculated by Patridge and Schwenke (eV) ');
%legend('E-Emin','E-Ecalc');
ylabel('Predicted monomer energy (eV)');
saveas(h1,['energy_correlation_',num2str(nc+5),'.fig']) 
hold off;
%--------------------------------

%plot ranges
HOHangle = [90 115];
roh1_norm = [0.8 2];
roh2_norm = roh1_norm;

%num of points to plot in each dimension
n = 20; 

b_angle = [min(HOHangle):((max(HOHangle)-min(HOHangle))/(n-1)):max(HOHangle)];
b_roh1 = [min(roh1_norm):((max(roh1_norm)-min(roh1_norm))/(n-1)):max(roh1_norm)+1];
b_roh2 = [min(roh2_norm):((max(roh2_norm)-min(roh2_norm))/(n-1)):max(roh2_norm)+1];

pathstr = '/home/dan/Dropbox/RESEARCH/monPIMD_project/PES_fitting';

%outfile = strcat(pathstr,'monomer_mesh_fit_',num2str(n),'.ANI');
%fid = fopen(outfile,'w');
%energyfile = strcat(pathstr,'ener_mesh_fit_',num2str(n),'.dat');
%fid2 = fopen(energyfile,'w');

rm1 = zeros(n,n); rm2 = rm1;
em2 = rm1; am = rm1; rm = rm1;
em = rm1; em0 = rm1;
    
for i_ang = 1:n   
    for i_roh1 = 1:n
        for i_roh2 = 1:n
            ang = b_angle(i_ang);
            ang = ang*(pi/180);
            roh1 = b_roh1(i_roh1);
            roh2 = b_roh2(i_roh2);
            rhh=sqrt(roh1^2+roh2^2-2*roh1*roh2*cos(ang));
            ra(:,1) = [0 0 0]';
            ra(:,2) = roh1*[-sin(ang/2) -cos(ang/2) 0]';
            ra(:,3) = roh2*[sin(ang/2) -cos(ang/2) 0]';
            x = xPack(ra);
            e_mesh = pot_nasa_fit(b,x);
            e_mesh0 = pot_nasa(ra);
            rm1(i_roh1,i_roh2) = roh1;
            rm2(i_roh1,i_roh2) = roh2;
            em2(i_roh1,i_roh2)  = e_mesh/(1.23981E-4); % eV -> (1/cm)
            if (roh1 == roh2)
                am(i_ang,i_roh1) = ang*(180/pi);
                rm(i_ang,i_roh1) = roh1;
                em(i_ang,i_roh1) = e_mesh/(1.23981E-4); % eV -> (1/cm)
                em0(i_ang,i_roh1) = e_mesh0/(1.23981E-4); % eV -> (1/cm)
            end
            %fprintf(fid2, '%13.9f \n', e_mesh);
            %fprintf(fid, '%i \n \n', 3);
            %fprintf(fid, 'O \t %13.9f %13.9f %13.9f  \n',...
            %ra(1,1), ra(2,1), ra(3,1));
            %for k = 1:2
            %    fprintf(fid, 'H \t %13.9f %13.9f %13.9f  \n',...
            %        ra(1,1+k), ra(2,1+k), ra(3,1+k));
            %end
        end
    end
end
%% ---------------- make plots ------------------------------------
if(1)
h=figure(3);
ylabel('Energy');
xlabel('r(ang)');
hold on; 
for i=1:2:20
    plot(rm(i,:),em(i,:))
    plot(rm(i,:),em0(i,:),'-*')
end   

%--------------------------------------
h=figure(4); hold on; 
contour(rm1,rm2,em2,'LevelStep',1000)
ylabel('r_oh2(ang)');
xlabel('r_oh1(ang)');
title(['Energy surface difference with ',num2str(nc+5),' parameters (1/cm) ']);
colorbar;

%--------------------------------------
h=figure(5);
contour(rm,am,em,'LevelStep',1000)
ylabel('theta');
xlabel('r(ang)');
title(['Adjusted Patridge and Schwenke energy surface with ',num2str(nc+5),' parameters (1/cm) ']);
colorbar;
%saveas(h,[pathstr,'energy_surface_ps_',num2str(nc+5),'.fig']) 

%--------------------------------------
h2=figure(6);
contour(rm,am,em0,'LevelStep',1000)
ylabel('theta');
xlabel('r(ang)');
title('Patridge and Schwenke energy surface (1/cm) ');
colorbar;
saveas(h2,[pathstr,'energy_surface_reference_ps_',num2str(nc+5),'.fig']) 

%--------------------------------------
h3=figure(7);
contour(rm,am,em-em0,'LevelStep',1000)
ylabel('theta');
xlabel('r(ang)');
title(['Energy surface difference with ',num2str(nc+5),' parameters (1/cm) ']);
colorbar;
saveas(h3,[pathstr,'energy_surface_diff_ps_',num2str(nc+5),'.fig']) 


end
end% if makeplots 


%write out parameters in format for .f90 import
writeParameters(c5zA_fit,roh,alphaoh,deohA,phh1A,phh2,filename)