function [e1,dr1] = pot_nasa_fit(b,x)
% calculates the energy (e1) and the derivatives (dr1) of a water monomer
% (coordinates in the array r1). The potential has been developed by
% "H. Partridge and D. W. Schwenke, J. Chem. Phys. 106, 4618 (1997)".
% Some extra code for the calculation of the energy derivatives has been
% added  by C. J. Burnham.
% translation to Matlab from Fortran by Michelle Fritz 

% global idx 
% mnasa_mod;
%ic = 75;
nij = 8;
nijk = 14;
reoh = 0.958649;
thetae = 104.3475;
b1 = 2.0;
f5z = 0.99967788500000;
%sizeb = size(b)
[c5zA_fit,roh,alphaoh,deohA,phh1A,phh2] = bUnpack(b);
% roh = roh
% alphaoh = alphaoh
% deohA = deohA
% phh1A = phh1A
% phh2 = phh2
% roh=roh
% sizec_potnasa = size(c5zA_fit);
r1 = xUnpack(x);
%sizer = size(r1)
% idx = zeros(75,3);
% %idx(1,1:3) = [1 1 1];
%idx = idx
idx0(:,1)= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,  ...
       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  ...
       2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  ...
       3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,  ...
       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  ...
       4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,  ...
       4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,  ...
       6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,  ...
       6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5,  ...
       5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,  ...
       7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6,  ...
       6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9,  ...
       9, 9, 9, 9, 9];
idx0(:,2)= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  ...
       1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  ...
       2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,  ...
       2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,  ...
       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  ...
       2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,  ...
       3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,  ...
       1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,  ...
       2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,  ...
       4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,  ...
       2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,  ...
       4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,  ...
       1, 1, 1, 1, 1];
idx0(:,3)= [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,  ...
       6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,  ...
      12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,  ...
       6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,  ...
       2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,  ...
      11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,  ...
       9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,  ...
       9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,  ...
       1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,  ...
       3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,  ...
       7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,  ...
       4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,  ...
       3, 4, 5, 6, 7];
   
 ij2 = 0;
 ij2 = ij2 + 1;
 idx(ij2,1:3) = [1 1 1];
 for  j=2:245
   inI = idx0(j,1)-1;
   %if (inI <= nc_fit)
    inJ = idx0(j,2)-1;
    inK = idx0(j,3)-1;
    if (((inI + inJ) <= nij) & (inK <= (nijk-(inI+inJ))))
   % if ((inJ <= nc_fit) & (inK <= nc_fit))
       ij2 = ij2 + 1;
       idx(ij2,1)=idx0(j,1);
       idx(ij2,2)=idx0(j,2);
       idx(ij2,3)=idx0(j,3);
   % end
   end
 end  
ic = ij2;
ROH1 = zeros(1,3);
fmat = zeros(16,3);

ROH1(1:3) = r1(1:3,2) - r1(1:3,1);
ROH2(1:3) = r1(1:3,3) - r1(1:3,1);
RHH(1:3)  = r1(1:3,2) - r1(1:3,3);
dROH1 = sqrt(dot(ROH1, ROH1));
dROH2 = sqrt(dot(ROH2, ROH2));
dRHH  = sqrt(dot(RHH , RHH ));
costh = (ROH1(1)*ROH2(1) + ROH1(2)*ROH2(2) + ROH1(3)*ROH2(3) ) / (dROH1*dROH2);


%c5z = f5z*c5zA + fbasis*cbasis + fcore*ccore + frest*crest;
c5z = f5z*c5zA_fit;
deoh = f5z*deohA;
phh1 = f5z*phh1A;
phh1 = phh1*exp(phh2);

costhe = -.24780227221366464506;
costhe = cos(thetae*(pi/180));

exp1 = exp(-alphaoh*(dROH1-roh));
exp2 = exp(-alphaoh*(dROH2-roh));
Va = deoh*(exp1*(exp1-2)+exp2*(exp2-2));
Vb  = phh1*exp(-phh2*dRHH);
dVa1= 2*alphaoh*deoh*exp1*(1-exp1)/(dROH1);
dVa2= 2*alphaoh*deoh*exp2*(1-exp2)/(dROH2);
dVb = -phh2*Vb/dRHH;
x1 = (dROH1-reoh)/(reoh);
x2 = (dROH2-reoh)/(reoh);
x3 = costh - costhe;
%fmat(0,1:3) = 0;
fmat(1,1:3) = 1;
for j=2:15
   fmat(j,1) = fmat(j-1,1)*x1;
   fmat(j,2) = fmat(j-1,2)*x2;
   fmat(j,3) = fmat(j-1,3)*x3;
end
fmat(16,1:3) = 0;

efac = exp(-b1*(  (dROH1-reoh)^2 + (dROH2-reoh)^2));

sum0 = 0; sum1 = 0; sum2 = 0; sum3 = 0;
for  j=2:ic
   inI = idx(j,1);
   inJ = idx(j,2);
   inK = idx(j,3);
   sum0 = sum0 + c5z(j) *  ( fmat(inI,1)*fmat(inJ,2) +   ...
        fmat(inJ,1)*fmat(inI,2)) *fmat(inK,3);
   findI = inI-1;
   if (findI == 0)
       findI = 16;
   end
   findJ = inJ-1;
   if (findJ == 0)
       findJ = 16;
   end
   findK = inK-1;
   if (findK == 0)
       findK = 16;
   end
   sum1 = sum1 + c5z(j) *  ( double(inI-1)*fmat(findI,1)*fmat(inJ,2) +   ...
        double(inJ-1)*fmat(findJ,1)*fmat(inI,2)  )*fmat(inK,3); 
   sum2 = sum2 + c5z(j) *  ( double(inJ-1)*fmat(inI,1)*fmat(findJ,2) +   ...
        double(inI-1)*fmat(inJ,1)*fmat(findI,2)  )*fmat(inK,3); 
   sum3 = sum3 + c5z(j) *  ( fmat(inI,1)*fmat(inJ,2) +   ...
        fmat(inJ,1)*fmat(inI,2)) * double(inK-1)*fmat(findK,3); 
end     

%.... energy..........
Vc= 2*c5z(1)+efac*sum0;
% c5z1 = c5z(1)
% efac = efac
% sum0 = sum0
%  Va = Va
%  Vb = Vb
%  Vc = Vc
e1 = Va+Vb+Vc;
%e1 = e1+0.44739574026257; % correction
%e1 = e1*0.00285914375100642899 % cm-1 --> Kcal/mol
e1 = e1*1.23981E-4; % cm-1 --> eV
%e1 = e1 + 4.5975;
%e1 = e1;
%.... derivatives .............
dVcdr1 = (-2*b1*efac*(dROH1-reoh)*sum0 + efac*sum1/(reoh))/(dROH1);
dVcdr2 = (-2*b1*efac*(dROH2-reoh)*sum0 + efac*sum2/(reoh))/(dROH2);
dVcdcth = efac*sum3;

dr1(1:3,2) = dVa1*ROH1 + dVb*RHH + dVcdr1*ROH1 + ...
       dVcdcth*(ROH2(1:3)/(dROH1*dROH2)-costh*ROH1(1:3)/(dROH1*dROH1));
dr1(1:3,3) = dVa2*ROH2-dVb*RHH +dVcdr2*ROH2  + ...
       dVcdcth*(ROH1(1:3)/(dROH1*dROH2)-costh*ROH2(1:3)/(dROH2*dROH2));
dr1(1:3,1) = -(dr1(1:3,2)+dr1(1:3,3));
%dr1 = dr1*.00285914375100642899;
dr1 = dr1*1.23981E-4; %cm-1*Ang^-1 --> eV/Ang

end % function pot_nasa