function [e1,dr1] = pot_nasa(r1)
% calculates the energy (e1) and the derivatives (dr1) of a water monomer
% (coordinates in the array r1). The potential has been developed by
% "H. Partridge and D. W. Schwenke, J. Chem. Phys. 106, 4618 (1997)".
% Some extra code for the calculation of the energy derivatives has been
% added  by C. J. Burnham.
% translation to Matlab from Fortran by Michelle Fritz 

global idx idxD coefD c5zA cbasis ccore crest idxm cmass reoh thetae b1 ...
    roh alphaoh deohA phh1A phh2 f5z fbasis fcore frest a b c0 c1 c2 b1D
mnasa_mod; % declare global variables
ROH1 = zeros(1,3);
fmat = zeros(16,3);


ROH1(1:3) = r1(1:3,2) - r1(1:3,1);
ROH2(1:3) = r1(1:3,3) - r1(1:3,1);
RHH(1:3)  = r1(1:3,2) - r1(1:3,3);
dROH1 = sqrt(dot(ROH1, ROH1));
dROH2 = sqrt(dot(ROH2, ROH2));
dRHH  = sqrt(dot(RHH , RHH ));
costh = (ROH1(1)*ROH2(1) + ROH1(2)*ROH2(2) + ROH1(3)*ROH2(3) ) / (dROH1*dROH2);


c5z = f5z*c5zA + fbasis*cbasis + fcore*ccore + frest*crest;
%c5z = f5z*c5zA;
deoh = f5z*deohA;
phh1 = f5z*phh1A;
phh1 = phh1*exp(phh2);

costhe = -.24780227221366464506;

exp1 = exp(-alphaoh*(dROH1-roh));
exp2 = exp(-alphaoh*(dROH2-roh));
Va = deoh*(exp1*(exp1-2)+exp2*(exp2-2));
Vb  = phh1*exp(-phh2*dRHH);
dVa1= 2*alphaoh*deoh*exp1*(1-exp1)/dROH1;
dVa2= 2*alphaoh*deoh*exp2*(1-exp2)/dROH2;
dVb = -phh2*Vb/dRHH;
x1 = (dROH1-reoh)/reoh;
x2 = (dROH2-reoh)/reoh;
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
for  j=2:245
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
e1 = Va+Vb+Vc;
e1 = e1+0.44739574026257; % correction
%e1 = e1*0.00285914375100642899 % cm-1 --> Kcal/mol
e1 = e1*1.23981E-4; % cm-1 --> eV
%.... derivatives .............
dVcdr1 = (-2*b1*efac*(dROH1-reoh)*sum0 + efac*sum1/reoh)/dROH1;
dVcdr2 = (-2*b1*efac*(dROH2-reoh)*sum0 + efac*sum2/reoh)/dROH2;
dVcdcth = efac*sum3;

dr1(1:3,2) = dVa1*ROH1 + dVb*RHH + dVcdr1*ROH1 + ...
       dVcdcth*(ROH2(1:3)/(dROH1*dROH2)-costh*ROH1(1:3)/(dROH1*dROH1));
dr1(1:3,3) = dVa2*ROH2-dVb*RHH +dVcdr2*ROH2  + ...
       dVcdcth*(ROH1(1:3)/(dROH1*dROH2)-costh*ROH2(1:3)/(dROH2*dROH2));
dr1(1:3,1) = -(dr1(1:3,2)+dr1(1:3,3));
%dr1 = dr1*.00285914375100642899;
dr1 = dr1*1.23981E-4; % cm-1/Ang --> eV/Ang

end % function pot_nasa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q3,gradq] = dms_nasa(r1)
%****
%**** calculates the charges (q3) on the Oxygen and the Hydrogen sites that
%**** reproduce the dipole moment surface  of a water monomer
%**** (coordinates in the array r1). The model has been developed by
%**** "H. Partridge and D. W. Schwenke, J. Chem. Phys. 106, 4618 (1997)"
%**** while some extra code for the calculation of the deivatives of the 
%**** charges with respect to the atomic displacemets (gradq) has
%**** added  by C. J. Burnham.
%*** ------------------------------------------------------------------------
%***  In the version 3.0 of the TTM3-F model some additional terms have been
%***  used for the calculation of the charges (See the comments)
%*** ------------------------------------------------------------------------

global vdwC vdwD vdwE aDD aCCaCD polarM polfacO polfacH polfacM ...
    dms_param1 dms_param2 dms_param3 gammaM MAXITER diptol dmix ...
    CHARGECON DEBYE 
global idx idxD coefD c5zA cbasis ccore crest idxm cmass reoh thetae b1 ...
    roh alphaoh deohA phh1A phh2 f5z fbasis fcore frest a b c0 c1 c2 b1D

mnasa_mod;
ttm3f_mod;

R1 = zeros(3,3);
q3 = zeros(1,3);
gradq = zeros(3,3,3);
ROH1 = zeros(1,3);
ROH2 = zeros(1,3);
RHH = zeros(1,3);
AxB = zeros(1,3);
fmat = zeros(16,3);

deoh = f5z*deohA;
phh1 = f5z*phh1A;
phh1 = phh1*exp(phh2);

ath0 = 1.82400520401572996557;
costhe = -.24780227221366464506;

ROH1 = R1(:,2) - R1(:,1);
ROH2 = R1(:,3) - R1(:,1);
RHH  = R1(:,2) - R1(:,3);
dROH1 = dsqrt(dot_product(ROH1, ROH1));
dROH2 = dsqrt(dot_product(ROH2, ROH2));
dRHH  = dsqrt(dot_product(RHH, RHH));
costh = dot_product(ROH1, ROH2)/(dROH1*dROH2);
efac = dexp(-b1D*(  (dROH1-reoh)^2 + (dROH2-reoh)^2));
   
x1 = (dROH1-reoh)/reoh;
x2 = (dROH2-reoh)/reoh;
x3 = costh - costhe;
fmat(0,1:3) = 0;
fmat(1,1:3) = 1;
for j=2:15
   fmat(j,1) = fmat(j-1,1)*x1;
   fmat(j,2) = fmat(j-1,2)*x2;
   fmat(j,3) = fmat(j-1,3)*x3;
end
fmat(16,1:3) = 0;
%    Calculate dipole moment

P1 = 0; 
P2 = 0;
PL1 = costh;
PL2 = 0.5*(3*PL1*PL1-1);
dp1dr1 = 0;
dp1dr2 = 0;
dp1dcabc = 0;
dp2dr1 = 0;
dp2dr2 = 0;
dp2dcabc = 0;
for j=2:84
   inI = idxD(j,1);
   inJ = idxD(j,2);
   inK = idxD(j,3);
   P1 = P1 + coefD(j) * fmat(inI ,1)*fmat(inJ, 2)*fmat(inK, 3);
   P2 = P2 + coefD(j) * fmat(inJ ,1)*fmat(inI, 2)*fmat(inK, 3);
   dp1dr1 =dp1dr1+coefD(j)*(double(inI-1)*fmat(inI-1,1)*fmat(inJ,2)*fmat(inK,3));
   dp1dr2 =dp1dr2+coefD(j)*(double(inJ-1)*fmat(inI,1)*fmat(inJ-1,2)*fmat(inK,3));
   dp1dcabc=dp1dcabc+coefD(j)*(double(inK-1)*fmat(inI,1)*fmat(inJ,2)*fmat(inK-1,3));
   dp2dr1 = dp2dr1+coefD(j)*( double(inJ-1)*fmat(inJ-1,1)*fmat(inI,2)*fmat(inK,3));
   dp2dr2 = dp2dr2+coefD(j)*( double(inI-1)*fmat(inJ,1)*fmat(inI-1,2)*fmat(inK,3));
   dp2dcabc=dp2dcabc+coefD(j)*(double(inK-1)*fmat(inJ,1)*fmat(inI,2)*fmat(inK-1,3));
end
dp1dr1=dp1dr1 / reoh*0.529177249;
dp1dr2=dp1dr2 / reoh*0.529177249;
dp2dr1=dp2dr1 / reoh*0.529177249;
dp2dr2=dp2dr2 / reoh*0.529177249;
PC0 = A*((dROH1^B)+(dROH2^B))*(C0+PL1*C1+PL2*C2);
dpc0dr1=a*(b*dROH1^(b-1))*(c0+pl1*c1+pl2*c2)*0.529177249*0.529177249;
dpc0dr2=a*(b*dROH2^(b-1))*(c0+pl1*c1+pl2*c2)*0.529177249*0.529177249;
dpc0dcabc=a*((dROH1^b)+(dROH2^b))*(c1+0.5*(6.*pl1)*c2)*0.529177249;
defacdr1 =-2*b1D*(dROH1-reoh)*efac*0.529177249;
defacdr2 =-2*b1D*(dROH2-reoh)*efac*0.529177249;
dp1dr1=dp1dr1*efac+p1*defacdr1+dpc0dr1;
dp1dr2=dp1dr2*efac+p1*defacdr2+dpc0dr2;
dp1dcabc=dp1dcabc*efac+dpc0dcabc;
dp2dr1=dp2dr1*efac+p2*defacdr1+dpc0dr1;
dp2dr2=dp2dr2*efac+p2*defacdr2+dpc0dr2;
dp2dcabc=dp2dcabc*efac+dpc0dcabc;
P1 = coefD(1)+P1*efac+PC0*0.529177249;
P2 = coefD(1)+P2*efac+PC0*0.529177249;
bfac=1/0.529177249;
q3(1) = -(P1+P2);   % Oxygen
q3(2) = P1;   % Hydrogen-1
q3(3) = P2;   %  Hydrogen-2
D = r1(1:3,1)*q3(1) + r1(1:3,2)*q3(2) + r1(1:3,3)*q3(3)
dp1dr1=dp1dr1*bfac;
dp1dr2=dp1dr2*bfac;
dp2dr1=dp2dr1*bfac;
dp2dr2=dp2dr2*bfac;
%--------------------------------------------------------------------------------
%............. Modification of the gas-phase dipole moment surface...........
%--------------------------------------------------------------------------------
   AxB(1) = ROH1(2)*ROH2(3) - ROH1(3)*ROH2(2);
   AxB(2) =-ROH1(1)*ROH2(3) + ROH1(3)*ROH2(1);
   AxB(3) = ROH1(1)*ROH2(2) - ROH1(2)*ROH2(1);
   dAxB = dsqrt(dot_product(AxB, AxB) );
   sinth = dAxB / (dROH1*dROH2);
   ang = atan2(sinth, costh);
   P1 = dms_param1*(dROH1 -dms_param2) + dms_param3*(ang-ath0);
   P2 = dms_param1*(dROH2 -dms_param2) + dms_param3*(ang-ath0);
   q3(1) = q3(1) - (P1+P2);
   q3(2) = q3(2) + P1;
   q3(3) = q3(3) + P2;
   dp1dr1 = dp1dr1 + dms_param1;
   dp2dr2 = dp2dr2 + dms_param1;
   dp1dcabc = dp1dcabc   - dms_param3/sinth;
   dp2dcabc = dp2dcabc   - dms_param3/sinth;
%--------------------------------------------------------------------------------

f1q1r13=(dp1dr1-(dp1dcabc*costh/dROH1))/dROH1;
f1q1r23=dp1dcabc/(dROH1*dROH2);
f2q1r23=(dp1dr2-(dp1dcabc*costh/dROH2))/dROH2;
f2q1r13=dp1dcabc/(dROH2*dROH1);
f1q2r13=(dp2dr1-(dp2dcabc*costh/dROH1))/dROH1;
f1q2r23=dp2dcabc/(dROH1*dROH2);
f2q2r23=(dp2dr2-(dp2dcabc*costh/dROH2))/dROH2;
f2q2r13=dp2dcabc/(dROH2*dROH1);


%gradient of charge h1 wrt displacement of h1
gradq(1,1,1)=f1q1r13*ROH1(1)+f1q1r23*ROH2(1);
gradq(1,1,2)=f1q1r13*ROH1(2)+f1q1r23*ROH2(2);
gradq(1,1,3)=f1q1r13*ROH1(3)+f1q1r23*ROH2(3);
%gradient of charge h1 wrt displacement of h2
gradq(2,1,1)=f2q1r13*ROH1(1)+f2q1r23*ROH2(1);
gradq(2,1,2)=f2q1r13*ROH1(2)+f2q1r23*ROH2(2);
gradq(2,1,3)=f2q1r13*ROH1(3)+f2q1r23*ROH2(3);
%gradient of charge h1 wrt displacement of O
gradq(3,1,1)=-(gradq(1,1,1)+gradq(2,1,1));
gradq(3,1,2)=-(gradq(1,1,2)+gradq(2,1,2));
gradq(3,1,3)=-(gradq(1,1,3)+gradq(2,1,3));
%gradient of charge h2 wrt displacement of h1
gradq(1,2,1)=f1q2r13*ROH1(1)+f1q2r23*ROH2(1);
gradq(1,2,2)=f1q2r13*ROH1(2)+f1q2r23*ROH2(2);
gradq(1,2,3)=f1q2r13*ROH1(3)+f1q2r23*ROH2(3);
%gradient of charge h2 wrt displacement of h2
gradq(2,2,1)=f2q2r13*ROH1(1)+f2q2r23*ROH2(1);
gradq(2,2,2)=f2q2r13*ROH1(2)+f2q2r23*ROH2(2);
gradq(2,2,3)=f2q2r13*ROH1(3)+f2q2r23*ROH2(3);
%gradient of charge h2 wrt displacement of O
gradq(3,2,1)=-(gradq(1,2,1)+gradq(2,2,1));
gradq(3,2,2)=-(gradq(1,2,2)+gradq(2,2,2));
gradq(3,2,3)=-(gradq(1,2,3)+gradq(2,2,3));
%gradient of charge O wrt displacement of h1
gradq(1,3,1)=-(gradq(1,1,1)+gradq(1,2,1));
gradq(1,3,2)=-(gradq(1,1,2)+gradq(1,2,2));
gradq(1,3,3)=-(gradq(1,1,3)+gradq(1,2,3));
%gradient of charge O wrt displacement of h2
gradq(2,3,1)=-(gradq(2,1,1)+gradq(2,2,1));
gradq(2,3,2)=-(gradq(2,1,2)+gradq(2,2,2));
gradq(2,3,3)=-(gradq(2,1,3)+gradq(2,2,3));
%gradient of charge O wrt displacement of O
gradq(3,3,1)=-(gradq(3,1,1)+gradq(3,2,1));
gradq(3,3,2)=-(gradq(3,1,2)+gradq(3,2,2));
gradq(3,3,3)=-(gradq(3,1,3)+gradq(3,2,3));
    
end % function dms_nasa