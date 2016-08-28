function mnasa_mod
    
% Parameters used for the calculation of the potential energy and the
% Dipole moment of the water monomer. Taken from:
% "H. Partridge and D. W. Schwenke, J. Chem. Phys. 106, 4618 (1997)"
% translation to Matlab from Fortran by Michelle Fritz 

global idx idxD coefD c5zA c5zA_fit cbasis ccore crest idxm cmass reoh thetae b1 ...
    roh alphaoh deohA phh1A phh2 f5z fbasis fcore frest a b c0 c1 c2 b1D

idx(:,1)= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,  ...
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
idx(:,2)= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  ...
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
idx(:,3)= [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,  ...
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

c5zA(1:245)= [4.2278462684916e+04, 4.5859382909906e-02, 9.4804986183058e+03, ...
       7.5485566680955e+02, 1.9865052511496e+03, 4.3768071560862e+02, ...
       1.4466054104131e+03, 1.3591924557890e+02,-1.4299027252645e+03, ...
       6.6966329416373e+02, 3.8065088734195e+03,-5.0582552618154e+02, ...
      -3.2067534385604e+03, 6.9673382568135e+02, 1.6789085874578e+03, ...
      -3.5387509130093e+03,-1.2902326455736e+04,-6.4271125232353e+03, ...
      -6.9346876863641e+03,-4.9765266152649e+02,-3.4380943579627e+03, ...
       3.9925274973255e+03,-1.2703668547457e+04,-1.5831591056092e+04, ...
       2.9431777405339e+04, 2.5071411925779e+04,-4.8518811956397e+04, ...
      -1.4430705306580e+04, 2.5844109323395e+04,-2.3371683301770e+03, ...
       1.2333872678202e+04, 6.6525207018832e+03,-2.0884209672231e+03, ...
      -6.3008463062877e+03, 4.2548148298119e+04, 2.1561445953347e+04, ...
      -1.5517277060400e+05, 2.9277086555691e+04, 2.6154026873478e+05, ...
      -1.3093666159230e+05,-1.6260425387088e+05, 1.2311652217133e+05, ...
      -5.1764697159603e+04, 2.5287599662992e+03, 3.0114701659513e+04, ...
      -2.0580084492150e+03, 3.3617940269402e+04, 1.3503379582016e+04, ...
      -1.0401149481887e+05,-6.3248258344140e+04, 2.4576697811922e+05, ...
       8.9685253338525e+04,-2.3910076031416e+05,-6.5265145723160e+04, ...
       8.9184290973880e+04,-8.0850272976101e+03,-3.1054961140464e+04, ...
      -1.3684354599285e+04, 9.3754012976495e+03,-7.4676475789329e+04, ...
      -1.8122270942076e+05, 2.6987309391410e+05, 4.0582251904706e+05, ...
      -4.7103517814752e+05,-3.6115503974010e+05, 3.2284775325099e+05, ...
       1.3264691929787e+04, 1.8025253924335e+05,-1.2235925565102e+04, ...
      -9.1363898120735e+03,-4.1294242946858e+04,-3.4995730900098e+04, ...
       3.1769893347165e+05, 2.8395605362570e+05,-1.0784536354219e+06, ...
      -5.9451106980882e+05, 1.5215430060937e+06, 4.5943167339298e+05, ...
      -7.9957883936866e+05,-9.2432840622294e+04, 5.5825423140341e+03, ...
       3.0673594098716e+03, 8.7439532014842e+04, 1.9113438435651e+05, ...
      -3.4306742659939e+05,-3.0711488132651e+05, 6.2118702580693e+05, ...
      -1.5805976377422e+04,-4.2038045404190e+05, 3.4847108834282e+05, ...
      -1.3486811106770e+04, 3.1256632170871e+04, 5.3344700235019e+03, ...
       2.6384242145376e+04, 1.2917121516510e+05,-1.3160848301195e+05, ...
      -4.5853998051192e+05, 3.5760105069089e+05, 6.4570143281747e+05, ...
      -3.6980075904167e+05,-3.2941029518332e+05,-3.5042507366553e+05, ...
       2.1513919629391e+03, 6.3403845616538e+04, 6.2152822008047e+04, ...
      -4.8805335375295e+05,-6.3261951398766e+05, 1.8433340786742e+06, ...
       1.4650263449690e+06,-2.9204939728308e+06,-1.1011338105757e+06, ...
       1.7270664922758e+06, 3.4925947462024e+05,-1.9526251371308e+04, ...
      -3.2271030511683e+04,-3.7601575719875e+05, 1.8295007005531e+05, ...
       1.5005699079799e+06,-1.2350076538617e+06,-1.8221938812193e+06, ...
       1.5438780841786e+06,-3.2729150692367e+03, 1.0546285883943e+04, ...
      -4.7118461673723e+04,-1.1458551385925e+05, 2.7704588008958e+05, ...
       7.4145816862032e+05,-6.6864945408289e+05,-1.6992324545166e+06, ...
       6.7487333473248e+05, 1.4361670430046e+06,-2.0837555267331e+05, ...
       4.7678355561019e+05,-1.5194821786066e+04,-1.1987249931134e+05, ...
       1.3007675671713e+05, 9.6641544907323e+05,-5.3379849922258e+05, ...
      -2.4303858824867e+06, 1.5261649025605e+06, 2.0186755858342e+06, ...
      -1.6429544469130e+06,-1.7921520714752e+04, 1.4125624734639e+04, ...
      -2.5345006031695e+04, 1.7853375909076e+05,-5.4318156343922e+04, ...
      -3.6889685715963e+05, 4.2449670705837e+05, 3.5020329799394e+05, ...
       9.3825886484788e+03,-8.0012127425648e+05, 9.8554789856472e+04, ...
       4.9210554266522e+05,-6.4038493953446e+05,-2.8398085766046e+06, ...
       2.1390360019254e+06, 6.3452935017176e+06,-2.3677386290925e+06, ...
      -3.9697874352050e+06,-1.9490691547041e+04, 4.4213579019433e+04, ...
       1.6113884156437e+05,-7.1247665213713e+05,-1.1808376404616e+06, ...
       3.0815171952564e+06, 1.3519809705593e+06,-3.4457898745450e+06, ...
       2.0705775494050e+05,-4.3778169926622e+05, 8.7041260169714e+03, ...
       1.8982512628535e+05,-2.9708215504578e+05,-8.8213012222074e+05, ...
       8.6031109049755e+05, 1.0968800857081e+06,-1.0114716732602e+06, ...
       1.9367263614108e+05, 2.8678295007137e+05,-9.4347729862989e+04, ...
       4.4154039394108e+04, 5.3686756196439e+05, 1.7254041770855e+05, ...
      -2.5310674462399e+06,-2.0381171865455e+06, 3.3780796258176e+06, ...
       7.8836220768478e+05,-1.5307728782887e+05,-3.7573362053757e+05, ...
       1.0124501604626e+06, 2.0929686545723e+06,-5.7305706586465e+06, ...
      -2.6200352535413e+06, 7.1543745536691e+06,-1.9733601879064e+04, ...
       8.5273008477607e+04, 6.1062454495045e+04,-2.2642508675984e+05, ...
       2.4581653864150e+05,-9.0376851105383e+05,-4.4367930945690e+05, ...
       1.5740351463593e+06, 2.4563041445249e+05,-3.4697646046367e+03, ...
      -2.1391370322552e+05, 4.2358948404842e+05, 5.6270081955003e+05, ...
      -8.5007851251980e+05,-6.1182429537130e+05, 5.6690751824341e+05, ...
      -3.5617502919487e+05,-8.1875263381402e+02,-2.4506258140060e+05, ...
       2.5830513731509e+05, 6.0646114465433e+05,-6.9676584616955e+05, ...
       5.1937406389690e+05, 1.7261913546007e+05,-1.7405787307472e+04, ...
      -3.8301842660567e+05, 5.4227693205154e+05, 2.5442083515211e+06, ...
      -1.1837755702370e+06,-1.9381959088092e+06,-4.0642141553575e+05, ...
       1.1840693827934e+04,-1.5334500255967e+05, 4.9098619510989e+05, ...
       6.1688992640977e+05, 2.2351144690009e+05,-1.8550462739570e+06, ...
       9.6815110649918e+03,-8.1526584681055e+04,-8.0810433155289e+04, ...
       3.4520506615177e+05, 2.5509863381419e+05,-1.3331224992157e+05, ...
      -4.3119301071653e+05,-5.9818343115856e+04, 1.7863692414573e+03, ...
       8.9440694919836e+04,-2.5558967650731e+05,-2.2130423988459e+04, ...
       4.4973674518316e+05,-2.2094939343618e+05];
   
c5zA_fit = c5zA;

%     expansion coefficients for basis correction

cbasis(1:245)= [6.9770019624764e-04,-2.4209870001642e+01, 1.8113927151562e+01, ...
       3.5107416275981e+01,-5.4600021126735e+00,-4.8731149608386e+01, ...
       3.6007189184766e+01, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -7.7178474355102e+01,-3.8460795013977e+01,-4.6622480912340e+01, ...
       5.5684951167513e+01, 1.2274939911242e+02,-1.4325154752086e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00,-6.0800589055949e+00, ...
       8.6171499453475e+01,-8.4066835441327e+01,-5.8228085624620e+01, ...
       2.0237393793875e+02, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       3.3525582670313e+02, 7.0056962392208e+01,-4.5312502936708e+01, ...
      -3.0441141194247e+02, 2.8111438108965e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-1.2983583774779e+02, 3.9781671212935e+01, ...
      -6.6793945229609e+01,-1.9259805675433e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-8.2855757669957e+02,-5.7003072730941e+01, ...
      -3.5604806670066e+01, 9.6277766002709e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 8.8645622149112e+02,-7.6908409772041e+01, ...
       6.8111763314154e+01, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       2.5090493428062e+02,-2.3622141780572e+02, 5.8155647658455e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 2.8919570295095e+03, ...
      -1.7871014635921e+02,-1.3515667622500e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-3.6965613754734e+03, 2.1148158286617e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00,-1.4795670139431e+03, ...
       3.6210798138768e+02, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -5.3552886800881e+03, 3.1006384016202e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 1.6241824368764e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 4.3764909606382e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 1.0940849243716e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 3.0743267832931e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00];

%     expansion coefficients for core correction

ccore(1:245)= [2.4332191647159e-02,-2.9749090113656e+01, 1.8638980892831e+01, ...
      -6.1272361746520e+00, 2.1567487597605e+00,-1.5552044084945e+01, ...
       8.9752150543954e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -3.5693557878741e+02,-3.0398393196894e+00,-6.5936553294576e+00, ...
       1.6056619388911e+01, 7.8061422868204e+01,-8.6270891686359e+01, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00,-3.1688002530217e+01, ...
       3.7586725583944e+01,-3.2725765966657e+01,-5.6458213299259e+00, ...
       2.1502613314595e+01, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       5.2789943583277e+02,-4.2461079404962e+00,-2.4937638543122e+01, ...
      -1.1963809321312e+02, 2.0240663228078e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-6.2574211352272e+02,-6.9617539465382e+00, ...
      -5.9440243471241e+01, 1.4944220180218e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-1.2851139918332e+03,-6.5043516710835e+00, ...
       4.0410829440249e+01,-6.7162452402027e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 1.0031942127832e+03, 7.6137226541944e+01, ...
      -2.7279242226902e+01, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -3.3059000871075e+01, 2.4384498749480e+01,-1.4597931874215e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 1.6559579606045e+03, ...
       1.5038996611400e+02,-7.3865347730818e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-1.9738401290808e+03,-1.4149993809415e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00,-1.2756627454888e+02, ...
       4.1487702227579e+01, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -1.7406770966429e+03,-9.3812204399266e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-1.1890301282216e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 2.3723447727360e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-1.0279968223292e+03, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 5.7153838472603e+02, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00];

%     expansion coefficients for v rest

crest(1:245)= [0.0000000000000e+00,-4.7430930170000e+00,-1.4422132560000e+01, ...
      -1.8061146510000e+01, 7.5186735000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
      -2.7962099800000e+02, 1.7616414260000e+01,-9.9741392630000e+01, ...
       7.1402447000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00,-7.8571336480000e+01, ...
       5.2434353250000e+01, 7.7696745000000e+01, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       1.7799123760000e+02, 1.4564532380000e+02, 2.2347226000000e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-4.3823284100000e+02,-7.2846553000000e+02, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00,-2.6752313750000e+02, 3.6170310000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00, ...
       0.0000000000000e+00, 0.0000000000000e+00];

%     expansion indicies for mass correction

       idxm = [1,2,1,1,3,2,1,2,1, ...
                 2,1,1,3,1,2,2,1,1, ...
                 1,1,2,1,1,1,2,2,3];

%     expansion coefficients for mass correction

       cmass= [-8.3554183e+00,3.7036552e+01,-5.2722136e+00, ...
            1.6843857e+01,-7.0929741e+01,5.5380337e+00,-2.9962997e+01, ...
            1.3637682e+02,-3.0530195e+00];

%     two body parameters

       reoh = 0.958649;
       thetae = 104.3475;
       b1 = 2.0;
       roh = 0.9519607159623009;
       alphaoh = 2.587949757553683;
       deohA = 42290.92019288289;
       phh1A = 16.94879431193463;
       phh2 = 12.66426998162947;

%     scaling factors for contributions to emperical potential

       f5z = 0.99967788500000;
       fbasis = 0.15860145369897;
       fcore = -1.6351695982132;
       frest = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxD(1:84,1) = [1, 1, 1, 2, 1, 1, 1, 2, 2, 3, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, ...
       1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 1, 1, 1, 1, 1, ...
       1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 1, 1, 1, 1, ...
       1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, ...
       5, 6, 6, 7];
idxD(1:84,2) =  [1, 1, 2, 1, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, ...
       1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 5, ...
       6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, ...
       5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, ...
       3, 1, 2, 1];
idxD(1:84,3) = [1, 2, 1, 1, 3, 2, 1, 2, 1, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, ...
       5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 6, 5, 4, 3, 2, ...
       1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 7, 6, 5, 4, ...
       3, 2, 1, 6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, ...
       1, 2, 1, 1];
coefD(1:84) = [-2.1689686086730e-03, 1.4910379754728e-02, 5.3546078430060e-02, ...
      -7.4055995388666e-02,-3.7764333017616e-03, 1.4089887256484e-01, ...
      -6.2584207687264e-02,-1.1260393113022e-01,-5.7824159269319e-02, ...
       1.4360743650655e-02,-1.5469680141070e-02,-1.3036350092795e-02, ...
       2.7515837781556e-02, 1.4098478875076e-01,-2.7663168397781e-02, ...
      -5.2378176254797e-03,-1.0237198381792e-02, 8.9571999265473e-02, ...
       7.2920263098603e-03,-2.6873260551686e-01, 2.0220870325864e-02, ...
      -7.0764766270927e-02, 1.2140640273760e-01, 2.0978491966341e-02, ...
      -1.9443840512668e-01, 4.0826835370618e-02,-4.5365190474650e-02, ...
       6.2779900072132e-02,-1.3194351021000e-01,-1.4673032718563e-01, ...
       1.1894031277247e-01,-6.4952851564679e-03, 8.8503610374493e-02, ...
       1.4899437409291e-01, 1.3962841511565e-01,-2.6459446720450D-02, ...
      -5.0128914532773D-02, 1.8329676428116D-01,-1.5559089125095D-01, ...
      -4.0176879767592D-02, 3.6192059996636D-01, 1.0202887240343D-01, ...
       1.9318668580051D-01,-4.3435977107932D-01,-4.2080828803311D-02, ...
       1.9144626027273D-01,-1.7851138969948D-01, 1.0524533875070D-01, ...
      -1.7954071602185D-02, 5.2022455612120D-02,-2.8891891146828D-01, ...
      -4.7452036576319D-02,-1.0939400546289D-01, 3.5916564473568D-01, ...
      -2.0162789820172D-01,-3.5838629543696D-01, 5.6706523551202D-03, ...
       1.3849337488211D-01,-4.1733982195604D-01, 4.1641570764241D-01, ...
      -1.2243429796296D-01, 4.7141730971228D-02,-1.8224510249551D-01, ...
      -1.8880981556620D-01,-3.1992359561800D-01,-1.8567550546587D-01, ...
       6.1850530431280D-01,-6.1142756235141D-02,-1.6996135584933D-01, ...
       5.4252879499871D-01, 6.6128603899427D-01, 1.2107016404639D-02, ...
      -1.9633639729189D-01, 2.7652059420824D-03,-2.2684111109778D-01, ...
      -4.7924491598635D-01, 2.4287790137314D-01,-1.4296023329441D-01, ...
       8.9664665907006D-02,-1.4003228575602D-01,-1.3321543452254D-01,...
      -1.8340983193745D-01, 2.3426707273520D-01, 1.5141050914514D-01];
b1D = 1;
a = 0.2999;
b = -0.6932;
c0 = 1.0099;
c1 = -0.1801;
c2 = 0.0892;

end %function
