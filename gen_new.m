%clear; clc;

noStates = 201;
uncondMean = 1.0;
coefLags = 0.949703;
stds = 0.023092;
bsiz = 100;
tsiz= 100;



[Z, Zprob] = tauchen1d(noStates, uncondMean, coefLags, stds, 4.2);

dlmwrite('Z.txt', Z, 'delimiter', '\t', 'precision', 18);
dlmwrite('Zprob.txt', Zprob, 'delimiter', '\t', 'precision', 18);

tgrid = linspace(0.0,0.60,tsiz);
dlmwrite('t.txt', tgrid, 'delimiter', '\t', 'precision', 18)

bgrid = linspace(-20,0,bsiz);
dlmwrite('b.txt', bgrid, 'delimiter', '\t', 'precision', 18);
