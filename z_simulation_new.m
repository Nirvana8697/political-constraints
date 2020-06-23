% This script is to simulate a discretized income process, using uniform
% distribution.


zsiz=201;
Z= dlmread('Z.txt');
Zprob=reshape(dlmread('Zprob.txt'),zsiz,zsiz);


initindex=1;
maxim=100000;

z_simul_loc=zeros(maxim,1);
z_simul_loc(1)=initindex;
iter=0;

for i=2:maxim
    trans_probs_1=Zprob(z_simul_loc(i-1),:);
    cumprob= cumsum(trans_probs_1);
    R=rand();
    k=find(cumprob>R,1);
    z_simul_loc(i)=k;
end
    dlmwrite('z_simul_loc.txt', z_simul_loc, 'delimiter', '\t', 'precision', 18);