ppc = monkey(1).day(1).correct.chan1PEC;
pfc = monkey(1).day(1).correct.chan3dPFC;

histogram(pfc(:,1343))


rng('default')
x = randn(30,4);
y = randn(30,4);
y(:,4) = sum(x,2); % introduce correlation

[r,p] = corr(pfc,ppc);
%corr(X,Y) returns a p1-by-p2 matrix containing the pairwise correlation 
%coefficient between each pair of columns in the n-by-p1 and n-by-p2 
%matrices X and Y.

sig = (p < 0.001);

imagesc(r)
colorbar