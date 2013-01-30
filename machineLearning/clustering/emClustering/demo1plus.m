% this is demo1.m
pp = [0.3333 0.3333];
mu1 = [0 2];
mu2 = [3 0];
mu3 = [-3 0];
mu = [mu1' mu2' mu3' ];
covar(:,:,1) = [1 0; 0 1];
covar(:,:,2) = [1 0; 0 1];
covar(:,:,3) = [1 0; 0 1];
y = genmix(900,mu,covar,pp);
clear mu mu1 mu2 mu3

iniCovar=zeros(2,2,20);
kmax=20;
for k=1:kmax
    iniCovar(:,:,k)=covar(:,:,1);
end

[bestk,bestpp,bestmu,bestcov,dl,countf]=mixtures4plus(y,kmax,iniCovar,4,'verbose',1);
