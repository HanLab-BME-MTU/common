function ljs = localJensenShannon(u,s,w)

n = numel(u);
ljs = nan(n,1);

for j = 1+w:n-w    
    ljs(j) = jensenShannonDivergenceSpec(u(j-w:j+w),s(j-w:j+w));
end