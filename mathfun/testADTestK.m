% Francois Aguet, 03/02/2012

function testADTestK()

% Example given in Table 4 of
% [1] Scholz and Stephens, J. Am. Stat. Assoc. 82(399), 1987

samples = {[38.7 41.5 43.8 44.5 45.5 46.0 47.7 58.0],...
           [39.2 39.3 39.7 41.4 41.8 42.9 43.3 45.8],...
           [34.0 35.0 39.0 40.0 43.0 43.0 44.0 45.0],...
           [34.0 34.8 34.8 35.4 37.2 37.8 41.2 42.8]};
[hval T cval] = adtestk(samples)
