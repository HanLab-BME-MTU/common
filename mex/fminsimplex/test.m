% File: test.m
% ------------
% This script finds local minima by feeding multiple initial guesses to the
% mexFunction fminsimplex.
%
% Modifications
% (1) Initial length is determined INSIDE the for loops.
% (2) Edge length vector of initial simplex is specified by the user, so 
%     that the size of initial simplex is constant for all initial guesses.
% (3) Parameters are normalized to the same scale.
%
% Pei-hsin Hsu
% April 28, 2009


load ref;

% reference model parameters
v = [4, -4];
a = [1, 1];
b = [10, 0];

% coordinate parameters
g2sSlope = 5:5:50;
s2gRate = 0.25:0.25:2.5;

% simulation parameters
model.simParam.time = 2000;                 % total simulation time [s]
model.simParam.dT = 0.10;                   % time step [s]
model.lenOffset = 0;
model.nStates = 2;

% growth
model.state(1).name = 'g';
model.state(1).speed = v(1);
model.state(1).nTransits = 1;
model.state(1).transit(1).dS = 2;
model.state(1).transit(1).rate(1) = a(1);
%model.state(1).transit(1).rate(2) = b(1);

% shrinkage
model.state(2).name = 's';
model.state(2).speed = v(2);
model.state(2).nTransits = 1;
model.state(2).transit(1).dS = 1;
%model.state(2).transit(1).rate(1) = a(2);
model.state(2).transit(1).rate(2) = b(2);

% initial condition
model.initCond.state = 1;

% mandatory fminsimplex parameters
model.nParam = 2;
model.objective = 'myfun';
model.ref = ref(3, 2); % ARMA(2, 1) model

model.scale = [20; 1];


% edge length vector of initial simplex
e = 0.5 * ones(2, 1);


for i = 1:length(g2sSlope)
    for j = 1:length(s2gRate)
        
        x = ones(2, 1);
        x(1) = g2sSlope(i);
        x(2) = s2gRate(j);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initial MT length is determined here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model.initCond.len = model.lenOffset + ...
            -(v(1)*x(2) + v(2)*a(1)) / (v(1)*b(2) + v(2)*x(1));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalize parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x(1) = x(1)/model.scale(1);
        x(2) = x(2)/model.scale(2);
        
        
        fprintf('\nStart point: [%f, %f]\n\n', x(1), x(2));

        [y, fval] = fminsimplex(x, model, e);

        fprintf('\nFinal point: [%f, %f], fval = %f\n\n', y(1), y(2), fval);
        
        init(i, j).x = x;
        init(i, j).y = y;
        inti(i, j).fval = fval;
    end
end

save('init', 'init');


