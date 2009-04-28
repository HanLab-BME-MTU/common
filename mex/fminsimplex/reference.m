% File: reference.m
% -----------------
% This script file is to generate a reference model.
%
% Pei-hsin Hsu
% April 27, 2009


v = [4, -4];
a = [1, 1];
b = [10, 0];

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
model.state(1).transit(1).rate(2) = b(1);

% shrinkage
model.state(2).name = 's';
model.state(2).speed = v(2);
model.state(2).nTransits = 1;
model.state(2).transit(1).dS = 1;
model.state(2).transit(1).rate(1) = a(2);
model.state(2).transit(1).rate(2) = b(2);

% initial conditions

model.initCond.len = model.lenOffset + ...
    -(v(1)*a(2) + v(2)*a(1)) / (v(1)*b(2) + v(2)*b(1));

model.initCond.state = 1;

% Seed random number generator
rand('twister', 100);

% Markov Chain simulation
[traj, state, model, flag] = mtMarkovLengthRate(model);

% Sample the trajectory at 1 sec interval
tr = sampleTraj(traj, 1);

plot(tr(:, 1), tr(:, 2));

% Fit the trajectory to ARMA model
ref = armaxFitKalman(tr(1:2000, 2), [], [0 2; 0 2; -1 -1], 'nl');

save('ref', 'ref');