function fval = myfun(x, model)
% MYFUN is a user-defined objective function for mexFuntion fminsimplex. 
% In this example, function evaluation is based on the linear combination 
% of - log (p value) of ARMA coefficients and WN variance comparison with
% the reference ARMA(2, 1) model.
%
% Pei-hsin Hsu
% April 28, 2009

fprintf('[%f, %f]: ', x(1), x(2));

model.state(1).transit(1).rate(2) = x(1) * model.scale(1);
model.state(2).transit(1).rate(1) = x(2) * model.scale(2);

if x(1) < 0 || x(2) < 0
    fval = 1.0E+3;
    fprintf('out of bounds\n');
    fprintf('fval: %f\n\n', fval);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('twister', 100);

% Markov Chain simulation
[traj, state, model] = mtMarkovLengthRate(model);

% Sample trajectory at 1 sec interval
tr = sampleTraj(traj, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit sampled trajectory to ARMA(2, 1) model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit = armaxFitKalman(tr(1:2000, 2), [], [2 2; 1 1; -1 -1], 'nl');

% Compare the fit model with the reference
[mLogPCoef, mLogPVar] = armaxModelComp(fit, model.ref);

fval = mLogPCoef + (5/12) * mLogPVar;

fprintf('fval: %f\n\n', fval);
