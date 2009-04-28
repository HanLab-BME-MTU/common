function [traj, stateIndex, model, errFlag] = mtMarkovLengthRate(model)

% MTMARKOVLENGTHRATE uses Markov Chain simulation to generate a MT trajectory.
%
% SYNOPSIS [traj, stateIndex, model, errFlag] = mtMarkovLengthRate(model)
%
% INPUT
%   
%   model
%
%       .nStates        : Number of states
%       .state          : Array of states. A state contains:
%           .name       : Name of the state
%           .speed      : Mean speed [um/minute]
%                         positive: growth, negative: shrinkage
%           .nTransits  : Number of legal transitions from the source state
%           .transit    : Array of transitions. A transition contains:
%               .rate(1): Transition frequency [1/s]
%               .rate(2): Slope [1/(s * um)], optional
%               .dS     : Integer representing the destination state
%       .lenOffset      : Offset length [um]
%       .simParam
%           .time       : Total simulation time [s]
%           .dT         : Time step [s]
%       .initCond       : Initial conditions
%           .len        : Initial length [um]
%           .state      : Initial state
%
% OUTPUT
%       traj            : Array comprising two columns 
%                         1st column: Simulation time [s]
%                         2nd column: MT length [um]
%       stateIndex      : State tranistion record (for debugging)
%       model           : original model with private fields added
%
% Recent modification:
% (line 121) if dist <= 0 && model.state(stateIndex(i-1)).name == 'g'
% The <= has replaced original < to take into account the machine epsilon.
%
% Pei-hsin Hsu
% April 23, 2009

time = model.simParam.time;
dT = model.simParam.dT;

traj = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 1
    fprintf('--mtMarkov: argument length is wrong.\n');
    errFlag = 1;
    return;
end

% Assign slope of transition rate to 0 if absent from input
for i = 1:model.nStates
    for j = 1:model.state(i).nTransits
        rate = model.state(i).transit(j).rate;
        if length(rate) == 1; model.state(i).transit(j).rate = [rate, 0]; end
    end
end


% Initialize state private field
% ------------------------------
% The following field of a state is on the implementation side. User should not
% have access to it.
%
% sS                    : Integer representing the source state
%                         Note: The destination state, dS, is on the client 
%                         side.

for i = 1:model.nStates; model.state(i).sS = i; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTimeSteps = ceil(time/dT); % total number of time steps
traj = zeros(nTimeSteps, 2);
traj(1, :) = [0, model.initCond.len];
stateIndex = zeros(nTimeSteps, 1);
r = zeros(nTimeSteps, 1);
stateIndex(1) = model.initCond.state;
r(1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Markov chain simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:nTimeSteps
    
    % Update the output trajectory
    traj(i, 1) = traj(i-1, 1) + dT;
    traj(i, 2) = traj(i-1, 2) + model.state(stateIndex(i-1)).speed * (dT/60);
    
    % Distance from current position to the length offset
    dist = traj(i, 2) - model.lenOffset;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the cumulative distribution function
    % of destination states on [0, 1]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nTransits = model.state(stateIndex(i-1)).nTransits;
    
    c = zeros(1, nTransits + 1);
    
    for j = 1:nTransits

        rate = model.state(stateIndex(i-1)).transit(j).rate(1) + ...
               model.state(stateIndex(i-1)).transit(j).rate(2) * dist;

        c(j) = rate * dT;
        
        if dist <= 0 && model.state(stateIndex(i-1)).name == 'g'
        %if dist < 0 && model.state(stateIndex(i-1)).name == 'g'
            c(1) = 0;
        end
        
        if c(j) < 0; c(j) = 0; end
        
        if j > 1; c(j) = c(j) + c(j-1); end

        if c(j) > 1; c(j) = 1; end
    end
    
    c(nTransits + 1) = 1; % upper bound of the source state
    
    model.state(stateIndex(i-1)).cdf = c;
    
    % Perform a one-step state transition
    stateIndex(i) = stateTransitMarkov(model.state(stateIndex(i-1)));
    
end

errFlag = 0;