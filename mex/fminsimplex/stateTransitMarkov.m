function dS = stateTransitMarkov(state)
% STATETRANSITMARKOV performs a one-step transition in Markov chain simulation.
% 
% SYNOPSIS dS = stateTransitMarkov(state)
%
% INPUT
%       state           : Structure comprising the following fields
%           Public fields
%           .name       : Name of the state
%           .speed      : Mean speed [um/minute]
%                         positive: growth, negative: shrinkage
%           .nTransitions : Number of legal transitions
%           .transit    : Array of transitions.
%               .dS     : Integer j denoting a destination state, model(j)
%
%           Private field
%           .sS         : Integer representing the source state
%           .cdf        : Cumulative distribution of the destination states
%                         on [0, 1]
%
% OUTPUT
%       dS              : Integer representing the destination state
%
% Pei-hsin Hsu
% April 20 2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the selection of destination state dS   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = find(rand - state.cdf < 0, 1, 'first');

if i == state.nTransits + 1      % Stay in the source state
    dS = state.sS;
else                             % Transit to the destination state
    dS = state.transit(i).dS;
end


