function [meanRow, covariance] = discrete_posterior_moments(states, probabilities)
%DISCRETE_POSTERIOR_MOMENTS Moments of a joint complex symbol posterior.
% states: number-of-states by number-of-users; probabilities sum to one.
% C(i,j) = E[(x_i-E[x_i])*conj(x_j-E[x_j])], as in paper (21)/(26).
probabilities = probabilities(:);
meanRow = probabilities.'*states;
centered = states-meanRow;
covariance = centered.'*(probabilities.*conj(centered));
end
