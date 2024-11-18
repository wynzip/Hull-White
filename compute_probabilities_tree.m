function [p_up, p_m, p_down] = compute_probabilities_tree(l_max, mu_hat)
% This function computes the probabilities og going up, down or staying
% level in a trinomial tree, following H-W model.
%
% INPUT
% l_max:        maximum number of vertical steps
% mu_hat:       parameter of tree hull white
%
% OUTPUT
% p_up:         vector of probabilities of going upwards
% p_m:          vector of probabilities of staying level
% p_down:       vector of probabilities of going downwards
%

l_min = -l_max;

% initialize vectors
p_up = zeros(l_max*2+1, 1);
p_m = zeros(l_max*2+1, 1);
p_down = zeros(l_max*2+1, 1);

% particular scenario corresponding to l_max: the tree doesn't allow to go
% upwards anymore.
p_up(1) = 0.5*(7/3 - 3*l_max*mu_hat + (l_max*mu_hat)^2);
p_m(1) = -1/3 + 2*l_max*mu_hat - (l_max*mu_hat)^2;
p_down(1) = 0.5*(1/3 - l_max*mu_hat + (l_max*mu_hat)^2);

% particular scenario corresponding to l_min: the tree doesn't allow to go
% downwards anymore.
p_up(2*l_max + 1) = 0.5*(1/3 + l_min*mu_hat + (l_min*mu_hat)^2);
p_m(2*l_max + 1) = -1/3 - 2*l_min*mu_hat - (l_min*mu_hat)^2;
p_down(2*l_max + 1) = 0.5*(7/3 + 3*l_min*mu_hat + (l_min*mu_hat)^2);

l_vect = [l_max:-1:-l_max];

%regular scenarios
for i = 2: 2*l_max
    l = l_vect(i);
    p_up(i) = 0.5 * (1/3 - l*mu_hat + (l*mu_hat)^2);
    p_m(i) = 2/3 - (l*mu_hat)^2;
    p_down(i) = 0.5 * (1/3 + l*mu_hat + (l*mu_hat)^2);
end


end