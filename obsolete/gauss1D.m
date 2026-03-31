function [x,w] = gauss1D(n)
% Gauss quadrature nodes and weights on [0,1]
switch n
    case 1
        x = 0.5;   w = 1;
    case 2
        x = [0.211324865405187, 0.788675134594813]';
        w = [0.5, 0.5]';
    case 3
        x = [0.112701665379258, 0.5, 0.887298334620742]';
        w = [5/18, 8/18, 5/18]';
    otherwise
        error('nq > 3 not implemented');
end
end