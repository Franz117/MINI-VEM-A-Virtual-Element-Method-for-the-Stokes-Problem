function [u1,u2,p,f1,f2,mu] = get_prob()

    syms x y

    % Define u1 and u2
    u1_sym = -cos(pi*2*x)*sin(pi*2*y) + sin(pi*2*y);
    u2_sym =  sin(pi*2*x)*cos(pi*2*y) - sin(pi*2*x);
    p_sym = cos(pi*2*x)*sin(pi*2*y);
    mu = 1;

    % Compute -Laplacian(u1) and -Laplacian(u2)
    f1_sym = -mu*diff(u1_sym, x, 2) - mu*diff(u1_sym, y, 2) +diff(p_sym, x, 1);
    f2_sym = -mu*diff(u2_sym, x, 2) - mu*diff(u2_sym, y, 2) +diff(p_sym, y, 1);

    % Convert to MATLAB function handles
    u1 = matlabFunction(u1_sym, 'Vars', [x y]);
    u2 = matlabFunction(u2_sym, 'Vars', [x y]);
    p = matlabFunction(p_sym, 'Vars', [x y]);
    f1 = matlabFunction(f1_sym, 'Vars', [x y]);
    f2 = matlabFunction(f2_sym, 'Vars', [x y]);


end