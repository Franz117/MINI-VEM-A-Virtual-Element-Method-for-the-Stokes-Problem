function [u1, u2, p, f1, f2, mu] = get_prob2()

    syms x y

    % Very simple linear fields
    u1_sym = x;
    u2_sym = y;
    p_sym  = x + y;
    mu = 1;

    % Compute -μΔu + ∇p (Stokes RHS)
    f1_sym = -mu * (diff(u1_sym, x, 2) + diff(u1_sym, y, 2)) + diff(p_sym, x);
    f2_sym = -mu * (diff(u2_sym, x, 2) + diff(u2_sym, y, 2)) + diff(p_sym, y);

    % Convert to functions
    u1 = matlabFunction(u1_sym, 'Vars', [x y]);
    u2 = matlabFunction(u2_sym, 'Vars', [x y]);
    p  = matlabFunction(p_sym,  'Vars', [x y]);
    f1 = matlabFunction(f1_sym, 'Vars', [x y]);
    f2 = matlabFunction(f2_sym, 'Vars', [x y]);

end