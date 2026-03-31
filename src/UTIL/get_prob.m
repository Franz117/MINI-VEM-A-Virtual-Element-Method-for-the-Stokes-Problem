function sol = get_prob(problem_name)

if nargin < 1
    problem_name = 'default';
end


switch problem_name
    case 'compute'
        syms x y
        mu = 1;
        u1_sym = -cos(2*pi*x)*sin(2*pi*y) + sin(2*pi*y);
        u2_sym =  sin(2*pi*x)*cos(2*pi*y) - sin(2*pi*x);
        p_sym  =  2*pi*(cos(2*pi*x) - cos(2*pi*y));

        du1dx_sym = diff(u1_sym,x); du1dy_sym = diff(u1_sym,y);
        du2dx_sym = diff(u2_sym,x); du2dy_sym = diff(u2_sym,y);

        f1_sym = -mu*(diff(u1_sym,x,2) + diff(u1_sym,y,2)) + diff(p_sym,x);
        f2_sym = -mu*(diff(u2_sym,x,2) + diff(u2_sym,y,2)) + diff(p_sym,y);

        % Convert to function handles
        sol.u1    = matlabFunction(u1_sym,'Vars',[x y]);
        sol.u2    = matlabFunction(u2_sym,'Vars',[x y]);
        sol.p     = matlabFunction(p_sym,'Vars',[x y]);
        sol.f1    = matlabFunction(f1_sym,'Vars',[x y]);
        sol.f2    = matlabFunction(f2_sym,'Vars',[x y]);
        sol.du1dx = matlabFunction(du1dx_sym,'Vars',[x y]);
        sol.du1dy = matlabFunction(du1dy_sym,'Vars',[x y]);
        sol.du2dx = matlabFunction(du2dx_sym,'Vars',[x y]);
        sol.du2dy = matlabFunction(du2dy_sym,'Vars',[x y]);

    case 'default'
        sol.u1 = @(x,y)sin(y.*pi.*2.0)-cos(x.*pi.*2.0).*sin(y.*pi.*2.0);
        sol.u2 = @(x,y)-sin(x.*pi.*2.0)+cos(y.*pi.*2.0).*sin(x.*pi.*2.0);
        sol.p  = @(x,y) -2.0*pi.*(cos(x.*pi.*2.0)-cos(y.*pi.*2.0));
        sol.f1 = @(x,y)pi.^2.*(sin(x.*pi.*2.0)+sin(y.*pi.*2.0)-cos(x.*pi.*2.0).*sin(y.*pi.*2.0).*2.0).*4.0;
        sol.f2 = @(x,y)pi.^2.*(sin(x.*pi.*2.0)+sin(y.*pi.*2.0)-cos(y.*pi.*2.0).*sin(x.*pi.*2.0).*2.0).*-4.0;
        sol.du1dx = @(x,y)pi.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*2.0;
        sol.du1dy = @(x,y)pi.*cos(y.*pi.*2.0).*(cos(x.*pi.*2.0)-1.0).*-2.0;
        sol.du2dx = @(x,y)pi.*cos(x.*pi.*2.0).*(cos(y.*pi.*2.0)-1.0).*2.0;
        sol.du2dy = @(x,y)pi.*sin(x.*pi.*2.0).*sin(y.*pi.*2.0).*-2.0;

    case 'patch_test'
        sol.u1    = @(x,y) 1 + 2*x + 3*y;
        sol.u2    = @(x,y) -1 + x - 2*y;
        sol.p     = @(x,y) x - y;
        sol.f1    = @(x,y) 1 +0*x + 0*y;
        sol.f2    = @(x,y) -1 + 0*x + 0*y;
        sol.du1dx = @(x,y) 2+ 0*x + 0*y;
        sol.du1dy = @(x,y) 3 + 0*x + 0*y;
        sol.du2dx = @(x,y) 1+ 0*x + 0*y;
        sol.du2dy = @(x,y) -2+ 0*x + 0*y;

    otherwise
        error('Unknown problem_name: %s', problem_name);
end