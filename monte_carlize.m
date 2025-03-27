if numel(MC.q)==2
    POLAR.q = randn(4,1); POLAR.q = POLAR.q/norm(POLAR.q)*sign(POLAR.q(4));
    FC.quat = POLAR.q+randn(4,1)*1e-2; FC.quat = FC.quat/norm(FC.quat);
end

if numel(MC.w)==2
    POLAR.w = MC.w(1) + rand(3,1)*diff(MC.w);
    FC.w = POLAR.omega+POLAR.bias;
end

if numel(MC.J)==2
    temp = POLAR.J - diag(diag(POLAR.J));
    while true
        D = rand(3,1);
        if all( ((ones(3)-2*eye(3))*D) > 0)
            break % This needs to be true for a moment of inertia matrix
        end
    end
    POLAR.J = temp + eye(3)*MC.J(1) + diff(MC.J)*diag(D);
    FC.J = POLAR.J;
end

