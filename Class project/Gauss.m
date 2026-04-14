function [v1, v2] = Gauss(r1, r2, TOF, short, mu)
%GAUSS Solve the two-point boundary-value transfer using Gauss/Lambert form.
%   [v1, v2] = Gauss(r1, r2, TOF, short, mu) returns the transfer velocity
%   vectors at r1 and r2 for a specified time of flight.
%
% Inputs:
%   r1, r2   Initial and final position vectors
%   TOF      Time of flight in the same canonical time units as mu
%   short    1 for short-way transfer, otherwise long-way transfer
%   mu       Gravitational parameter

    t = 0;
    tol = 1e-10;

    rdot = dot(r1, r2);
    normr1 = norm(r1);
    normr2 = norm(r2);

    % Transfer angle between the two position vectors.
    if short == 1
        delTheta = acos(rdot / (normr1 * normr2));
    else
        delThetaShort = acos(rdot / (normr1 * normr2));
        delTheta = 2*pi - delThetaShort;
    end

    % Standard Gauss/Lambert auxiliary quantity A.
    A = sqrt(normr1 * normr2) * sin(delTheta) / sqrt(1 - cos(delTheta));

    % Initial guess for the universal variable z.
    z = 1;

    while abs(TOF - t) >= tol
        % Stumpff functions S(z) and C(z).
        if z > 0
            S = (sqrt(z) - sin(sqrt(z))) / sqrt(z^3);
            C = (1 - cos(sqrt(z))) / z;
        elseif z < 0
            % Preserve the original implementation, but present it more clearly.
            S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
            C = (cosh(sqrt(-z)) - 1)/(-z);
        else
            S = 1/factorial(3) - z/factorial(5) + z^2/factorial(7);
            C = 1/factorial(2) - z/factorial(4) + z^2/factorial(6);
        end

        % Intermediate Gauss/Lambert quantities.
        y = normr2 + normr1 - A * (1 - z*S) / sqrt(C);
        x = sqrt(y / C);

        % Time-of-flight function and its derivative.
        t = (x^3 * S + A * sqrt(y)) / sqrt(mu);
        ds = (C - 3*S) / (2*z);
        dc = (1 - z*S - 2*C) / (2*z);
        dt = (x^3 * (ds - (3*S*dc)/(2*C)) + A/8 * ((3*S*sqrt(y))/C + A/x)) / sqrt(mu);

        % Newton iteration on z.
        z = z + (TOF - t) / dt;
    end

    % Lagrange f and g coefficients.
    f    = 1 - (y / normr1);
    g    = A * sqrt(y / mu);
    fdot = (-sqrt(mu) * x) / (normr2 * normr1) * (1 - z*S);
    gdot = 1 - (y / normr2);

    if abs(f*gdot - fdot*g - 1) > 1e-8
        fprintf('Universal TOF error: f and g consistency check failed.\n')
    end

    % Transfer velocities at the boundary points.
    v1 = (r2 - f*r1) / g;
    v2 = (gdot*r2 - r1) / g;
end
