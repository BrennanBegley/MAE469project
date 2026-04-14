function [r_ijk_new, v_ijk_new] = universalTOF(mu, dt, r0_vec, v0_vec)
%UNIVERSALTOF Propagate a Keplerian orbit using the universal-variable form.
%   [r_new, v_new] = universalTOF(mu, dt, r0_vec, v0_vec) propagates the
%   initial state (r0_vec, v0_vec) by a time dt under two-body dynamics.
%
% Inputs:
%   mu      Gravitational parameter
%   dt      Time of flight
%   r0_vec  Initial position vector
%   v0_vec  Initial velocity vector

    r0 = norm(r0_vec);
    v0 = norm(v0_vec);

    specificEnergy = v0^2/2 - mu/r0;
    a = -mu / (2*specificEnergy);

    e_vec = ((norm(v0_vec)^2 - mu/r0) .* r0_vec - dot(r0_vec, v0_vec) .* v0_vec) / mu;
    ecc = norm(e_vec);

    % Guard against cases that may struggle to converge with this setup.
    if specificEnergy >= 0 && abs(dt) > 1e5
        fprintf('Energy >= 0 -> Hyperbolic/Parabolic orbit.\n');
        fprintf('Long time-of-flight may not converge.\n\n');
        return;
    end

    % Initial guess for the universal anomaly x.
    x = sqrt(mu) * dt / a;
    tol = 1e-10;
    err = 1;

    while abs(err) > tol
        z = x^2 / a;

        % Stumpff functions.
        if ecc < 1
            % Elliptic case.
            if z == 0
                C = 1/2;
                S = 1/6;
            else
                sqrt_z = sqrt(z);
                S = (sqrt_z - sin(sqrt_z)) / (sqrt_z^3);
                C = (1 - cos(sqrt_z)) / z;
            end
        elseif ecc > 1
            % Hyperbolic case.
            sqrt_neg_z = sqrt(-z);
            S = (sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z^3);
            C = (1 - cosh(sqrt_neg_z)) / z;
        else
            % Parabolic limit using a series expansion.
            C = 1/2 - z/24 + z^2/720;
            S = 1/6 - z/120 + z^2/5040;
        end

        % Universal Kepler equation.
        t = x^3*S + (dot(r0_vec, v0_vec)/sqrt(mu))*x^2*C + r0*x*(1 - z*S);

        % Radius magnitude at the propagated time.
        r = x^2*C + (dot(r0_vec, v0_vec)/sqrt(mu))*x*(1 - z*S) + r0*(1 - z*C);

        % Newton update.
        dtdx = r / sqrt(mu);
        x_new = x + (dt - t) / dtdx;
        err = x_new - x;
        x = x_new;
    end

    % Lagrange f and g functions.
    f = 1 - (x^2/r0) * C;
    g = t - (x^3/sqrt(mu)) * S;

    % Propagated position.
    r_ijk_new = f*r0_vec + g*v0_vec;

    % Lagrange derivatives for the propagated velocity.
    fdot = sqrt(mu)*x/(r0*r) * (z*S - 1);
    gdot = 1 - (x^2/r) * C;
    v_ijk_new = fdot*r0_vec + gdot*v0_vec;

    % Sanity check for f-g consistency.
    if abs(f*gdot - fdot*g - 1) > 1e-5
        fprintf('Universal TOF error: f-g check is not correct.\n')
        disp(f*gdot - fdot*g - 1)
    end
end
