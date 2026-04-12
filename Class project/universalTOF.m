% Use universal TOF with S and C variation to calculate the new R and V at
% some time after the j2000 date with mu ,desired TOF, R vec, V vec 
function [runiIJK,vuniIJK] = universalTOF(mu, dt, r0_vec, v0_vec)
r0=norm(r0_vec);
v0=norm(v0_vec);
energy = v0^2/2 - mu/r0;
a = -mu/(2*energy);
e_vec = ((norm(v0_vec)^2 - mu/r0).*r0_vec - dot(r0_vec,v0_vec).*v0_vec)/mu;
ecc = norm(e_vec);
% Handle no-solution case
if energy >= 0 && abs(dt) > 1e5
    fprintf("Energy >= 0 → Hyperbolic/Parabolic orbit.\n");
    fprintf("Long time-of-flight may not converge.\n\n");
    return;
end

x = sqrt(mu)*dt/a;
tol = 1e-10;
err = 1;
i = 0;

while abs(err) > tol
    z = x^2 / a;
    if ecc < 1   % Elliptical
        if z == 0
            C = 1/2; S = 1/6;
        else
            sqrt_z = sqrt(z);
            S = (sqrt_z - sin(sqrt_z)) / (sqrt_z^3);
            C = (1 - cos(sqrt_z)) / z;
        end

    elseif ecc > 1   % Hyperbolic
        sqrt_neg_z = sqrt(-z);
        S = (sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z^3);
        C = (1 - cosh(sqrt_neg_z)) / z;

    else   % Parabolic (series)
        C = 1/2 - z/24 + z^2/720;
        S = 1/6 - z/120 + z^2/5040;
    end

    % Time equation
    t = x^3*S + (dot(r0_vec,v0_vec)/sqrt(mu))*x^2*C + r0*x*(1 - z*S);

    % Radius equation
    r = x^2*C + (dot(r0_vec,v0_vec)/sqrt(mu))*x*(1 - z*S) + r0*(1 - z*C);

    dtdx = r / sqrt(mu);

    x_new = x + (dt - t)/dtdx;
    err = x_new - x;
    x = x_new;
    i = i + 1;
end

f = 1 - (x^2/r0)*C;
g = t - (x^3/sqrt(mu))*S;

runiIJK = f*r0_vec + g*v0_vec;


fdot = sqrt(mu)*x/(r0*r)*(z*S - 1);
gdot = 1 - (x^2/r)*C;

vuniIJK = fdot*r0_vec + gdot*v0_vec;

if abs(f*gdot - fdot*g - 1) > 1e-8
    fprintf('Universal TOF error f g check not correct')
end
end