function [r_ijk, v_ijk] = posandvelvector(planet, mu)
%POSANDVELVECTOR Convert orbital elements to inertial state vectors.
%   [r_ijk, v_ijk] = posandvelvector(planet, mu) builds the position and
%   velocity vectors from the classical orbital elements stored in planet.
%
% Required fields in planet:
%   a, e, inc, OMEGA, omega, theta
%
% Angles are expected in radians.

    a       = planet.a;
    e       = planet.e;
    inc     = planet.inc;
    OMEGA   = planet.OMEGA;
    omega   = planet.omega;
    trueAnom = planet.theta;

    % Semilatus rectum and perifocal-frame state.
    p = a * (1 - e^2);
    r = p / (1 + e*cos(trueAnom));

    r_pqw = [r*cos(trueAnom); r*sin(trueAnom); 0];
    v_pqw = sqrt(mu/p) * [-sin(trueAnom); e + cos(trueAnom); 0];

    % Rotation from perifocal coordinates to inertial IJK.
    R_pqw_to_ijk = [ ...
        cos(OMEGA)*cos(omega) - sin(OMEGA)*sin(omega)*cos(inc), ...
       -cos(OMEGA)*sin(omega) - sin(OMEGA)*cos(omega)*cos(inc), ...
        sin(OMEGA)*sin(inc); ...
        sin(OMEGA)*cos(omega) + cos(OMEGA)*sin(omega)*cos(inc), ...
       -sin(OMEGA)*sin(omega) + cos(OMEGA)*cos(omega)*cos(inc), ...
       -cos(OMEGA)*sin(inc); ...
        sin(omega)*sin(inc), cos(omega)*sin(inc), cos(inc)];

    % Inertial state vector.
    r_ijk = R_pqw_to_ijk * r_pqw;
    v_ijk = R_pqw_to_ijk * v_pqw;
end
