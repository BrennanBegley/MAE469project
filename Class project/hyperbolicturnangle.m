function [turnAngle, vInfOut_IJK, vInfIn_AU_TU, scDepartureHeliocentric_v] = ...
    hyperbolicturnangle(scArrival_v, planetArrival_v, AUTUtoKms, flybyAltitude_km, planetStruct, SunStruct, SUNAUtokm)
%HYPERBOLICTURNANGLE Compute a ballistic flyby turn and outgoing velocity.
%   This function estimates the turn angle for an unpowered planetary flyby,
%   rotates the incoming planet-relative hyperbolic excess velocity in the
%   perifocal frame, and converts the result back to heliocentric velocity.
%
% Inputs:
%   scArrival_v            Spacecraft heliocentric arrival velocity [AU/TU]
%   planetArrival_v        Planet heliocentric arrival velocity [AU/TU]
%   AUTUtoKms              AU/TU to km/s conversion factor
%   flybyAltitude_km       Flyby periapsis altitude above the planet [km]
%   planetStruct           Planet data structure
%   SunStruct              Sun data structure
%   SUNAUtokm              AU to km conversion factor

    % Sphere of influence radius in AU.
    planetSOI_AU = norm(planetStruct.a) * (planetStruct.mass / SunStruct.mass)^(2/5);

    % Approximate incoming asymptote position direction in the planet frame.
    scSOIarrival_r = -planetSOI_AU * scArrival_v / norm(scArrival_v);

    % Incoming hyperbolic excess velocity in the planet-centered frame.
    vInfIn_AU_TU = scArrival_v - planetArrival_v;

    % Build the incoming hyperbolic orbit using planet-centered units.
    [a, ~, eMag, inc, RAAN, argPeri, trueAnom] = ...
        orbitalelementscalc(scSOIarrival_r * SUNAUtokm, vInfIn_AU_TU * AUTUtoKms, planetStruct.mu);

    % Semilatus rectum and instantaneous PQW state for the incoming branch.
    p_sc = a * (1 - eMag^2);      % [km]
    r_pqw_mag = p_sc / (1 + eMag*cos(trueAnom));
    r_pqw = [r_pqw_mag*cos(trueAnom); r_pqw_mag*sin(trueAnom); 0];

    % Velocity in the perifocal frame.
    P_term = [-sin(trueAnom); 0; 0];
    Q_term = [0; eMag + cos(trueAnom); 0];
    vInf_pqw_in = sqrt(planetStruct.mu / p_sc) * (P_term + Q_term); % [km/s]

    % Ballistic flyby turn angle from hyperbolic geometry.
    periapsisRadius_km = planetStruct.r + flybyAltitude_km;
    turnAngle = 2 * asin(1 / (1 + periapsisRadius_km * norm(vInfIn_AU_TU * AUTUtoKms)^2 / planetStruct.mu));

    % Rotate the asymptotic velocity vector by the flyby turn angle.
    H = [cos(turnAngle), -sin(turnAngle), 0; ...
         sin(turnAngle),  cos(turnAngle), 0; ...
         0,               0,              1];
    vInf_pqw_out = H * vInf_pqw_in;

    % PQW -> IJK rotation matrix.
    R_pqw_to_ijk = [ ...
        cos(RAAN)*cos(argPeri)-sin(RAAN)*sin(argPeri)*cos(inc), ...
       -cos(RAAN)*sin(argPeri)-sin(RAAN)*cos(argPeri)*cos(inc), ...
        sin(RAAN)*sin(inc); ...
        sin(RAAN)*cos(argPeri)+cos(RAAN)*sin(argPeri)*cos(inc), ...
       -sin(RAAN)*sin(argPeri)+cos(RAAN)*cos(argPeri)*cos(inc), ...
       -cos(RAAN)*sin(inc); ...
        sin(argPeri)*sin(inc), cos(argPeri)*sin(inc), cos(inc)];

    % Outgoing asymptotic velocity in the inertial frame [km/s].
    vInfOut_IJK = R_pqw_to_ijk * vInf_pqw_out;

    % Convert back to heliocentric spacecraft velocity [AU/TU].
    scDepartureHeliocentric_v = (vInfOut_IJK + planetArrival_v * AUTUtoKms) / AUTUtoKms;

    % Helpful diagnostic output for debugging/reporting.
    fprintf('Flyby semilatus rectum p [km]: %g\n', p_sc);
    fprintf('Incoming PQW position [km]: [%g, %g, %g]\n', r_pqw);
    fprintf('Incoming PQW velocity [km/s]: [%g, %g, %g]\n', vInf_pqw_in);
    fprintf('Outgoing V-infinity in IJK [km/s]: [%g, %g, %g]\n', vInfOut_IJK);
end
