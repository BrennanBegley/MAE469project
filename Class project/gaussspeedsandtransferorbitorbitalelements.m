function [transferScV1, transferScV2, transferOrbitStruct, transferOrbit_r1, transferOrbit_r2] = ...
    gaussspeedsandtransferorbitorbitalelements(departurePos_r1, arrivalPos_r2, TOF, sunMu, shortLong)
%GAUSSSPEEDSANDTRANSFERORBITORBITALELEMENTS Solve transfer and package results.
%   This is a convenience wrapper around Gauss() and orbitalelementscalc().
%   It returns transfer endpoint velocities plus a structure containing the
%   orbital elements of the transfer orbit.

    transferOrbit_r1 = departurePos_r1;
    transferOrbit_r2 = arrivalPos_r2;

    % Solve for transfer velocities.
    [transferScV1, transferScV2] = Gauss(departurePos_r1, arrivalPos_r2, TOF, shortLong, sunMu);

    % Compute orbital elements from the departure state on the transfer orbit.
    [a, ~, eMag, inc, RAAN, argPeri, trueAnom] = orbitalelementscalc(departurePos_r1, transferScV1, sunMu);

    transferOrbitStruct = struct( ...
        'a',     a, ...
        'e',     eMag, ...
        'inc',   inc, ...
        'OMEGA', RAAN, ...
        'omega', argPeri, ...
        'theta', trueAnom);
end
