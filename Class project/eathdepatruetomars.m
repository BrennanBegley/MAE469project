function [earthDepart_r, earthDepart_v, marsDepart_r, marsDepart_v, ...
          marsArrival_r, marsArrival_v, transferOrbitStruct, ...
          earthDepartureDeltaV, marsTurnAngle, marsVinfOut, marsVinfIn, ...
          marsFlybyDepartureScVelocity, transferScV1, transferScV2] = ...
          eathdepatruetomars(SUNAUtokm, SunStruct, shortLong, ...
          earthJ2000_r, earthJ2000_v, marsJ2000_r, marsJ2000_v, ...
          departDateJulian, transferOrbit1TOF, earthStruct, marsStruct, ...
          earthCircularOrbitAlt, sunMuCanonical, AUTUtoKms, ...
          j2000DateJulian, hyperbolicRadius)
%EATHDEPATRUETOMARS Build the first Earth-to-Mars leg of the mission.
%   This helper computes the heliocentric states of Earth and Mars at the
%   selected departure date, solves Lambert/Gauss for the transfer, then
%   computes Earth departure delta-V and the Mars flyby turn.
%
% Notes:
%   - The original filename/function name is preserved for compatibility.
%   - Heliocentric quantities are mostly in AU and AU/TU.
%   - Planet-centered flyby quantities are converted to km and km/s.

    % Convert departure epoch from Julian date to canonical solar TU.
    departTOF = (departDateJulian - j2000DateJulian) / 58.13;

    % Propagate Earth and Mars from J2000 to the chosen departure date.
    [earthDepart_r, earthDepart_v] = universalTOF(sunMuCanonical, departTOF, earthJ2000_r, earthJ2000_v);
    [marsDepart_r,  marsDepart_v]  = universalTOF(sunMuCanonical, departTOF, marsJ2000_r,  marsJ2000_v);

    % Propagate Mars further to the spacecraft arrival date.
    [marsArrival_r, marsArrival_v] = universalTOF(sunMuCanonical, transferOrbit1TOF, marsDepart_r, marsDepart_v);

    % Solve Gauss/Lambert for the transfer trajectory from Earth departure
    % position to Mars arrival position over the selected time of flight.
    [transferScV1, transferScV2, transferOrbitStruct, ~, ~] = ...
        gaussspeedsandtransferorbitorbitalelements(earthDepart_r, marsArrival_r, transferOrbit1TOF, sunMuCanonical, shortLong);

    % Hyperbolic excess velocity at Earth departure, in km/s.
    earthDepartureVinf = (transferScV1 - earthDepart_v) * AUTUtoKms;

    % Burn needed to leave Earth circular parking orbit onto the transfer.
    earthDepartureDeltaV = deltavforcirculartohyperbolic(earthCircularOrbitAlt, earthStruct, earthDepartureVinf);

    % Compute the Mars flyby geometry and outgoing heliocentric velocity.
    [marsTurnAngle, marsVinfOut, marsVinfIn, marsFlybyDepartureScVelocity] = ...
        hyperbolicturnangle(transferScV2, marsArrival_v, AUTUtoKms, ...
        hyperbolicRadius, marsStruct, SunStruct, SUNAUtokm);
end
