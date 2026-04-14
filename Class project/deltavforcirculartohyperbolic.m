function deltaV = deltavforcirculartohyperbolic(circularOrbitAlt_km, planetStruct, hyperbolicExcessVelocity)
%DELTAVFORCIRCULARTOHYPERBOLIC Compute burn from circular orbit to escape hyperbola.
%   deltaV = deltavforcirculartohyperbolic(circularOrbitAlt_km, planetStruct, Vinf)
%   returns the required impulsive burn [km/s] to leave a circular parking
%   orbit and enter a hyperbolic trajectory with a specified V-infinity.
%
% Inputs:
%   circularOrbitAlt_km      Parking-orbit altitude above the planet surface [km]
%   planetStruct             Structure with fields:
%                               mu : gravitational parameter [km^3/s^2]
%                               r  : planet radius [km]
%   hyperbolicExcessVelocity Hyperbolic excess velocity vector or scalar [km/s]
%
% Output:
%   deltaV                   Required burn magnitude [km/s]

    % If the caller provides a vector, only its magnitude matters here.
    vinf = norm(hyperbolicExcessVelocity);

    % Radius of the circular parking orbit measured from the planet center.
    parkingRadius_km = planetStruct.r + circularOrbitAlt_km;

    % Specific orbital energy of the desired departure hyperbola.
    hyperbolicSpecificEnergy = vinf^2 / 2;  % [km^2/s^2]

    % Speed required at periapsis of the hyperbola.
    hyperbolicSpeedAtPeriapsis = sqrt(2 * (hyperbolicSpecificEnergy + planetStruct.mu / parkingRadius_km));

    % Speed of the initial circular parking orbit.
    circularOrbitSpeed = sqrt(planetStruct.mu / parkingRadius_km);

    % Impulsive burn needed to transition from circular to hyperbolic motion.
    deltaV = hyperbolicSpeedAtPeriapsis - circularOrbitSpeed;
end
