function plotSinglePlanetOneOrbit(r0, v0, Planet_TOF_res, muSunCanonical, planetName)
% Plot one colored planet marker + exactly one full orbit

    % Planet marker color map
    switch lower(planetName)
        case 'mercury'
            planetColor = [0.7 0.7 0.7];   % gray
        case 'venus'
            planetColor = [1 0.8 0.2];     % yellow-orange
        case 'earth'
            planetColor = [0.2 0.5 1];     % blue
        case 'mars'
            planetColor = [1 0.3 0.1];     % red-orange
        case 'jupiter'
            planetColor = [0.9 0.7 0.5];   % tan
        case 'saturn'
            planetColor = [0.95 0.85 0.55]; % pale gold
        case 'uranus'
            planetColor = [0.4 1 1];       % cyan
        case 'neptune'
            planetColor = [0.2 0.4 1];     % deep blue
        case 'pluto'
            planetColor = [0.8 0.6 0.6];   % dusty rose
        otherwise
            planetColor = [1 1 1];
    end

    % Planet marker
    plot3(r0(1), r0(2), r0(3), 'o', ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', planetColor, ...
        'MarkerEdgeColor', 'w', ...
        'DisplayName', planetName)

    % Get orbital elements
    [a, ~, eMag, ~, ~, ~, ~] = orbitalelementscalc(r0, v0, muSunCanonical);

    if isfinite(a) && a > 0 && eMag < 1
        orbitalPeriod_TU = 2*pi*sqrt(a^3 / muSunCanonical);
    else
        error('%s does not have a closed elliptical orbit.', planetName);
    end

    TOF = linspace(0, orbitalPeriod_TU, Planet_TOF_res);

    orbit_r = zeros(3, length(TOF));

    for k = 1:length(TOF)
        [orbit_r(:,k), ~] = universalTOF(muSunCanonical, TOF(k), r0, v0);
    end

    % Keep orbit traces white
    plot3(orbit_r(1,:), orbit_r(2,:), orbit_r(3,:), ...
        '-w', 'LineWidth', 1.0, ...
        'HandleVisibility', 'off')
end