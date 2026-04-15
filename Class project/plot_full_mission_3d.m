function plot_full_mission_3d( ...
    Earth_dep_r, Earth_dep_v, ...
    Mars_dep_r, Mars_dep_v, ...
    Jupiter_dep_r, Jupiter_dep_v, ...
    SC_leg1_r0, SC_leg1_v0, TOF_leg1, ...
    SC_leg2_r0, SC_leg2_v0, TOF_leg2, ...
    SC_leg3_r0, SC_leg3_v0, TOF_leg3, ...
    SC_leg4_r0, SC_leg4_v0, TOF_leg4, ...
    nPoints, figurenumber)
% PLOT_FULL_MISSION_3D
% Plots the full interplanetary mission in 3D and marks planetary positions
% at each major mission event.

    if nargin < 19 || isempty(nPoints)
        nPoints = 500;
    end

    if nargin < 20 || isempty(figurenumber)
        figurenumber = 5;
    end

    sunMu = 1;  % solar canonical units

    % -----------------------------
    % Total mission duration
    % -----------------------------
    totalTOF = TOF_leg1 + TOF_leg2 + TOF_leg3 + TOF_leg4;
    tPlanet = linspace(0, totalTOF, nPoints);

    % -----------------------------
    % Preallocate planet trajectories
    % -----------------------------
    Earth_r   = zeros(3, nPoints);
    Mars_r    = zeros(3, nPoints);
    Jupiter_r = zeros(3, nPoints);

    % -----------------------------
    % Propagate planets over full mission
    % -----------------------------
    for k = 1:nPoints
        [Earth_r(:,k), ~]   = universalTOF(sunMu, tPlanet(k), Earth_dep_r,   Earth_dep_v);
        [Mars_r(:,k), ~]    = universalTOF(sunMu, tPlanet(k), Mars_dep_r,    Mars_dep_v);
        [Jupiter_r(:,k), ~] = universalTOF(sunMu, tPlanet(k), Jupiter_dep_r, Jupiter_dep_v);
    end

    % -----------------------------
    % Build time arrays for each spacecraft leg
    % -----------------------------
    t1 = linspace(0, TOF_leg1, nPoints);
    t2 = linspace(0, TOF_leg2, nPoints);
    t3 = linspace(0, TOF_leg3, nPoints);
    t4 = linspace(0, TOF_leg4, nPoints);

    % -----------------------------
    % Preallocate spacecraft trajectories
    % -----------------------------
    SC1_r = zeros(3, nPoints);
    SC2_r = zeros(3, nPoints);
    SC3_r = zeros(3, nPoints);
    SC4_r = zeros(3, nPoints);

    % -----------------------------
    % Propagate spacecraft legs
    % -----------------------------
    for k = 1:nPoints
        [SC1_r(:,k), ~] = universalTOF(sunMu, t1(k), SC_leg1_r0, SC_leg1_v0);
        [SC2_r(:,k), ~] = universalTOF(sunMu, t2(k), SC_leg2_r0, SC_leg2_v0);
        [SC3_r(:,k), ~] = universalTOF(sunMu, t3(k), SC_leg3_r0, SC_leg3_v0);
        [SC4_r(:,k), ~] = universalTOF(sunMu, t4(k), SC_leg4_r0, SC_leg4_v0);
    end

    % -----------------------------
    % Spacecraft event locations
    % -----------------------------
    Earth_departure = SC_leg1_r0;
    Mars_flyby      = SC1_r(:,end);
    DSM_point       = SC2_r(:,end);
    Earth_flyby     = SC3_r(:,end);
    Jupiter_arrival = SC4_r(:,end);

    % -----------------------------
    % Mission event times
    % -----------------------------
    eventTimes = [ ...
        0, ...
        TOF_leg1, ...
        TOF_leg1 + TOF_leg2, ...
        TOF_leg1 + TOF_leg2 + TOF_leg3, ...
        totalTOF];

    eventNames = { ...
        'Earth Departure', ...
        'Mars Flyby', ...
        'DSM', ...
        'Earth Flyby', ...
        'Jupiter Arrival'};

    % -----------------------------
    % Planet positions at each event
    % -----------------------------
    Earth_event_r   = zeros(3, numel(eventTimes));
    Mars_event_r    = zeros(3, numel(eventTimes));
    Jupiter_event_r = zeros(3, numel(eventTimes));

    for k = 1:numel(eventTimes)
        [Earth_event_r(:,k), ~]   = universalTOF(sunMu, eventTimes(k), Earth_dep_r,   Earth_dep_v);
        [Mars_event_r(:,k), ~]    = universalTOF(sunMu, eventTimes(k), Mars_dep_r,    Mars_dep_v);
        [Jupiter_event_r(:,k), ~] = universalTOF(sunMu, eventTimes(k), Jupiter_dep_r, Jupiter_dep_v);
    end

    % -----------------------------
    % Plot setup
    % -----------------------------
    figure(figurenumber);
    clf;
    hold on;
    grid on;
    axis equal;
    view(3);

    set(gcf, 'Color', [0.06 0.06 0.06]);
    ax = gca;
    ax.Color = [0.03 0.03 0.03];
    ax.XColor = [0.95 0.95 0.95];
    ax.YColor = [0.95 0.95 0.95];
    ax.ZColor = [0.95 0.95 0.95];
    ax.GridColor = [0.5 0.5 0.5];

    % -----------------------------
    % Colors
    % -----------------------------
    earthColor   = [0.2 0.6 1.0];
    marsColor    = [1.0 0.4 0.1];
    jupiterColor = [0.9 0.8 0.6];

    % -----------------------------
    % Plot Sun
    % -----------------------------
    plot3(0, 0, 0, 'yo', ...
        'MarkerFaceColor', 'y', ...
        'MarkerSize', 10, ...
        'DisplayName', 'Sun');

    % -----------------------------
    % Plot planet orbits
    % -----------------------------
    plot3(Earth_r(1,:),   Earth_r(2,:),   Earth_r(3,:),   '-', ...
        'Color', earthColor, 'LineWidth', 1.2, 'DisplayName', 'Earth Orbit');

    plot3(Mars_r(1,:),    Mars_r(2,:),    Mars_r(3,:),    '-', ...
        'Color', marsColor, 'LineWidth', 1.2, 'DisplayName', 'Mars Orbit');

    plot3(Jupiter_r(1,:), Jupiter_r(2,:), Jupiter_r(3,:), '-', ...
        'Color', jupiterColor, 'LineWidth', 1.2, 'DisplayName', 'Jupiter Orbit');

    % -----------------------------
    % Plot spacecraft trajectory legs
    % -----------------------------
    plot3(SC1_r(1,:), SC1_r(2,:), SC1_r(3,:), 'w-', 'LineWidth', 2.0, ...
        'DisplayName', 'Leg 1: Earth to Mars');
    plot3(SC2_r(1,:), SC2_r(2,:), SC2_r(3,:), 'c-', 'LineWidth', 2.0, ...
        'DisplayName', 'Leg 2: Mars to DSM');
    plot3(SC3_r(1,:), SC3_r(2,:), SC3_r(3,:), 'm-', 'LineWidth', 2.0, ...
        'DisplayName', 'Leg 3: DSM to Earth');
    plot3(SC4_r(1,:), SC4_r(2,:), SC4_r(3,:), 'g-', 'LineWidth', 2.0, ...
        'DisplayName', 'Leg 4: Earth to Jupiter');

    % -----------------------------
    % Plot spacecraft major mission events
    % -----------------------------
    plot3(Earth_departure(1), Earth_departure(2), Earth_departure(3), ...
        'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'SC Earth Departure');

    plot3(Mars_flyby(1), Mars_flyby(2), Mars_flyby(3), ...
        'o', 'Color', marsColor, 'MarkerFaceColor', marsColor, ...
        'MarkerSize', 8, 'DisplayName', 'SC Mars Flyby');

    plot3(DSM_point(1), DSM_point(2), DSM_point(3), ...
        'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8, 'DisplayName', 'SC DSM');

    plot3(Earth_flyby(1), Earth_flyby(2), Earth_flyby(3), ...
        'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'DisplayName', 'SC Earth Flyby');

    plot3(Jupiter_arrival(1), Jupiter_arrival(2), Jupiter_arrival(3), ...
        'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'SC Jupiter Arrival');

    % -----------------------------
    % Plot planets at each major event
    % -----------------------------
    for k = 1:numel(eventTimes)
        % Earth
        plot3(Earth_event_r(1,k), Earth_event_r(2,k), Earth_event_r(3,k), ...
            'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', earthColor, ...
            'MarkerEdgeColor', 'w', ...
            'HandleVisibility', 'off');

        % Mars
        plot3(Mars_event_r(1,k), Mars_event_r(2,k), Mars_event_r(3,k), ...
            'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', marsColor, ...
            'MarkerEdgeColor', 'w', ...
            'HandleVisibility', 'off');

        % Jupiter
        plot3(Jupiter_event_r(1,k), Jupiter_event_r(2,k), Jupiter_event_r(3,k), ...
            'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', jupiterColor, ...
            'MarkerEdgeColor', 'w', ...
            'HandleVisibility', 'off');

        % Labels
        text(Earth_event_r(1,k), Earth_event_r(2,k), Earth_event_r(3,k), ...
            sprintf('  E%d', k), 'Color', earthColor, 'FontSize', 8, 'FontWeight', 'bold');

        text(Mars_event_r(1,k), Mars_event_r(2,k), Mars_event_r(3,k), ...
            sprintf('  M%d', k), 'Color', marsColor, 'FontSize', 8, 'FontWeight', 'bold');

        text(Jupiter_event_r(1,k), Jupiter_event_r(2,k), Jupiter_event_r(3,k), ...
            sprintf('  J%d', k), 'Color', jupiterColor, 'FontSize', 8, 'FontWeight', 'bold');
    end

    % -----------------------------
    % Label spacecraft event points
    % -----------------------------
    text(Earth_departure(1), Earth_departure(2), Earth_departure(3), ...
        '  Earth Departure', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

    text(Mars_flyby(1), Mars_flyby(2), Mars_flyby(3), ...
        '  Mars Flyby', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

    text(DSM_point(1), DSM_point(2), DSM_point(3), ...
        '  DSM', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

    text(Earth_flyby(1), Earth_flyby(2), Earth_flyby(3), ...
        '  Earth Flyby', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

    text(Jupiter_arrival(1), Jupiter_arrival(2), Jupiter_arrival(3), ...
        '  Jupiter Arrival', 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

    % -----------------------------
    % Optional event index guide
    % -----------------------------
    annotationText = sprintf([ ...
        'Planet markers at event times:\n' ...
        '1 = %s\n' ...
        '2 = %s\n' ...
        '3 = %s\n' ...
        '4 = %s\n' ...
        '5 = %s'], eventNames{:});

    text(ax.XLim(1), ax.YLim(2), ax.ZLim(2), annotationText, ...
        'Color', [0.95 0.95 0.95], 'FontSize', 8, ...
        'VerticalAlignment', 'top');

    xlabel('x (AU)');
    ylabel('y (AU)');
    zlabel('z (AU)');
    title('Full Interplanetary Mission: Earth \rightarrow Mars \rightarrow DSM \rightarrow Earth \rightarrow Jupiter', ...
        'Color', [0.95 0.95 0.95]);

    legend('TextColor', [0.95 0.95 0.95], ...
           'Color', [0.12 0.12 0.12], ...
           'Location', 'bestoutside');

    hold off;
end