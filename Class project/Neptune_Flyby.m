clc
clear


%% Interplanetary Mission Project (Readable Version)
% This script does two main things:
%   1) Propagates all planets from the J2000 epoch to 25-Dec-2026 and
%      reports heliocentric state vectors and true anomalies.
%   2) Builds one example Earth -> Mars flyby -> DSM -> Earth flyby ->
%      Jupiter mission using patched-conic style calculations.
%
% Notes on units used throughout this script:
%   - Distances in heliocentric propagation are typically in AU.
%   - Velocities in heliocentric propagation are typically in AU/TU.
%   - Planet-centered flyby/circular-orbit calculations use km and km/s.
%   - Angles stored in structures are in radians.

%% ------------------------------------------------------------------------
%  CONSTANTS AND PLANET DATA
%  ------------------------------------------------------------------------
muSunCanonical = 1;          % canonical solar gravitational parameter [AU^3/TU^2]
solarTU_days   = 58.13;      % 1 solar time unit in days
AU_to_km       = 1.495e8;    % astronomical unit to km
AU_TU_to_kms   = 29.79;      % canonical velocity conversion [km/s]

Sun = struct( ...
    'mass', 1.989e30, ...    % [kg]
    'mu',   1.327e11);       % [km^3/s^2]

Mercury = struct('a',0.387099,'e',0.205631,'inc',deg2rad(7.00487), ...
    'OMEGA',deg2rad(48.33167),'omega',deg2rad(29.12478), ...
    'theta',deg2rad(174.7944),'mu',22031.868551);

Venus = struct('a',0.723332,'e',0.006773,'inc',deg2rad(3.39471), ...
    'OMEGA',deg2rad(76.68069),'omega',deg2rad(54.85229), ...
    'theta',deg2rad(50.44675),'mu',324858.592);

Earth = struct('a',1.0,'e',0.01671,'inc',deg2rad(0.00005), ...
    'OMEGA',deg2rad(-11.26064),'omega',deg2rad(114.20783), ...
    'theta',deg2rad(-2.48284),'mu',398600.435507,'r',6378, ...
    'mass',5.972e24);

Mars = struct('a',1.523662,'e',0.093412,'inc',deg2rad(1.85061), ...
    'OMEGA',deg2rad(49.57854),'omega',deg2rad(286.4623), ...
    'theta',deg2rad(19.41248),'mu',4.305e4,'r',3380, ...
    'mass',6.41693e23);

Jupiter = struct('a',5.203363,'e',0.048393,'inc',deg2rad(1.3053), ...
    'OMEGA',deg2rad(100.55615),'omega',deg2rad(-85.8023), ...
    'theta',deg2rad(19.55053),'mu',1.26713e8,'r',71370, ...
    'mass',1.89852e27);

Saturn = struct('a',9.537070,'e',0.054151,'inc',deg2rad(2.48446), ...
    'OMEGA',deg2rad(113.71504),'omega',deg2rad(-21.2831), ...
    'theta',deg2rad(-42.4876),'mu',3.79406e7);

Uranus = struct('a',19.19126,'e',0.047168,'inc',deg2rad(0.76986), ...
    'OMEGA',deg2rad(74.22988),'omega',deg2rad(96.73436), ...
    'theta',deg2rad(142.2679),'mu',5.79456e6);

Neptune = struct('a',30.06896,'e',0.008586,'inc',deg2rad(1.76917), ...
    'OMEGA',deg2rad(131.72169),'omega',deg2rad(-86.75034), ...
    'theta',deg2rad(259.9087),'mu',6.83653e6,'r',22320);

Pluto = struct('a',39.48169,'e',0.248808,'inc',deg2rad(17.14175), ...
    'OMEGA',deg2rad(110.30347),'omega',deg2rad(113.76329), ...
    'theta',deg2rad(14.86205),'mu',975.500);

J2000epoch = juliandate(datetime(2000,1,1,11,58,0));

%% ------------------------------------------------------------------------
%  PLANET STATES AT J2000
%  ------------------------------------------------------------------------
[mercuryJ2000_r, mercuryJ2000_v] = posandvelvector(Mercury, muSunCanonical);
[venusJ2000_r,   venusJ2000_v]   = posandvelvector(Venus,   muSunCanonical);
[earthJ2000_r,   earthJ2000_v]   = posandvelvector(Earth,   muSunCanonical);
[marsJ2000_r,    marsJ2000_v]    = posandvelvector(Mars,    muSunCanonical);
[jupiterJ2000_r, jupiterJ2000_v] = posandvelvector(Jupiter, muSunCanonical);
[saturnJ2000_r,  saturnJ2000_v]  = posandvelvector(Saturn,  muSunCanonical);
[uranusJ2000_r,  uranusJ2000_v]  = posandvelvector(Uranus,  muSunCanonical);
[neptuneJ2000_r, neptuneJ2000_v] = posandvelvector(Neptune, muSunCanonical);
[plutoJ2000_r,   plutoJ2000_v]   = posandvelvector(Pluto,   muSunCanonical);

%% ------------------------------------------------------------------------
%  PART 1: PLANET TABLE ON 25-DEC-2026 23:12 UTC
%  ------------------------------------------------------------------------
reportDate = juliandate(datetime(2026,12,25,23,12,0));
timeSinceJ2000_TU = (reportDate - J2000epoch) / solarTU_days;

[Mercury_r, Mercury_v] = universalTOF(muSunCanonical, timeSinceJ2000_TU, mercuryJ2000_r, mercuryJ2000_v);
[Venus_r,   Venus_v]   = universalTOF(muSunCanonical, timeSinceJ2000_TU, venusJ2000_r,   venusJ2000_v);
[Earth_r,   Earth_v]   = universalTOF(muSunCanonical, timeSinceJ2000_TU, earthJ2000_r,   earthJ2000_v);
[Mars_r,    Mars_v]    = universalTOF(muSunCanonical, timeSinceJ2000_TU, marsJ2000_r,    marsJ2000_v);
[Jupiter_r, Jupiter_v] = universalTOF(muSunCanonical, timeSinceJ2000_TU, jupiterJ2000_r, jupiterJ2000_v);
[Saturn_r,  Saturn_v]  = universalTOF(muSunCanonical, timeSinceJ2000_TU, saturnJ2000_r,  saturnJ2000_v);
[Uranus_r,  Uranus_v]  = universalTOF(muSunCanonical, timeSinceJ2000_TU, uranusJ2000_r,  uranusJ2000_v);
[Neptune_r, Neptune_v] = universalTOF(muSunCanonical, timeSinceJ2000_TU, neptuneJ2000_r, neptuneJ2000_v);
[Pluto_r,   Pluto_v]   = universalTOF(muSunCanonical, timeSinceJ2000_TU, plutoJ2000_r,   plutoJ2000_v);

[~,~,~,~,~,~,thetaMercury] = orbitalelementscalc(Mercury_r, Mercury_v, muSunCanonical);
[~,~,~,~,~,~,thetaVenus]   = orbitalelementscalc(Venus_r,   Venus_v,   muSunCanonical);
[~,~,~,~,~,~,thetaEarth]   = orbitalelementscalc(Earth_r,   Earth_v,   muSunCanonical);
[~,~,~,~,~,~,thetaMars]    = orbitalelementscalc(Mars_r,    Mars_v,    muSunCanonical);
[~,~,~,~,~,~,thetaJupiter] = orbitalelementscalc(Jupiter_r, Jupiter_v, muSunCanonical);
[~,~,~,~,~,~,thetaSaturn]  = orbitalelementscalc(Saturn_r,  Saturn_v,  muSunCanonical);
[~,~,~,~,~,~,thetaUranus]  = orbitalelementscalc(Uranus_r,  Uranus_v,  muSunCanonical);
[~,~,~,~,~,~,thetaNeptune] = orbitalelementscalc(Neptune_r, Neptune_v, muSunCanonical);
[~,~,~,~,~,~,thetaPluto]   = orbitalelementscalc(Pluto_r,   Pluto_v,   muSunCanonical);

planetNames = ["Mercury"; "Venus"; "Earth"; "Mars"; "Jupiter"; ...
               "Saturn"; "Uranus"; "Neptune"; "Pluto"];

R = [Mercury_r(:)'; Venus_r(:)'; Earth_r(:)'; Mars_r(:)'; Jupiter_r(:)'; ...
     Saturn_r(:)'; Uranus_r(:)'; Neptune_r(:)'; Pluto_r(:)'];

V = [Mercury_v(:)'; Venus_v(:)'; Earth_v(:)'; Mars_v(:)'; Jupiter_v(:)'; ...
     Saturn_v(:)'; Uranus_v(:)'; Neptune_v(:)'; Pluto_v(:)'];

trueAnomaly_deg = rad2deg([thetaMercury; thetaVenus; thetaEarth; thetaMars; ...
                           thetaJupiter; thetaSaturn; thetaUranus; ...
                           thetaNeptune; thetaPluto]);

planetStateTable = table(planetNames, ...
    R(:,1), R(:,2), R(:,3), V(:,1), V(:,2), V(:,3), trueAnomaly_deg, ...
    'VariableNames', {'Planet','rx_AU','ry_AU','rz_AU', ...
                      'vx_AU_per_TU','vy_AU_per_TU','vz_AU_per_TU', ...
                      'TrueAnomaly_deg'});

disp(planetStateTable)
plotStart_TU = 0;
plotResolution = 500;
Planet_TOF_final=1500;
figure(7)
planet_plot(Mercury_r, Mercury_v, Venus_r, Venus_v, Earth_r, Earth_v, ...
    Mars_r, Mars_v, Jupiter_r, Jupiter_v, Saturn_r, Saturn_v, Uranus_r, Uranus_v, ...
    Neptune_r, Neptune_v, Pluto_r, Pluto_v, plotStart_TU, Planet_TOF_final, plotResolution)

%% ------------------------------------------------------------------------
%  PART 2: EARTH -> MARS -> DSM -> EARTH -> JUPITER MISSION
%  ------------------------------------------------------------------------
selectedDepartureDate = datetime(2026,9,28, 23,12,00);
departureJulianDate   = juliandate(selectedDepartureDate);
departureFromJ2000_TU = (departureJulianDate - J2000epoch) / solarTU_days;
marsarivaltime=juliandate(datetime(2027,8,14,23,12,00));
transfer1_TOF_TU=(marsarivaltime-departureJulianDate)/58.13;

% User-selected mission design values.
marsToDSM_TOF_TU     = 50  / solarTU_days;  % Mars flyby -> DSM
DSMtoEarth_TOF_TU    = 400/ solarTU_days;  % DSM -> Earth flyby
EarthToJupiter_TOF_TU = 740 / solarTU_days; % Earth flyby -> Jupiter
JupiterToDSM_TOF_TU=1000 / solarTU_days; % Jupiter flyby -> DSM
DSMToNeptune_TOF_TU=800 / solarTU_days; % DSM -> Neptune


earthParkingOrbitAlt_km   = 500;
marsFlybyAltitude_km      = 100;
earthFlybyAltitude_km     = 100;
jupiterParkingOrbitAlt_km = 10000;
jupiterFlybyAltitude_km= 10000;
useShortWayTransfer       = 0;
useShortWayTransferDSM=1;
neptuneParkingOrbitAlt_km= 21000

% Jupiter state at the chosen Earth departure date.
[Jupiter_departure_r, Jupiter_departure_v] = universalTOF( ...
    muSunCanonical, departureFromJ2000_TU, jupiterJ2000_r, jupiterJ2000_v);
% Neptune state at the chosen Earth departure date.
[Neptune_departure_r, Neptune_departure_v] = universalTOF( ...
    muSunCanonical, departureFromJ2000_TU, neptuneJ2000_r, neptuneJ2000_v);

%% Leg 1: Earth departure to Mars flyby
[earthDepart_r, earthDepart_v, marsDepart_r, marsDepart_v, ...
    marsArrival_r, marsArrival_v, transfer1Struct, departureDeltaV_kms, ...
    marsTurnAngle_rad, marsVinfOut_kms, marsVinfIn_kms, ...
    scAfterMarsFlyby_v, transfer1_v1, transfer1_v2] = ...
    eathdepatruetomars(AU_to_km, Sun, useShortWayTransfer, ...
    earthJ2000_r, earthJ2000_v, marsJ2000_r, marsJ2000_v, ...
    departureJulianDate, transfer1_TOF_TU, Earth, Mars, ...
    earthParkingOrbitAlt_km, muSunCanonical, AU_TU_to_kms, ...
    J2000epoch, marsFlybyAltitude_km);

fprintf('\n--- Leg 1: Earth to Mars ---\n');
fprintf('Earth departure position [AU]:\n'); disp(earthDepart_r);
fprintf('Earth departure velocity [AU/TU]:\n'); disp(earthDepart_v);
fprintf('Mars arrival position [AU]:\n'); disp(marsArrival_r);
fprintf('Mars arrival velocity [AU/TU]:\n'); disp(marsArrival_v);
fprintf('Transfer departure velocity [AU/TU]:\n'); disp(transfer1_v1);
fprintf('Transfer arrival velocity [AU/TU]:\n'); disp(transfer1_v2);
fprintf('Earth departure delta-V [km/s]: %g\n', departureDeltaV_kms);
fprintf('Mars flyby turn angle [deg]: %g\n', rad2deg(marsTurnAngle_rad));
fprintf('Mars incoming V-infinity [km/s]:\n'); disp(marsVinfIn_kms);
fprintf('Mars outgoing V-infinity [km/s]:\n'); disp(marsVinfOut_kms);
fprintf('Post-Mars-flyby heliocentric spacecraft velocity [km/s]:\n');
disp(scAfterMarsFlyby_v * AU_TU_to_kms);

%% Plot launch-to-Mars leg
plotStart_TU = 0;
plotResolution = 500;

[earthAfterLeg1_r, earthAfterLeg1_v, marsAfterLeg1_r, marsAfterLeg1_v, ...
    jupiterAfterLeg1_r, jupiterAfterLeg1_v,neptuneAfterLeg1_r, neptuneAfterLeg1_v, scAfterLeg1_r, scAfterLeg1_v] = ...
    planet_plotearthsmall(earthDepart_r, earthDepart_v, marsDepart_r, ...
    marsDepart_v, Jupiter_departure_r, Jupiter_departure_v, Neptune_departure_r, Neptune_departure_v,...
    plotStart_TU, transfer1_TOF_TU, plotResolution, earthDepart_r, transfer1_v1);

%% Leg 2: Mars flyby to DSM
[transfer2_a, transfer2_eVec, transfer2_e, transfer2_i, ...
    transfer2_RAAN, transfer2_argPeri, transfer2_trueAnom] = ...
    orbitalelementscalc(marsArrival_r, scAfterMarsFlyby_v, muSunCanonical);

fprintf('\n--- Leg 2 Orbit After Mars Flyby ---\n');
fprintf('Semi-major axis [AU]: %g\n', transfer2_a);
fprintf('Eccentricity vector:\n'); disp(transfer2_eVec);
fprintf('Eccentricity magnitude: %g\n', transfer2_e);
fprintf('Inclination [deg]: %g\n', rad2deg(transfer2_i));
fprintf('RAAN [deg]: %g\n', rad2deg(transfer2_RAAN));
fprintf('Argument of periapsis [deg]: %g\n', rad2deg(transfer2_argPeri));
fprintf('True anomaly [deg]: %g\n', rad2deg(transfer2_trueAnom));

[earthAtDSM_r, earthAtDSM_v] = universalTOF(muSunCanonical, marsToDSM_TOF_TU, earthAfterLeg1_r, earthAfterLeg1_v);
[spacecraftAtDSM_r, spacecraftAtDSM_v] = universalTOF(muSunCanonical, marsToDSM_TOF_TU, scAfterLeg1_r, scAfterMarsFlyby_v);
[JupiterAtDSM_r, JupiterAtDSM_v] = universalTOF(muSunCanonical, marsToDSM_TOF_TU, jupiterAfterLeg1_r, jupiterAfterLeg1_v);
[MarsAtDSM_r, MarsAtDSM_v] = universalTOF(muSunCanonical, marsToDSM_TOF_TU, marsAfterLeg1_r, marsAfterLeg1_v);
[NeptuneAtDSM_r, NeptuneAtDSM_v] = universalTOF(muSunCanonical, marsToDSM_TOF_TU, neptuneAfterLeg1_r, neptuneAfterLeg1_v);


%% Leg 3: DSM to Earth flyby
[EarthFlyby_r, EarthFlyby_v] = universalTOF(muSunCanonical, DSMtoEarth_TOF_TU, earthAtDSM_r, earthAtDSM_v);

[transfer3_v1, transfer3_v2, transfer3Struct, transfer3_r1, transfer3_r2] = ...
    gaussspeedsandtransferorbitorbitalelements(spacecraftAtDSM_r, EarthFlyby_r, ...
    DSMtoEarth_TOF_TU, muSunCanonical, useShortWayTransferDSM);

fprintf('\n--- Leg 3: DSM to Earth Flyby ---\n');
fprintf('Velocity at DSM [AU/TU]: [%g, %g, %g]\n', transfer3_v1);
fprintf('Velocity at Earth flyby [AU/TU]: [%g, %g, %g]\n', transfer3_v2);
fprintf('Semi-major axis [AU]: %g\n', transfer3Struct.a);
fprintf('Eccentricity: %g\n', transfer3Struct.e);
fprintf('Inclination [deg]: %g\n', rad2deg(transfer3Struct.inc));
fprintf('RAAN [deg]: %g\n', rad2deg(transfer3Struct.OMEGA));
fprintf('Argument of periapsis [deg]: %g\n', rad2deg(transfer3Struct.omega));
fprintf('True anomaly at r1 [deg]: %g\n', rad2deg(transfer3Struct.theta));

[earthAfterDSM_r, earthAfterDSM_v, marsAfterDSM_r, marsAfterDSM_v, ...
    jupiterAfterDSM_r, jupiterAfterDSM_v,neptuneAfterDSM_r, neptuneAfterDSM_v, scAfterDSM_r, scAfterDSM_v] = ...
    planet_plotNeptunesmall(earthAfterLeg1_r, earthAfterLeg1_v, marsAfterLeg1_r, ...
    marsAfterLeg1_v, jupiterAfterLeg1_r, jupiterAfterLeg1_v, neptuneAfterLeg1_r, neptuneAfterLeg1_v,...
    plotStart_TU, marsToDSM_TOF_TU, plotResolution, scAfterLeg1_r, ...
    scAfterMarsFlyby_v, 2);

DSM_deltaV_kms = norm(transfer3_v1 - spacecraftAtDSM_v) * AU_TU_to_kms;
fprintf('\nDSM delta-V [km/s]: %g\n', DSM_deltaV_kms);

[earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, marsAfterFlybyPlot_r, ...
    marsAfterFlybyPlot_v, jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, neptuneAfterFlybyPlot_r, neptuneAfterFlybyPlot_v,...
    scBeforeEarthFlyby_r, scBeforeEarthFlyby_v] = ...
    planet_plotNeptunesmall(earthAfterDSM_r, earthAfterDSM_v, marsAfterDSM_r, ...
    marsAfterDSM_v, jupiterAfterDSM_r, jupiterAfterDSM_v, neptuneAfterDSM_r, neptuneAfterDSM_v,plotStart_TU, ...
    DSMtoEarth_TOF_TU, plotResolution, scAfterDSM_r, transfer3_v1, 3);

%% Earth flyby and final transfer to Jupiter
[earthTurnAngle_rad, earthVinfOut_kms, earthVinfIn_kms, scAfterEarthFlyby_v] = ...
    hyperbolicturnangle(transfer3_v2, EarthFlyby_v, AU_TU_to_kms, ...
    earthFlybyAltitude_km, Earth, Sun, AU_to_km);

[transfer4_a, transfer4_eVec, transfer4_e, transfer4_i, ...
    transfer4_RAAN, transfer4_argPeri, transfer4_trueAnom] = ...
    orbitalelementscalc(transfer3_r2, scAfterEarthFlyby_v, muSunCanonical);

fprintf('\n--- Leg 4 Orbit After Earth Flyby ---\n');
fprintf('Semi-major axis [AU]: %g\n', transfer4_a);
fprintf('Eccentricity vector:\n'); disp(transfer4_eVec);
fprintf('Eccentricity magnitude: %g\n', transfer4_e);
fprintf('Inclination [deg]: %g\n', rad2deg(transfer4_i));
fprintf('RAAN [deg]: %g\n', rad2deg(transfer4_RAAN));
fprintf('Argument of periapsis [deg]: %g\n', rad2deg(transfer4_argPeri));
fprintf('True anomaly [deg]: %g\n', rad2deg(transfer4_trueAnom));

%% Earth->Jupiter 
[Earth_final_r, Earth_final_v, Mars_final_r, Mars_final_v, ...
    Jupiter_final_r, Jupiter_final_v,Neptune_final_r, Neptune_final_v, SC_final_r, SC_final_v] = ...
    planet_plotNeptunesmall(earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, ...
    marsAfterFlybyPlot_r, marsAfterFlybyPlot_v, ...
    jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, neptuneAfterFlybyPlot_r, neptuneAfterFlybyPlot_v, ...
    plotStart_TU, EarthToJupiter_TOF_TU, plotResolution, ...
    scBeforeEarthFlyby_r, scAfterEarthFlyby_v, 4);
jupiterVinf_kms = (SC_final_v - Jupiter_final_v) * AU_TU_to_kms;

jupiterArrivalDeltaV_kms = deltavforcirculartohyperbolic( ...
    jupiterParkingOrbitAlt_km, Jupiter, jupiterVinf_kms);

%% Jupiter Flyby
[jupiterTurnAngle_rad, jupiterVinfOut_kms, jupiterVinfIn_kms, scAfterjupiterFlyby_v] = ...
    hyperbolicturnangle(SC_final_v, Jupiter_final_v, AU_TU_to_kms, ...
    jupiterFlybyAltitude_km, Jupiter, Sun, AU_to_km);


[transfer5_a, transfer5_eVec, transfer5_e, transfer5_i, ...
    transfer5_RAAN, transfer5_argPeri, transfer5_trueAnom] = ...
    orbitalelementscalc(SC_final_r, scAfterjupiterFlyby_v, muSunCanonical);

[Earth_AfterJ_r, Earth__AfterJ_v, ...
          Mars__AfterJ_r, Mars__AfterJ_v, ...
          Jupiter__AfterJ_r, Jupiter__AfterJ_v,Neptune__AfterJ_r, Neptune__AfterJ_v, ...
          SC__AfterJ_r, SC__AfterJ_v] = ...
          planet_plotNeptunesmall(Earth_final_r, Earth_final_v, Mars_final_r, Mars_final_v, ...
         Jupiter_final_r, Jupiter_final_v, Neptune_final_r, Neptune_final_v,  plotStart_TU, JupiterToDSM_TOF_TU, ...
          plotResolution, SC_final_r, SC_final_v, 8)

[earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, marsAfterFlybyPlot_r, ...
    marsAfterFlybyPlot_v, jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, neptuneAfterFlybyPlot_r, neptuneAfterFlybyPlot_v, ...
    scBeforeEarthFlyby_r, scBeforeEarthFlyby_v] = ...
    planet_plotNeptunesmall(earthAfterDSM_r, earthAfterDSM_v, marsAfterDSM_r, ...
    marsAfterDSM_v, jupiterAfterDSM_r, jupiterAfterDSM_v, neptuneAfterDSM_r, neptuneAfterDSM_v,plotStart_TU, ...
    DSMtoEarth_TOF_TU, plotResolution, scAfterDSM_r, transfer3_v1, 6);


[earthAtDSM2_r, earthAtDSM2_v] = universalTOF(muSunCanonical, JupiterToDSM_TOF_TU, earthAfterLeg1_r, earthAfterLeg1_v);
[spacecraftAtDSM2_r, spacecraftAtDSM2_v] = universalTOF(muSunCanonical,  JupiterToDSM_TOF_TU, scAfterLeg1_r, scAfterMarsFlyby_v);
[JupiterAtDSM2_r, JupiterAtDSM2_v] = universalTOF(muSunCanonical,  JupiterToDSM_TOF_TU, jupiterAfterLeg1_r, jupiterAfterLeg1_v);
[MarsAtDSM2_r, MarsAtDSM2_v] = universalTOF(muSunCanonical,  JupiterToDSM_TOF_TU, marsAfterLeg1_r, marsAfterLeg1_v);
[NeptuneAtDSM2_r, NeptuneAtDSM2_v] = universalTOF(muSunCanonical,  JupiterToDSM_TOF_TU, neptuneAfterLeg1_r, neptuneAfterLeg1_v);

totalMissionDeltaV_kms = departureDeltaV_kms + DSM_deltaV_kms + jupiterArrivalDeltaV_kms;

fprintf('\n--- Jupiter Arrival ---\n');
fprintf('Jupiter hyperbolic excess velocity [km/s]:\n'); disp(jupiterVinf_kms);
fprintf('Jupiter arrival delta-V [km/s]: %g\n', jupiterArrivalDeltaV_kms);

fprintf('\n--- Mission Summary ---\n');
fprintf('Leg 1 TOF Earth->Mars [TU]: %g\n', transfer1_TOF_TU);
fprintf('Leg 2 TOF Mars->DSM [TU]: %g\n', marsToDSM_TOF_TU);
fprintf('Leg 3 TOF DSM->Earth [TU]: %g\n', DSMtoEarth_TOF_TU);
fprintf('Leg 4 TOF Earth->Jupiter [TU]: %g\n', EarthToJupiter_TOF_TU);
fprintf('Earth departure delta-V [km/s]: %g\n', departureDeltaV_kms);
fprintf('DSM delta-V [km/s]: %g\n', DSM_deltaV_kms);
fprintf('Jupiter arrival delta-V [km/s]: %g\n', jupiterArrivalDeltaV_kms);
fprintf('Total mission delta-V [km/s]: %g\n', totalMissionDeltaV_kms);


%% ========================================================================
% FULL MISSION 3D PLOT
% ========================================================================
plot_full_mission_3d( ...
    earthDepart_r, earthDepart_v, ...                  % Earth at launch
    marsDepart_r, marsDepart_v, ...                    % Mars at launch epoch
    Jupiter_departure_r, Jupiter_departure_v,... 
    Neptune_departure_r, Neptune_departure_v,...      % Jupiter at launch epoch
    earthDepart_r, transfer1_v1, transfer1_TOF_TU, ...            % Leg 1: Earth -> Mars
    marsArrival_r, scAfterMarsFlyby_v, marsToDSM_TOF_TU, ...      % Leg 2: Mars -> DSM
    spacecraftAtDSM_r, transfer3_v1, DSMtoEarth_TOF_TU, ...       % Leg 3: DSM -> Earth
    transfer3_r2, scAfterEarthFlyby_v, EarthToJupiter_TOF_TU, ... % Leg 4: Earth -> Jupiter
    1000,5);
%% ========================================================================
% RESULTS TABLES FOR REPORT
% ========================================================================

%% -------------------------------
% Table 1: Key Mission Event States
% -------------------------------
eventNames = ["Earth Departure";
              "Mars Flyby";
              "DSM";
              "Earth Flyby";
              "Jupiter Arrival"];

eventR = [earthDepart_r(:)';
          marsArrival_r(:)';
          spacecraftAtDSM_r(:)';
          transfer3_r2(:)';
          Jupiter_final_r(:)'];

eventV = [earthDepart_v(:)';
          marsArrival_v(:)';
          spacecraftAtDSM_v(:)';
          EarthFlyby_v(:)';
          Jupiter_final_v(:)'];

missionEventTable = table( ...
    eventNames, ...
    eventR(:,1), eventR(:,2), eventR(:,3), ...
    eventV(:,1), eventV(:,2), eventV(:,3), ...
    'VariableNames', {'Event','rx_AU','ry_AU','rz_AU', ...
                      'vx_AU_per_TU','vy_AU_per_TU','vz_AU_per_TU'});

disp('================ Mission Event Table ================')
disp(missionEventTable)

%% -------------------------------
% Table 2: Transfer Orbit Elements
% -------------------------------
transferNames = ["Leg 1 Earth->Mars";
                 "Leg 2 Mars->DSM";
                 "Leg 3 DSM->Earth";
                 "Leg 4 Earth->Jupiter"];

orbitElementTable = table( ...
    transferNames, ...
    [transfer1Struct.a; transfer2_a; transfer3Struct.a; transfer4_a], ...
    [transfer1Struct.e; transfer2_e; transfer3Struct.e; transfer4_e], ...
    rad2deg([transfer1Struct.inc; transfer2_i; transfer3Struct.inc; transfer4_i]), ...
    rad2deg([transfer1Struct.OMEGA; transfer2_RAAN; transfer3Struct.OMEGA; transfer4_RAAN]), ...
    rad2deg([transfer1Struct.omega; transfer2_argPeri; transfer3Struct.omega; transfer4_argPeri]), ...
    'VariableNames', {'TransferLeg','a_AU','e','inc_deg','RAAN_deg','argPeri_deg'});

disp('================ Transfer Orbit Elements ================')
disp(orbitElementTable)

%% -------------------------------
% Table 3: TOF and Delta-V Summary
% -------------------------------
tofDays = [transfer1_TOF_TU;
           marsToDSM_TOF_TU;
           DSMtoEarth_TOF_TU;
           EarthToJupiter_TOF_TU] * solarTU_days;

deltaVSummaryTable = table( ...
    ["Leg 1 Earth->Mars";
     "Leg 2 Mars->DSM";
     "Leg 3 DSM->Earth";
     "Leg 4 Earth->Jupiter";
     "TOTAL"], ...
    [tofDays; sum(tofDays)], ...
    [departureDeltaV_kms;
     0;
     DSM_deltaV_kms;
     jupiterArrivalDeltaV_kms;
     totalMissionDeltaV_kms], ...
    'VariableNames', {'MissionSegment','TOF_days','DeltaV_kms'});

disp('================ TOF and Delta-V Summary ================')
disp(deltaVSummaryTable)
%% ========================================================================
%  LOCAL PLOTTING FUNCTIONS
%  ========================================================================
function [Earth_final_r, Earth_final_v, ...
          Mars_final_r, Mars_final_v, ...
          Jupiter_final_r, Jupiter_final_v, Neptune_final_r, Neptune_final_v,...
          SC_final_r, SC_final_v] = ...
          planet_plotearthsmall(Earth_r, Earth_v, Mars_r, Mars_v, ...
          Jupiter_r, Jupiter_v,Neptune_r, Neptune_v, Planet_TOF_initial, Planet_TOF_final, ...
          Planet_TOF_res, Spacecraft_r, Spacecraft_v)
% Plot Earth, Mars, Jupiter, and an optional spacecraft over one time span.

    figure(1)
    clf
    set(gcf, 'Color', [0.06 0.06 0.06]);

    muSunCanonical = 1;
    TOF = linspace(Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res);

    Earth_orbit_r   = zeros(3, length(TOF));
    Mars_orbit_r    = zeros(3, length(TOF));
    Jupiter_orbit_r = zeros(3, length(TOF));
    Neptune_orbit_r = zeros(3, length(TOF));
    SC_orbit_r      = zeros(3, length(TOF));

    Earth_orbit_v   = zeros(3, length(TOF));
    Mars_orbit_v    = zeros(3, length(TOF));
    Jupiter_orbit_v = zeros(3, length(TOF));
    Neptune_orbit_v = zeros(3, length(TOF));
    SC_orbit_v      = zeros(3, length(TOF));

    for k = 1:length(TOF)
        [Earth_orbit_r(:,k), Earth_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Earth_r, Earth_v);
        [Mars_orbit_r(:,k),  Mars_orbit_v(:,k)]  = universalTOF(muSunCanonical, TOF(k), Mars_r, Mars_v);
        [Jupiter_orbit_r(:,k), Jupiter_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Jupiter_r, Jupiter_v);
        [Neptune_orbit_r(:,k), Neptune_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Neptune_r, Neptune_v);

        if nargin >= 11 && ~isempty(Spacecraft_r) && ~isempty(Spacecraft_v)
            [SC_orbit_r(:,k), SC_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Spacecraft_r, Spacecraft_v);
        end
    end

    makeFourBodyPlot(Earth_r, Mars_r, Jupiter_r, Neptune_r,Earth_orbit_r, Mars_orbit_r, ...
        Jupiter_orbit_r,Neptune_orbit_r, SC_orbit_r, Spacecraft_r, 1);

    Earth_final_r = Earth_orbit_r(:,end); Earth_final_v = Earth_orbit_v(:,end);
    Mars_final_r = Mars_orbit_r(:,end);   Mars_final_v = Mars_orbit_v(:,end);
    Jupiter_final_r = Jupiter_orbit_r(:,end); Jupiter_final_v = Jupiter_orbit_v(:,end);
    Neptune_final_r = Neptune_orbit_r(:,end); Neptune_final_v = Neptune_orbit_v(:,end);

    if any(SC_orbit_r(:))
        SC_final_r = SC_orbit_r(:,end);
        SC_final_v = SC_orbit_v(:,end);
    else
        SC_final_r = [];
        SC_final_v = [];
    end
end


function [Earth_final_r, Earth_final_v, ...
          Mars_final_r, Mars_final_v, ...
          Jupiter_final_r, Jupiter_final_v,Neptune_final_r, Neptune_final_v, ...
          SC_final_r, SC_final_v] = ...
          planet_plotNeptunesmall(Earth_r, Earth_v, Mars_r, Mars_v, ...
          Jupiter_r, Jupiter_v, Neptune_r, Neptune_v, Planet_TOF_initial, Planet_TOF_final, ...
          Planet_TOF_res, Spacecraft_r, Spacecraft_v, figNumber)
% Same plotting helper as above, but allows caller to choose figure number.

    figure(figNumber)
    clf
    set(gcf, 'Color', [0.06 0.06 0.06]);

    muSunCanonical = 1;
    TOF = linspace(Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res);

    Earth_orbit_r   = zeros(3, length(TOF));
    Mars_orbit_r    = zeros(3, length(TOF));
    Jupiter_orbit_r = zeros(3, length(TOF));
    Neptune_orbit_r = zeros(3, length(TOF));
    SC_orbit_r      = zeros(3, length(TOF));

    Earth_orbit_v   = zeros(3, length(TOF));
    Mars_orbit_v    = zeros(3, length(TOF));
    Neptune_orbit_v = zeros(3, length(TOF));
    Jupiter_orbit_v = zeros(3, length(TOF));
    SC_orbit_v      = zeros(3, length(TOF));

    for k = 1:length(TOF)
        [Earth_orbit_r(:,k), Earth_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Earth_r, Earth_v);
        [Mars_orbit_r(:,k),  Mars_orbit_v(:,k)]  = universalTOF(muSunCanonical, TOF(k), Mars_r, Mars_v);
        [Jupiter_orbit_r(:,k), Jupiter_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Jupiter_r, Jupiter_v);
        [Neptune_orbit_r(:,k), Neptune_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Neptune_r, Neptune_v);

        if nargin >= 11 && ~isempty(Spacecraft_r) && ~isempty(Spacecraft_v)
            [SC_orbit_r(:,k), SC_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Spacecraft_r, Spacecraft_v);
        end
    end

    makeFourBodyPlot(Earth_r, Mars_r, Jupiter_r, Neptune_r, Earth_orbit_r, Mars_orbit_r, ...
        Jupiter_orbit_r,Neptune_orbit_r, SC_orbit_r, Spacecraft_r, figNumber);

    Earth_final_r = Earth_orbit_r(:,end); Earth_final_v = Earth_orbit_v(:,end);
    Mars_final_r = Mars_orbit_r(:,end);   Mars_final_v = Mars_orbit_v(:,end);
    Jupiter_final_r = Jupiter_orbit_r(:,end); Jupiter_final_v = Jupiter_orbit_v(:,end);
    Neptune_final_r = Neptune_orbit_r(:,end); Neptune_final_v = Neptune_orbit_v(:,end);

    if any(SC_orbit_r(:))
        SC_final_r = SC_orbit_r(:,end);
        SC_final_v = SC_orbit_v(:,end);
    else
        SC_final_r = [];
        SC_final_v = [];
    end
end

function makeFourBodyPlot(Earth_r, Mars_r, Jupiter_r, Neptune_r, Earth_orbit_r, Mars_orbit_r, ...
    Jupiter_orbit_r, Neptune_orbit_r, SC_orbit_r, Spacecraft_r, figNumber)
% Shared plotting utility used by both local plotting functions.

    figure(figNumber)
    ax = axes('Parent', gcf, ...
        'Color', [0.03 0.03 0.03], ...
        'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
        'ZColor', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
    hold(ax, 'on')

    sunColor = [1.0, 0.85, 0.0];
    earthColor = [0.2, 0.6, 1.0];
    marsColor = [1.0, 0.5, 0.2];
    jupiterColor = [1.0, 0.85, 0.7];
    neptuneColor= [0.2, 0.6, 1.0];
    orbitColor = [0.7, 0.7, 0.7];
    scColor = [1.0, 0.2, 0.2];
    finalMarkerColor = [0.9, 0.9, 0.9];

    plot3(ax, 0, 0, 0, 'o', 'MarkerFaceColor', sunColor, 'MarkerEdgeColor', sunColor, 'MarkerSize', 8, 'DisplayName', 'Sun')

    plot3(ax, Earth_r(1), Earth_r(2), Earth_r(3), 'o', 'MarkerFaceColor', earthColor, 'MarkerEdgeColor', earthColor, 'DisplayName', 'Earth')
    plot3(ax, Mars_r(1), Mars_r(2), Mars_r(3), 'o', 'MarkerFaceColor', marsColor, 'MarkerEdgeColor', marsColor, 'DisplayName', 'Mars')
    plot3(ax, Jupiter_r(1), Jupiter_r(2), Jupiter_r(3), 'o', 'MarkerFaceColor', jupiterColor, 'MarkerEdgeColor', jupiterColor, 'DisplayName', 'Jupiter')
    plot3(ax, Neptune_r(1), Neptune_r(2), Neptune_r(3), 'o', 'MarkerFaceColor', neptuneColor, 'MarkerEdgeColor', neptuneColor, 'DisplayName', 'Neptune')

    plot3(ax, Earth_orbit_r(1,:), Earth_orbit_r(2,:), Earth_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')
    plot3(ax, Mars_orbit_r(1,:), Mars_orbit_r(2,:), Mars_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')
    plot3(ax, Jupiter_orbit_r(1,:), Jupiter_orbit_r(2,:), Jupiter_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')
    plot3(ax, Neptune_orbit_r(1,:), Neptune_orbit_r(2,:), Neptune_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')

    if any(SC_orbit_r(:))
        plot3(ax, SC_orbit_r(1,:), SC_orbit_r(2,:), SC_orbit_r(3,:), '-', 'Color', scColor, 'LineWidth', 1.5, 'DisplayName', 'Spacecraft')
        plot3(ax, Spacecraft_r(1), Spacecraft_r(2), Spacecraft_r(3), 's', 'MarkerFaceColor', scColor, 'MarkerEdgeColor', scColor, 'DisplayName', 'SC Start')
        SC_final_r = SC_orbit_r(:,end);
        plot3(ax, SC_final_r(1), SC_final_r(2), SC_final_r(3), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'MarkerSize', 8, 'LineWidth', 1.2, 'DisplayName', 'SC Final')
    end

    plot3(ax, Earth_orbit_r(1,end), Earth_orbit_r(2,end), Earth_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')
    plot3(ax, Mars_orbit_r(1,end), Mars_orbit_r(2,end), Mars_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')
    plot3(ax, Jupiter_orbit_r(1,end), Jupiter_orbit_r(2,end), Jupiter_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')
    plot3(ax, Neptune_orbit_r(1,end), Neptune_orbit_r(2,end), Jupiter_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')

    legend(ax, 'TextColor', [0.95 0.95 0.95], 'Color', [0.12 0.12 0.12])
    title(ax, 'Solar System (Earth, Mars, Jupiter)', 'Color', [0.95 0.95 0.95])
    xlabel(ax, 'x (AU)', 'Color', [0.95 0.95 0.95])
    ylabel(ax, 'y (AU)', 'Color', [0.95 0.95 0.95])
    zlabel(ax, 'z (AU)', 'Color', [0.95 0.95 0.95])
    axis(ax, 'vis3d')
    daspect(ax, [1 1 1])
    grid(ax, 'on')
    view(ax, 3)
    hold(ax, 'off')
end
function [bestTOF_days, bestTOF_TU, bestSeparation_AU, ...
          bestEarth_r, bestEarth_v, bestMars_r, bestMars_v, ...
          bestJupiter_r, bestJupiter_v, bestSC_r, bestSC_v] = ...
    iterate_jupiter_tof_match( ...
    earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, ...
    marsAfterFlybyPlot_r, marsAfterFlybyPlot_v, ...
    jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, ...
    plotStart_TU, plotResolution, ...
    scBeforeEarthFlyby_r, scAfterEarthFlyby_v, ...
    TOF_days_start, TOF_days_end, TOF_days_step, ...
    solarTU_days, makePlots, figNumber, JupiterStruct, SunStruct)
% ITERATE_JUPITER_TOF_MATCH
% Sweep Earth->Jupiter time of flight in days and stop when the spacecraft
% ends inside Jupiter's sphere of influence.
%
% Tolerance is computed internally as Jupiter's SOI in AU:
%   r_SOI = a * (m_planet / m_sun)^(2/5)

    if nargin < 17 || isempty(makePlots)
        makePlots = false;
    end

    if nargin < 18 || isempty(figNumber)
        figNumber = 4;
    end

    % Jupiter sphere of influence in AU
    % matchTolerance_AU = JupiterStruct.a * (JupiterStruct.mass / SunStruct.mass)^(2/5);
    matchTolerance_AU=0.05; %AU

    fprintf('\nUsing Jupiter SOI as tolerance:\n');
    fprintf('Jupiter SOI = %.6f AU\n', matchTolerance_AU);

    bestSeparation_AU = inf;
    bestTOF_days = NaN;
    bestTOF_TU = NaN;

    bestEarth_r = [];
    bestEarth_v = [];
    bestMars_r = [];
    bestMars_v = [];
    bestJupiter_r = [];
    bestJupiter_v = [];
    bestSC_r = [];
    bestSC_v = [];

    for tof_days = TOF_days_start:TOF_days_step:TOF_days_end
        tof_TU = tof_days / solarTU_days;

        [Earth_final_r, Earth_final_v, Mars_final_r, Mars_final_v, ...
            Jupiter_final_r, Jupiter_final_v, SC_final_r, SC_final_v] = ...
            planet_plotmarssmall(earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, ...
            marsAfterFlybyPlot_r, marsAfterFlybyPlot_v, ...
            jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, ...
            plotStart_TU, tof_TU, plotResolution, ...
            scBeforeEarthFlyby_r, scAfterEarthFlyby_v, figNumber);

        separation_AU = norm(SC_final_r - Jupiter_final_r);

        fprintf('TOF = %8.2f days | Separation = %.6f AU\n', ...
            tof_days, separation_AU);

        if separation_AU < bestSeparation_AU
            bestSeparation_AU = separation_AU;
            bestTOF_days = tof_days;
            bestTOF_TU = tof_TU;

            bestEarth_r = Earth_final_r;
            bestEarth_v = Earth_final_v;
            bestMars_r = Mars_final_r;
            bestMars_v = Mars_final_v;
            bestJupiter_r = Jupiter_final_r;
            bestJupiter_v = Jupiter_final_v;
            bestSC_r = SC_final_r;
            bestSC_v = SC_final_v;
        end

        if separation_AU <= matchTolerance_AU
            fprintf('\nSpacecraft entered Jupiter SOI.\n');
            fprintf('Best TOF = %.2f days (%.6f TU)\n', bestTOF_days, bestTOF_TU);
            fprintf('Final separation = %.6f AU\n', bestSeparation_AU);
            return;
        end

        if ~makePlots
            clf(figNumber);
        end
    end

    fprintf('\nNo TOF entered Jupiter SOI in the tested range.\n');
    fprintf('Best TOF tested = %.2f days (%.6f TU)\n', bestTOF_days, bestTOF_TU);
    fprintf('Minimum separation = %.6f AU\n', bestSeparation_AU);

    % Replot best case found
    planet_plotmarssmall(earthAfterFlybyPlot_r, earthAfterFlybyPlot_v, ...
        marsAfterFlybyPlot_r, marsAfterFlybyPlot_v, ...
        jupiterAfterFlybyPlot_r, jupiterAfterFlybyPlot_v, ...
        plotStart_TU, bestTOF_TU, plotResolution, ...
        scBeforeEarthFlyby_r, scAfterEarthFlyby_v, figNumber);
end


function [Earth_final_r, Earth_final_v, ...
          Mars_final_r, Mars_final_v, ...
          Jupiter_final_r, Jupiter_final_v, ...
          SC_final_r, SC_final_v] = ...
          planet_plotmarssmall(Earth_r, Earth_v, Mars_r, Mars_v, ...
          Jupiter_r, Jupiter_v, Planet_TOF_initial, Planet_TOF_final, ...
          Planet_TOF_res, Spacecraft_r, Spacecraft_v, figNumber)
% Same plotting helper as above, but allows caller to choose figure number.

    figure(figNumber)
    clf
    set(gcf, 'Color', [0.06 0.06 0.06]);

    muSunCanonical = 1;
    TOF = linspace(Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res);

    Earth_orbit_r   = zeros(3, length(TOF));
    Mars_orbit_r    = zeros(3, length(TOF));
    Jupiter_orbit_r = zeros(3, length(TOF));
    SC_orbit_r      = zeros(3, length(TOF));

    Earth_orbit_v   = zeros(3, length(TOF));
    Mars_orbit_v    = zeros(3, length(TOF));
    Jupiter_orbit_v = zeros(3, length(TOF));
    SC_orbit_v      = zeros(3, length(TOF));

    for k = 1:length(TOF)
        [Earth_orbit_r(:,k), Earth_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Earth_r, Earth_v);
        [Mars_orbit_r(:,k),  Mars_orbit_v(:,k)]  = universalTOF(muSunCanonical, TOF(k), Mars_r, Mars_v);
        [Jupiter_orbit_r(:,k), Jupiter_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Jupiter_r, Jupiter_v);

        if nargin >= 11 && ~isempty(Spacecraft_r) && ~isempty(Spacecraft_v)
            [SC_orbit_r(:,k), SC_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Spacecraft_r, Spacecraft_v);
        end
    end

    makeThreeBodyPlot(Earth_r, Mars_r, Jupiter_r, Earth_orbit_r, Mars_orbit_r, ...
        Jupiter_orbit_r, SC_orbit_r, Spacecraft_r, figNumber);

    Earth_final_r = Earth_orbit_r(:,end); Earth_final_v = Earth_orbit_v(:,end);
    Mars_final_r = Mars_orbit_r(:,end);   Mars_final_v = Mars_orbit_v(:,end);
    Jupiter_final_r = Jupiter_orbit_r(:,end); Jupiter_final_v = Jupiter_orbit_v(:,end);

    if any(SC_orbit_r(:))
        SC_final_r = SC_orbit_r(:,end);
        SC_final_v = SC_orbit_v(:,end);
    else
        SC_final_r = [];
        SC_final_v = [];
    end
end
function makeThreeBodyPlot(Earth_r, Mars_r, Jupiter_r, Earth_orbit_r, Mars_orbit_r, ...
    Jupiter_orbit_r, SC_orbit_r, Spacecraft_r, figNumber)
% Shared plotting utility used by both local plotting functions.

    figure(figNumber)
    ax = axes('Parent', gcf, ...
        'Color', [0.03 0.03 0.03], ...
        'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
        'ZColor', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
    hold(ax, 'on')

    sunColor = [1.0, 0.85, 0.0];
    earthColor = [0.2, 0.6, 1.0];
    marsColor = [1.0, 0.5, 0.2];
    jupiterColor = [1.0, 0.85, 0.7];
    orbitColor = [0.7, 0.7, 0.7];
    scColor = [1.0, 0.2, 0.2];
    finalMarkerColor = [0.9, 0.9, 0.9];

    plot3(ax, 0, 0, 0, 'o', 'MarkerFaceColor', sunColor, 'MarkerEdgeColor', sunColor, 'MarkerSize', 8, 'DisplayName', 'Sun')

    plot3(ax, Earth_r(1), Earth_r(2), Earth_r(3), 'o', 'MarkerFaceColor', earthColor, 'MarkerEdgeColor', earthColor, 'DisplayName', 'Earth')
    plot3(ax, Mars_r(1), Mars_r(2), Mars_r(3), 'o', 'MarkerFaceColor', marsColor, 'MarkerEdgeColor', marsColor, 'DisplayName', 'Mars')
    plot3(ax, Jupiter_r(1), Jupiter_r(2), Jupiter_r(3), 'o', 'MarkerFaceColor', jupiterColor, 'MarkerEdgeColor', jupiterColor, 'DisplayName', 'Jupiter')

    plot3(ax, Earth_orbit_r(1,:), Earth_orbit_r(2,:), Earth_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')
    plot3(ax, Mars_orbit_r(1,:), Mars_orbit_r(2,:), Mars_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')
    plot3(ax, Jupiter_orbit_r(1,:), Jupiter_orbit_r(2,:), Jupiter_orbit_r(3,:), 'Color', orbitColor, 'LineWidth', 1, 'HandleVisibility', 'off')

    if any(SC_orbit_r(:))
        plot3(ax, SC_orbit_r(1,:), SC_orbit_r(2,:), SC_orbit_r(3,:), '-', 'Color', scColor, 'LineWidth', 1.5, 'DisplayName', 'Spacecraft')
        plot3(ax, Spacecraft_r(1), Spacecraft_r(2), Spacecraft_r(3), 's', 'MarkerFaceColor', scColor, 'MarkerEdgeColor', scColor, 'DisplayName', 'SC Start')
        SC_final_r = SC_orbit_r(:,end);
        plot3(ax, SC_final_r(1), SC_final_r(2), SC_final_r(3), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'MarkerSize', 8, 'LineWidth', 1.2, 'DisplayName', 'SC Final')
    end

    plot3(ax, Earth_orbit_r(1,end), Earth_orbit_r(2,end), Earth_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')
    plot3(ax, Mars_orbit_r(1,end), Mars_orbit_r(2,end), Mars_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')
    plot3(ax, Jupiter_orbit_r(1,end), Jupiter_orbit_r(2,end), Jupiter_orbit_r(3,end), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', finalMarkerColor, 'HandleVisibility', 'off')

    legend(ax, 'TextColor', [0.95 0.95 0.95], 'Color', [0.12 0.12 0.12])
    title(ax, 'Solar System (Earth, Mars, Jupiter)', 'Color', [0.95 0.95 0.95])
    xlabel(ax, 'x (AU)', 'Color', [0.95 0.95 0.95])
    ylabel(ax, 'y (AU)', 'Color', [0.95 0.95 0.95])
    zlabel(ax, 'z (AU)', 'Color', [0.95 0.95 0.95])
    axis(ax, 'vis3d')
    daspect(ax, [1 1 1])
    grid(ax, 'on')
    view(ax, 3)
    hold(ax, 'off')
end
