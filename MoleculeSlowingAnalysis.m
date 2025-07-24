colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
    [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
    [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

moleculeName = "SrF";
molecule = def_molecule(moleculeName);
vib_repump = 0; %1 if vibrational repump used, 0 otherwise
bfield = "5"; %static b-field used
num_lasers_used = "2"; %number of lasers used in WLS. 1 for WLS no repump, 2 for WLS w/repump OR WLS no repump + push, 3 for WLS w/repump + push
data_timestamp = "20250311_2259";
base_folder = "saveData/";
if vib_repump == 0
    dataFolder = strcat(base_folder, moleculeName, "WhiteLightSlowingbFieldSettingStaticBFieldGauss", bfield, "ForceTwoDNumLasers", num_lasers_used, "Date", data_timestamp);
else
    dataFolder = strcat(base_folder, moleculeName, "WhiteLightSlowingVibRepumpbFieldSettingStaticBFieldGauss", bfield, "ForceTwoDNumLasers", num_lasers_used, "Date", data_timestamp);
end

pos = {'0.01'};

for i=1:length(pos)
    posForPlot(i) = str2num(pos{i});
end

testFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{1},'MMSameDir.dat');
testData = readtable(testFile);

LongSpeeds = zeros(size(testData, 1),length(pos));
azs = zeros(size(testData,1),length(pos));
excitedPop = zeros(size(testData,1),length(pos));
for i=1:length(pos)
    currFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{i},'MMSameDir.dat');
    currData = readtable(currFile);
    %vels = currData.Speed;
    LongSpeeds(:,i) = currData.LongSpeed;
    azs(:,i) = currData.az/1e3;
    excitedPop(:,i) = currData.PExc;
end


% Filter out data where LongSpeeds <= 0
validIdx = LongSpeeds > 0;
azs = azs(validIdx);
LongSpeeds = LongSpeeds(validIdx);

% Sort the data for proper interpolation
[LongSpeedsSorted, sortIdx] = sort(LongSpeeds);
azsSorted = azs(sortIdx);

% Perform 1D interpolation
LongSpeedsInterp = linspace(min(LongSpeedsSorted), max(LongSpeedsSorted), 200); % Fine grid
azsInterp = interp1(LongSpeedsSorted, azsSorted, LongSpeedsInterp, 'spline'); % Spline interpolation

% Plot the data points and the interpolation
figure(1);
%plot(LongSpeeds, azs, 'o', 'DisplayName', 'Data Points'); % Data points
%hold on;
plot(LongSpeedsInterp, azsInterp, '-', 'LineWidth', 2, 'DisplayName', '1D Interpolation'); % Interpolated curve
%hold off;

% Customize the plot
xlabel('v_z (m/s)'); % x-axis label
ylabel('a_z(v_z) (mm/ms^2)'); % y-axis label
title('White Light Slowing of SrF molecules');
%legend('show', 'Location','northeastoutside');
grid on;


%now, given an a_z(v_z) curve, and initial condition at beam source,
%propagate molecule motion to the MOT region and see where it ends up.

slowingAccel = @(v) interp1(LongSpeedsSorted, azsSorted, v, 'spline');
%define anonymous function for computing slowingAccel at some speed v. Note v_z is in units of m/s and a_z in units of mm/ms^2

odefun = @(t,p) [p(2); slowingAccel(p(2))];
%define system. [dx/dt=v; dv/dt=slowingAccel(v)]


%initial parameters
tspan = [0 20]; %solve from t=0 ms to t=20 ms.
initialPos = 0.0; %in mm
initialVel = 140; %in mm/ms (equivalently m/s)
slowingLength = 700.0; % in mm
y0 = [initialPos; initialVel]; %vector of initial conditions


[t, y] = ode45(odefun, tspan, y0); %solve system



% Plot z vs t trajectory
figure(2);
h1 = plot(t, y(:,1), 'b', 'LineWidth', 1.5); % Plot position vs time
hold on;

% Find the first time when z crosses 1000 mm
z_target = slowingLength; % Target position in mm
crossingIdx = find(y(:,1) >= z_target, 1, 'first'); % First index where z >= 1000 mm

if ~isempty(crossingIdx)
    t_cross = t(crossingIdx); % Time when z first reaches 1000 mm

    % Plot horizontal and vertical lines
    yline(z_target, 'k--', 'LineWidth', 1.2); 
    xline(t_cross, 'r--', 'LineWidth', 1.2);

    % Mark the intersection point
    plot(t_cross, z_target, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
end

% Labels and title
xlabel('t (ms)', 'FontSize', 12);
ylabel('z (mm)', 'FontSize', 12);
title(['SrF Molecule Slowing Trajectory for v_0 = ' num2str(initialVel) ' m/s'], 'FontSize', 14);

% Add legend
legend(h1, 'Trajectory', 'Location', 'best');

hold off;


% Plot v vs t trajectory
figure(3);
h2 = plot(t, y(:,2), 'b', 'LineWidth', 1.5); % Plot velocity vs time
hold on;

% Find the time when v crosses zero
zeroCrossingIdx = find(y(:,2) <= 0, 1, 'first'); % First index where v <= 0
if ~isempty(zeroCrossingIdx)
    t_cross2 = t(zeroCrossingIdx); % Time when v crosses zero

    % Plot horizontal and vertical lines
    yline(0, 'k--', 'LineWidth', 1.2); 
    xline(t_cross2, 'r--', 'LineWidth', 1.2);

    % Mark the intersection point
    plot(t_cross2, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
end

% Labels and title
xlabel('t (ms)', 'FontSize', 12);
ylabel('v_z (m/s)', 'FontSize', 12);
title(['SrF Molecule Slowing Trajectory for v_0 = ' num2str(initialVel) ' m/s'], 'FontSize', 14);

% Add legend
legend(h2, 'Trajectory', 'Location', 'best');

hold off;

%plot v_z vs z trajectory

% Define target position z_target
z_target = slowingLength; % Example value, replace with your desired z_target

% Plot v_z vs z trajectory
figure(4);
h3 = plot(y(:,1), y(:,2), 'b', 'LineWidth', 1.5); % Plot velocity vs position
hold on;

% Find the first index where z >= z_target
crossingIdx = find(y(:,1) >= z_target, 1, 'first'); % First index where z >= z_target

if ~isempty(crossingIdx)
    v_at_z_target = y(crossingIdx,2); % Velocity when z first reaches z_target

    % Plot horizontal line at v_z where z = z_target
    yline(v_at_z_target, 'k--', 'LineWidth', 1.2); % Dashed black horizontal line

    % Plot vertical line at z = z_target
    xline(z_target, 'r--', 'LineWidth', 1.2); % Dashed red vertical line

    % Mark the intersection point
    plot(z_target, v_at_z_target, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
end

% Labels and title
xlabel('z (mm)', 'FontSize', 12);
ylabel('v_z (m/s)', 'FontSize', 12);
title(['SrF Molecule Slowing Trajectory for v_0 = ' num2str(initialVel) ' m/s'], 'FontSize', 14);

% Add legend
legend(h3, 'Trajectory', 'Location', 'best');

hold off;

