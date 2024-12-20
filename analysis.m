colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
    [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
    [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

simTrapDynamics =1; %change to 1 if you want to simulate particle trajectory, with random photon scatter, to get 'true' size and temperature
moleculeName = "SrF";
molecule = def_molecule(moleculeName);
vib_repump = 0; %1 if vibrational repump used, 0 otherwise
bfield_grad = "8.8"; %b-field gradient used
num_lasers_used = "5"; %number of lasers used in MOT
data_timestamp = "20240910_2054";
base_folder = "saveData/";
if vib_repump == 0
    dataFolder = strcat(base_folder, moleculeName, "RedMOTbFieldSettingThreeDBGradGPerCM", bfield_grad, "ForceThreeDNumLasers", num_lasers_used, "Date", data_timestamp);
else
    dataFolder = strcat(base_folder, moleculeName, "RedMOTVibRepumpbFieldSettingThreeDBGradGPerCM", bfield_grad, "ForceThreeDNumLasers", num_lasers_used, "Date", data_timestamp);
end
%dataFolder = 'saveData/SrFRedMOTbFieldSettingThreeDBGradGPerCM8.0ForceThreeDNumLasers5Date20240911_0927';

% pos = {'0.5','1.0','1.5','2.0','2.5','3.0'};
%pos = {'0.5','1.5','3.0','4.5','6.0','7.5'};
%pos = {'1', '2', '3', '5', '7', '9', '11', '14', '17'};
pos = {'0.5', '1.5', '3.0', '4.5', '6.0', '7.5', '9.0', '10.5', '12.0', '13.5', '15.0'};
%pos = {'0.5', '1.0', '1.5', '2.0', '3.0', '4.0', '5.0', '6.0', '7.5', '9.0', '10.5', '12.0', '14.0', '16.0'};
for i=1:length(pos)
    posForPlot(i) = str2num(pos{i});
end
posInMM = posForPlot.*1e-3;
testFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{1},'MMSameDir.dat');
testData = readtable(testFile);

accelsInVelDirection = zeros(size(testData,1),length(pos));
accelsInPosDirection = zeros(size(testData,1),length(pos));
excitedPop = zeros(size(testData,1),length(pos));
for i=1:length(pos)
    currFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{i},'MMSameDir.dat');
    currData = readtable(currFile);
    vels = currData.Speed;
    accelsInVelDirection(:,i) = currData.av;
    accelsInPosDirection(:,i) = currData.ar;
    excitedPop(:,i) = currData.PExc;
end

%reverse sign of accel for negative velocities
for i=1:length(vels)
    if vels(i)<0
        accelsInVelDirection(i,:) = accelsInVelDirection(i,:).*-1;
    end
end

%add opposite positions
accelsInVelDirectionFull = [-flipud(fliplr(accelsInVelDirection)),accelsInVelDirection];
excitedPopFull = [flipud(fliplr(excitedPop)),excitedPop];

sortedPos = sort([-posForPlot,posForPlot]);
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v);%in mm,ms units


%simulate capture with linear interpolation (spline acts up here for some
%reason.  Probably not a big deal to use linear)
diffEqVals = @(t,p) [p(2);oneAxisAccel(p(1),p(2))];
vsToTryForCapture = [1:0.1:30]; %step size is 0.1 m/s
for i=1:length(vsToTryForCapture)
    currV = vsToTryForCapture(i);
    [ts2,ps2] = ode23(diffEqVals,[0 20],[min(sortedPos);currV]);
    %[ts2,ps2] = ode23(diffEqVals,[0 50],[min(sortedPos);currV]); %try 50 ms capture time
    %[ts2,ps2] = ode23(diffEqVals,[0 1e-1],[-0*1e-3;currV]);
    if isnan(ps2(end,1))
        break;
    end
end

capVel = currV-0.2;

%{
% Plot multiple phase space trajectories on same plot
figure;
hold on; % This ensures all plots are on the same figure

% Initialize an array to store plot handles
plotHandles = gobjects(currV, 1); % Preallocate for better performance

% Initialize a cell array to store legend labels
legendLabels = cell(currV, 1);

% Loop from currV down to 1, in steps of -3
for v = currV:-3:1
    % Solve the ODE with the current velocity
    [ts2, ps2] = ode23(diffEqVals, [0 50], [min(sortedPos); v]);
    
    % Plot the entire trajectory and store the plot handle
    plotHandles(v) = plot(ps2(:, 1), ps2(:, 2), 'LineWidth',2);
    
    % Store the legend label for the current plot
    legendLabels{v} = ['v_0 = ' num2str(v) ' m/s'];
end

% Set labels and title
xlabel('d (mm)');
ylabel('v_{d} (m/s)');
title('Particle Trajectories for Different Initial Velocities');

% Create the legend
legend(plotHandles(currV:-3:1), legendLabels(currV:-3:1), 'Location', 'best');
%saveas(gcf, 'phase_space_MOT_capture_highDPI.png');
%print(gcf, 'phase_space_MOT_capture_highDPI', '-dpng', '-r300'); % 300 DPI

hold off;
%}
%{
print('HighDPIPlot', '-dpng', '-r300'); % Saves as 'HighDPIPlot.png'


%}


%plot the highest value of v for which we have capture (evolve out to t=100 ms)
[ts2,ps2] = ode23(diffEqVals,[0 100],[min(sortedPos);currV-0.2]);
figure(1);
plot(ps2(:,1),ps2(:,2),'LineWidth',2);
xlabel('d (mm)', fontsize=12);
ylabel('v_{d} (m/s)', fontsize=12);
title(strcat('Molecule Trajectory for v_{Cap}=',num2str(currV-0.2),' m/s'), fontsize=14)
%write to csv
%writematrix(ps2, strcat('MoleculeMOTTrials/CaptureVelocityTrajectory', data_timestamp, '.csv'));

%{
%plot lower velocity trajectories as well, and for 50 ms instead of 20 ms
[ts2,ps2] = ode23(diffEqVals,[0 50],[min(sortedPos);currV-5]);
figure(1);

% Filter the time points between t = 4 and t = 50
timeIndices = (ts2 >= 4) & (ts2 <= 50);

% Plot the filtered data
plot(ps2(timeIndices, 1), ps2(timeIndices, 2));
xlabel('z(mm)');
ylabel('v (m/s)')
title('Particle Trajectory from t = 4 to t = 50');

%plot trajectory for a speed lower than capture velocity
%{
plot(ps2(:,1),ps2(:,2));
xlabel('x(mm)');
ylabel('v (m/s)')
title(strcat('particle trajectory for v_{Cap}=',num2str(currV-8),' m/s'))
%}
%}


%now plot a(x,v) heat map
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v,'spline');%in mm,ms units


%make heat map
vsForHeatMap = [min(vels):.1:max(vels)];
xsForHeatMap = [min(sortedPos):.1:max(sortedPos)];
for i=1:length(vsForHeatMap)
    for j=1:length(xsForHeatMap)
        heatMap(i,j) = oneAxisAccel(xsForHeatMap(j),vsForHeatMap(i));
    end
end
figure(2);
imagesc(xsForHeatMap,vsForHeatMap,heatMap)
colorbar
xlabel('d (mm)');
ylabel('v_{d} (m/s)')
title('Heat Map of MOT Dynamics, CaF molecule w/Vibrational Repump')
xlim([min(sortedPos) max(sortedPos)])
h=colorbar;
h.Title.String = "a_{d} (mm/ms^{2})"

%make x,v plots from matrix
vMaxToInt = 4;
xMaxToInt = 7;
[~,minCol] = min(abs(xsForHeatMap+xMaxToInt));
[~,maxCol] = min(abs(xsForHeatMap-xMaxToInt));
[~,minRow] = min(abs(vsForHeatMap+vMaxToInt));
[~,maxRow] = min(abs(vsForHeatMap-vMaxToInt));
accVsPosForPlot = trapz(vsForHeatMap(minRow:maxRow),heatMap(minRow:maxRow,:),1)./(2*vMaxToInt);
accVsVelForPlot = trapz(xsForHeatMap(minCol:maxCol),heatMap(:,minCol:maxCol),2)./(2*xMaxToInt);
figure(3);
hold all;
plot(xsForHeatMap,accVsPosForPlot,'Linewidth',2);
xlabel('d (mm)');
ylabel('a_{d} (mm/ms^{2})')
%xlim([0 15])
%ylim([-1 0.5])
title('Acceleration vs Position, SrF MOT')
figure(4);
hold all;
plot(vsForHeatMap,accVsVelForPlot,'LineWidth',2);
xlabel('v_{d} (m/s)');
ylabel('a_{d} (mm/ms^{2})')
%xlim([0 20])
%ylim([-5 1])
title('Acceleration vs Velocity, SrF MOT')

%write to csv
%writematrix([xsForHeatMap; accVsPosForPlot]', strcat('MoleculeMOTTrials/AccelVsPos', data_timestamp, '.csv'));
%writematrix([vsForHeatMap; accVsVelForPlot']', strcat('MoleculeMOTTrials/AccelVsVel', data_timestamp, '.csv'));


[~,minCol] = min(abs(sortedPos+xMaxToInt));
[~,maxCol] = min(abs(sortedPos-xMaxToInt));
[~,minRow] = min(abs(vels+vMaxToInt));
[~,maxRow] = min(abs(vels-vMaxToInt));
meanExcPop = mean(mean(excitedPopFull(minRow:maxRow,minCol:maxCol)));


%{
%get temperature and rms size
maxTime=40;
gam = molecule.gam;
kA = molecule.kA;
kRepump = molecule.kRepump;
mass = molecule.mass;
hbar = 1.05e-34;
kb = 1.38e-23;
if simTrapDynamics==1
    %     scatterRateData = dlmread(strcat(folder,'forceVsVCenter.dat'));
    %     scatterRate = scatterRateData(11,3).*gamSrF;
    scatterRate = meanExcPop*gam*1e-3;%in 1/ms
    tKick = 1/scatterRate;
    v(1) = 1;%mm/ms
    r(1) = 1;%mm
    vKick = hbar*kA/mass;
    for i=1:round(maxTime/tKick)
        if mod(i,50)==0
            vKick = hbar * kRepump / mass;
        else
            vKick = hbar * kA / mass;
        end
        randPhi1 = 2*pi*rand;
        randPhi2 = 2*pi*rand;
        v(i+1) = v(i)+vKick*(cos(randPhi1)+cos(randPhi2))+oneAxisAccel(r(i),v(i))*tKick;
        r(i+1) = r(i)+v(i)*tKick;
    end
    simTimes = 0:(tKick):maxTime;
    startTime = maxTime/2;
    startInd = maxTime/2/tKick;
    endInd= i;
    vSq=mean(v(startInd:endInd).^2);
    sigma=sqrt(mean(r(startInd:endInd).^2));
    temp = vSq*mass/kb;
end



%plot MOT trajectory including random photon scatter

figure;
hold on;

% Solve the differential equation and plot the trajectory (without photon scatter)
[ts2,ps2] = ode23(diffEqVals,[0 100],[1;1]); %1;1 indicates initial v and r
plot(ps2(:,1), ps2(:,2), 'LineWidth', 2, 'DisplayName', 'Without Photon Scatter'); % Thicker line and legend entry

% Use the trap dynamics algorithm to plot the trajectory with photon scatter
maxTime = 100;
scatterRate = meanExcPop * gam * 1e-3; % in 1/ms
tKick = 1 / scatterRate;
v(1) = 1; % mm/ms
r(1) = 1; % mm
vKick = hbar * kA / mass;

for i = 1:round(maxTime / tKick)
    if mod(i,50)==0
        vKick = hbar * kRepump / mass;
    else
        vKick = hbar * kA / mass;
    end
    randPhi1 = 2 * pi * rand;
    randPhi2 = 2 * pi * rand;
    v(i + 1) = v(i) + vKick * (cos(randPhi1) + cos(randPhi2)) + oneAxisAccel(r(i), v(i)) * tKick;
    r(i + 1) = r(i) + v(i) * tKick;
end

% Plot the trajectory with photon scatter
plot(r, v, '-', 'LineWidth', 2, 'DisplayName', 'With Photon Scatter'); % Thicker line and legend entry

% Add grid, labels, and title
grid on;
xlabel('d (mm)');
ylabel('v_{d} (mm/ms)');
title('Particle Trajectory for Captured Molecule in MOT for 100 ms');

% Create a legend
legend('Location', 'best');

%write to csv
%writematrix(ps2, strcat('MoleculeMOTTrials/MOTTrajectoryNoScatter', data_timestamp, '.csv'));
%writematrix([r; v]', strcat('MoleculeMOTTrials/MOTTrajectoryYesScatter', data_timestamp, '.csv'));

%}
% Save the figure as a 300 DPI image
%saveas(gcf, 'particle_trajectory_MOT.png');
%print(gcf, 'particle_trajectory_MOT', '-dpng', '-r300'); % 300 DPI
%}