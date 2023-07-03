clc;
close all;

imuFs = 160;
gpsFs = 1;

startPosition = [42.2825 -72.3430 53.0352];
magneticField = [19.5281 -5.0741 48.0067];

secondsToSimulate = 50;

%% Loading ground truth trajectory
load LoggedQuadcopter.mat trajData;
trajOrient = trajData.Orientation;
trajVel = trajData.Velocity;
trajPos = trajData.Position;
trajAcc = trajData.Acceleration;
trajAngVel = trajData.AngularVelocity;

%% Generating sensor readings

rng(1);

gps = gpsSensor('UpdateRate', gpsFs);
gps.ReferenceLocation = startPosition;
gps.DecayFactor = 0.5;              % Random walk noise parameter
gps.HorizontalPositionAccuracy = 1.6;
gps.VerticalPositionAccuracy =  1.6;
gps.VelocityAccuracy = 0.1;

imu = imuSensor('accel-gyro-mag', 'SampleRate', imuFs);
imu.MagneticField = magneticField;

% Accelerometer
imu.Accelerometer.MeasurementRange =  19.6133;
imu.Accelerometer.Resolution = 0.0023928;
imu.Accelerometer.ConstantBias = 0.19;
imu.Accelerometer.NoiseDensity = 0.0012356;

% Gyroscope
imu.Gyroscope.MeasurementRange = deg2rad(250);
imu.Gyroscope.Resolution = deg2rad(0.0625);
imu.Gyroscope.ConstantBias = deg2rad(3.125);
imu.Gyroscope.AxesMisalignment = 1.5;
imu.Gyroscope.NoiseDensity = deg2rad(0.025);

% Magnetometer
imu.Magnetometer.MeasurementRange = 1000;
imu.Magnetometer.Resolution = 0.1;
imu.Magnetometer.ConstantBias = 100;
imu.Magnetometer.NoiseDensity = 0.3 / sqrt(50);

%% Setting up filter

load filterParams.mat;

margFilter = insfilterMARG;
margFilter.IMUSampleRate = imuFs;
margFilter.ReferenceLocation = startPosition;

initstate = zeros(22, 1);
initstate(1:4) = compact(meanrot(trajOrient(1:100)));
initstate(5:7) = mean(trajPos(1:100, :), 1);
initstate(8:10) = mean(trajVel(1:100, :), 1);
initstate(11:13) =  imu.Gyroscope.ConstantBias ./ imuFs;
initstate(14:16) =  imu.Accelerometer.ConstantBias ./ imuFs;
initstate(17:19) =  imu.MagneticField;
initstate(20:22) = imu.Magnetometer.ConstantBias;

margFilter.State = initstate;

% Process noises
margFilter.AccelerometerBiasNoise =  filterParams.AccelerometerBiasNoise;
margFilter.AccelerometerNoise = filterParams.AccelerometerNoise;
margFilter.GyroscopeBiasNoise = filterParams.GyroscopeBiasNoise;
margFilter.GyroscopeNoise =  filterParams.GyroscopeNoise;
margFilter.MagnetometerBiasNoise = filterParams.MagnetometerBiasNoise;
margFilter.GeomagneticVectorNoise = filterParams.GeomagneticVectorNoise;

% Initial error covariance
%margFilter.StateCovariance = 1e-9 * eye(22);
margFilter.StateCovariance = filterParams.StateCovariance;

%% Main Loop

assert(rem(imuFs, gpsFs) == 0, "GPS sampling rate must be an integer factor of IMU sampling rate."); % for simplicity
imuSamplesPerGPS = (imuFs / gpsFs);
numsamples = secondsToSimulate * imuFs;

loopBound = floor(numsamples);
loopBound = floor(loopBound / imuFs) * imuFs; % ensure enough IMU Samples

% Log data for final metric computation.
pqorient = quaternion.zeros(loopBound, 1);
pqpos = zeros(loopBound, 3);

fcnt = 1;

while fcnt <= loopBound
    % |predict| loop at IMU update frequency.
    for ff = 1:imuSamplesPerGPS
        % Simulate the IMU data from the current pose.
        [accel, gyro, mag] = imu(trajAcc(fcnt, :), trajAngVel(fcnt, :), trajOrient(fcnt));

        % Use the |predict| method to estimate the filter state based
        % on the simulated accelerometer and gyroscope signals.
        predict(margFilter, accel, gyro);

        % Acquire the current estimate of the filter states.
        [fusedPos, fusedOrient] = pose(margFilter);

        % Save the position and orientation for post processing.
        pqorient(fcnt) = fusedOrient;
        pqpos(fcnt, :) = fusedPos;

        fcnt = fcnt + 1;
    end

    % This next step happens at the GPS sample rate.
    % Simulate the GPS output based on the current pose.
    [lla, gpsvel] = gps(trajPos(fcnt, :), trajVel(fcnt, :));

    % Correct the filter states based on the GPS data and magnetic
    % field measurements.
    fusegps(margFilter, lla, filterParams.gpsPositionMesurmentNoise, gpsvel, filterParams.gpsVelocityMesurmentNoise);
    fusemag(margFilter, mag, filterParams.MagnetometerMesurmentNoise);
end

positionError = pqpos(1:loopBound, :) - trajPos(1:loopBound, :);
orientationError = rad2deg(dist(pqorient(1:loopBound), trajOrient(1:loopBound)));

%% Displaying results
fprintf('\n\nEnd-to-End Simulation Position RMS Error\n');

msep = sqrt(mean(positionError.^2));
fprintf('\tX: %.2f , Y: %.2f, Z: %.2f   (meters)\n\n', msep(1), msep(2), msep(3));

fprintf('End-to-End Quaternion Distance RMS Error (degrees) \n');

fprintf('\t%.2f (degrees)\n\n', sqrt(mean(orientationError.^2)));

figure();
tiledlayout('flow');

time = 0:1 / imuFs:loopBound * 1 / imuFs;
time = time(:,1:loopBound)';

nexttile;
plot(time, [positionError(:, 1), positionError(:, 2), positionError(:, 3)]);
title("Position error");
grid on;
xlabel('time');
legend("x", "y", "z");

nexttile;
plot(time, orientationError);
title("Orientation error");
grid on;
xlabel('time');
