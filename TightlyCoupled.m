clc
close all

imuFs = 160;
gpsFs = 1;

startPosition = [42.2825 -72.3430 53.0352]; % latitude longitude altitude

secondsToSimulate = 10;

%% Loading ground truth trajectory
load LoggedQuadcopter.mat trajData;
trajOrient = trajData.Orientation;
trajVel = trajData.Velocity;
trajPos = trajData.Position;
trajAcc = trajData.Acceleration;
trajAngVel = trajData.AngularVelocity;

%% Setting up sensors

rng(1);

% GPS sattelite constalation data
data = rinexread("GODS00USA_R_20211750000_01D_GN.rnx");
[~,idx] = unique(data.GPS.SatelliteID);
navmsg = data.GPS(idx,:);
t0 = navmsg.Time(1);

imu = imuSensor('accel-gyro', 'SampleRate', imuFs);
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

%% Setting up filter

load TunedValuesTightlyCoupled.mat;

gnss = exampleHelperINSGNSS;
gnss.ReferenceLocation = startPosition;
accel = insAccelerometer;
gyro = insGyroscope;

filt = insEKF(accel,gyro,gnss,insMotionPose);

filt.AdditiveProcessNoise = TunedValues.AdditiveProcessNoise;

stateparts(filt,"Orientation",compact(trajOrient(1)));
stateparts(filt,"Position",trajPos(1,:));
stateparts(filt,"Velocity",trajVel(1,:));
stateparts(filt,"Accelerometer_Bias",imu.Accelerometer.ConstantBias ./ imuFs);
stateparts(filt,"Gyroscope_Bias",imu.Gyroscope.ConstantBias ./ imuFs);

% statecovparts(filt,"Orientation", (eye(4) * 1e-3) );
filt.StateCovariance = TunedValues.StateCovariance;

%% Main loop

assert(rem(imuFs, gpsFs) == 0, "GPS sampling rate must be an integer factor of IMU sampling rate."); % for simplicity
imuSamplesPerGPS = (imuFs / gpsFs);

numsamples = secondsToSimulate * imuFs;
loopBound = floor(numsamples / imuFs) * imuFs;

% Data for final metric computation.
pqorient = quaternion.zeros(loopBound, 1);
pqpos = zeros(loopBound, 3);

fcnt = 1;

while fcnt <= loopBound
    
    % |predict| loop at IMU update frequency.
    for ff = 1:imuSamplesPerGPS

        [accelMeas, gyroMeas] = imu(trajAcc(fcnt, :), trajAngVel(fcnt, :), trajOrient(fcnt)); % simulating imu data

        fuse(filt,accel,accelMeas, TunedValues.AccelerometerNoise);
        fuse(filt,gyro,gyroMeas, TunedValues.GyroscopeNoise);

        pqpos(fcnt,:) = stateparts(filt,"Position");
        pqorient(fcnt,:) = quaternion(stateparts(filt,"Orientation"));

        predict(filt,1/imuFs);

        fcnt = fcnt + 1;
    end

    t = t0 + fcnt*seconds(1/imuFs);
    [satPos,satVel,satIDs] = gnssconstellation(t,"RINEXData",navmsg);
    [az,el,vis] = lookangles(trajPos(fcnt, :),satPos);
    [p,pdot] = pseudoranges(trajPos(fcnt, :),satPos(vis,:),trajVel(fcnt, :),satVel(vis,:));
    z = [p; pdot];

    gnss.SatellitePosition = satPos(vis,:);
    gnss.SatelliteVelocity = satVel(vis,:);

    pseudorangeNoise = TunedValues.exampleHelperINSGNSSNoise(1);
    pseudorangeRateNoise = TunedValues.exampleHelperINSGNSSNoise(2);

    fuse(filt,gnss,z, [pseudorangeNoise*ones(1,numel(p)),pseudorangeRateNoise*ones(1,numel(pdot))]);

end

positionError = pqpos(1:loopBound, :) - trajPos(1:loopBound, :);
orientationError = rad2deg(dist(pqorient(1:loopBound), trajOrient(1:loopBound)));

%% Displaying results
fprintf('\n\nEnd-to-End Simulation Position RMS Error\n');

msep = sqrt(mean(positionError.^2));
fprintf('\tX: %.2f , Y: %.2f, Z: %.2f   (meters)\n\n', msep(1), ...
    msep(2), msep(3));

fprintf('End-to-End Quaternion Distance RMS Error (degrees) \n');

fprintf('\t%.2f (degrees)\n\n', sqrt(mean(orientationError.^2)));

figure();
tiledlayout('flow');

time = 0:1 / imuFs:loopBound * 1 / imuFs;
time = time(:,1:loopBound)';

nexttile;
plot(time, [positionError(:, 1), positionError(:, 2), positionError(:, 3)]);
title("Ошибка определение положения");
grid on;
xlabel('время');
ylabel("м.");
legend("x", "y", "z");

nexttile;
plot(time, orientationError);
title("Ошибка определенния ориентации");
grid on;
xlabel('время');
ylabel("град.");
