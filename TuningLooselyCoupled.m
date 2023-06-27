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

gps = gpsSensor('UpdateRate', gpsFs);
gps.ReferenceLocation = startPosition;
gps.DecayFactor = 0.5;
gps.HorizontalPositionAccuracy = 1.6;
gps.VerticalPositionAccuracy =  1.6;
gps.VelocityAccuracy = 0.1;

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

%% Generating sensor readings

assert(rem(imuFs, gpsFs) == 0, "GPS sampling rate must be an integer factor of IMU sampling rate."); % for simplicity
imuSamplesPerGPS = (imuFs / gpsFs);

numsamples = secondsToSimulate * imuFs;
loopBound = floor(numsamples / imuFs) * imuFs;

% Sensor readings
Accelerometer = zeros(loopBound, 3);
Gyroscope = zeros(loopBound, 3);
GPS = zeros(loopBound/imuSamplesPerGPS, 6);

fcnt = 1;
while fcnt <= loopBound

    for ff = 1:imuSamplesPerGPS
        [accel, gyro] = imu(trajAcc(fcnt, :), trajAngVel(fcnt, :), trajOrient(fcnt));
        Accelerometer(fcnt,:) = accel;
        Gyroscope(fcnt,:) = gyro;
        fcnt = fcnt + 1;
    end

    [lla, gpsvel] = gps(trajPos(fcnt, :), trajVel(fcnt, :));
    GPS((fcnt-1)/imuSamplesPerGPS,:) = [lla, gpsvel];
end

%% Packing ground-truth and sensor readings
imuData =  timetable(Accelerometer, Gyroscope, 'SampleRate', imuFs);
gpsData = timetable(GPS, 'SampleRate', gpsFs);
sensorData = synchronize(imuData,gpsData);

Orientation = trajOrient(1:loopBound,:);
Position = trajPos(1:loopBound,:);
groundTruth = timetable(Orientation, Position, 'SampleRate', imuFs);

%% Setting up filter

load TunedValuesLooselyCoupled.mat;

gpssens = insGPS;
gpssens.ReferenceLocation = startPosition;
accel = insAccelerometer;
gyro = insGyroscope;

filt = insEKF(accel,gyro,gpssens,insMotionPose);

filt.AdditiveProcessNoise = TunedValues.AdditiveProcessNoise;

stateparts(filt,"Orientation",compact(trajOrient(1)));
stateparts(filt,"Position",trajPos(1,:));
stateparts(filt,"Velocity",trajVel(1,:));
stateparts(filt,"Accelerometer_Bias",imu.Accelerometer.ConstantBias ./ imuFs);
stateparts(filt,"Gyroscope_Bias",imu.Gyroscope.ConstantBias ./ imuFs);

filt.StateCovariance = TunedValues.StateCovariance;

%% Tuning
config = tunerconfig(filt,'MaxIterations',10);
config.StepForward = 5;
config.StepBackward = 0.99;

measNoise = tunernoise(filt);
measNoise.AccelerometerNoise = TunedValues.AccelerometerNoise;
measNoise.GyroscopeNoise = TunedValues.GyroscopeNoise;
measNoise.GPSNoise = TunedValues.GPSNoise;

tunedParams = tune(filt,measNoise,sensorData,groundTruth,config);
