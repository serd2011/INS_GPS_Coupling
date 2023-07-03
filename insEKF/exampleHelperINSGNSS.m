classdef exampleHelperINSGNSS < positioning.INSSensorModel
%EXAMPLEHELPERINSGNSS Model raw GNSS reading for sensor fusion
%   This sensor model converts the filter states to the predicted GNSS raw
%   measurements (pseudoranges and pseudorange rates). These are used in
%   the update phase of the Kalman filter. 
%
%   This model also tracks the GNSS receiver clock bias and drift. These
%   additional states are added to the Kalman filter and updated during the
%   prediction phase.

%   This class is for internal use only. It may be removed in the future.

%   Copyright 2021 The MathWorks, Inc.

    properties
        %SATELLITEPOSITION Visible satellite positions in ECEF coordinates
        %   Specify the visible satellite positions as an N-by-3 matrix in
        %   meters in the Earth-Centered-Earth-Fixed (ECEF) coordinate
        %   system. N is the number of visible satellites.
        SatellitePosition;
        %SATELLITEVELOCITY Visible satellite velocities in ECEF coordinates
        %   Specify the visible satellite velocities as an N-by-3 matrix in
        %   meters per second in the Earth-Centered-Earth-Fixed (ECEF)
        %   coordinate system. N is the number of visible satellites.
        SatelliteVelocity;

        %ReferenceLocation Local origin in geodetic coordinates
        %   Specify the origin of the local coordinate system as a
        %   3-element vector in geodetic coordinates (latitude, longitude,
        %   and altitude). Altitude is the height above the reference
        %   ellipsoid model, WGS84. The reference location is in [degrees
        %   degrees meters]. The default value is [0 0 0].
        ReferenceLocation = [0 0 0];

    end

    properties (Access = private)
        %LOSVector Line-of-sight vector from receiver to satellites.
        %   This is a private property used to avoid recalculating the
        %   line-of-sight vector for the measurement Jacobian matrix. It is
        %   an N-by-3 matrix of vectors. N is the number of visible
        %   satellites.
        LOSVector
    end

    methods
        function statesdot = stateTransition(sensor, filt, dt)
            clockDrift = stateparts(filt, sensor, "ClockDrift");
            statesdot = struct( ...
                "ClockBias", clockDrift * dt, ...
                "ClockDrift", 0);
        end

        function z = measurement(sensor, filt)
            refFrame = ...
                fusion.internal.frames.ReferenceFrame.getMathObject( ...
                filt.ReferenceFrame);
            % Convert state from NED to ECEF.
            lla0 = sensor.ReferenceLocation;
            
            lla = refFrame.frame2lla(stateparts(filt, "Position"), ...
                lla0);

            % Convert position state.
            recPos = fusion.internal.frames.lla2ecef(lla);

            % Convert velocity state.
            recVel = refFrame.frame2ecefv( ...
                stateparts(filt, "Velocity"), lla(1), lla(2));

            % Get pseudorange and pseudorange rate from filter state.
            satPos = sensor.SatellitePosition;
            satVel = sensor.SatelliteVelocity;
            [p, pdot, losVector] ...
                = nav.internal.gnss.calculatePseudoranges( ...
                satPos, satVel, recPos, recVel);
            % Add predicted clock bias to predicted pseudorange.
            p = p + stateparts(filt, sensor, "ClockBias");
            % Add predicted clock drift to predicted pseudorange rate.
            pdot = pdot + stateparts(filt, sensor, "ClockDrift");

            z = [p(:); pdot(:)];

            % Convert line-of-sight vector from ECEF to NED.
            R = [0 1 0; 1 0 0; 0 0 -1] * fusion.internal.frames.ecef2enurotmat(lla(1), lla(2));
            for ii = 1:size(losVector,1)
                losVector(ii,:) = R * (losVector(ii,:).');
            end
            % Store line-of-sight (LOS) vector to be used with the
            % measurement Jacobian matrix. 
            sensor.LOSVector = -losVector;
        end

        function dhdx = measurementJacobian(sensor, filt)
            losVector = sensor.LOSVector;
            info = stateinfo(filt);
            % Compute measurement Jacobian.
            numMeas = size(sensor.SatellitePosition, 1);
            dhdx = zeros(2*numMeas, numel(filt.State));
            dhdx(1:numMeas,info.Position) = losVector;
            dhdx(1:numMeas,info.exampleHelperINSGNSS_ClockBias) = 1;
            dhdx(numMeas+(1:numMeas),info.Velocity) = losVector;
            dhdx(numMeas+(1:numMeas),info.exampleHelperINSGNSS_ClockDrift) = 1;
        end
    end

    methods
        function s = sensorstates(filt, opts) %#ok<INUSD> 
            s = struct( ...
                "ClockBias", 0, ...
                "ClockDrift", 0);
        end
    end
end