function OrbitPointing
% OrbitPointing
% -------------------------
% Epoch and orbit setup
% -------------------------
y = 2025; m = 4; d = 24;            % date (year, month, day)
hu = 8; mi = 0; sec = 0;           % time (hours, minutes, seconds)
date0 = juliandate(y,m,d,hu,mi,sec);% epoch (Julian date)
sma = 9000;                         % semi-major axis [km]
tarray = [0 90 50 140];                    % array of true anomalies [deg]
iarray = [45 45 135 90];                  % inclination array [deg]
RAarray = [0 0 50 50];                   % RAAN array [deg]

% Ground station coordinates (can be multiple)
long = [31,22,-100,-50];            % longitudes [deg]
lat  = [29,-20,40,-5];              % latitudes [deg]

d1t = 70;                            % time step between frames [s]
minAng = 50;                         % minimum pointing angle threshold [deg]
meshF = 'marsSat1.stl';              % satellite mesh file for poseplot

warning('off')

% -------------------------
% Physical constants
% -------------------------
gmst20 = 4.89496121282306;         % initial GMST (rad)
mu = 398600.4418;                   % Earth gravitational parameter [km^3/s^2]
Re = 6378.137;                       % Earth radius [km]

% -------------------------
% Figure / Earth texture
% -------------------------
hff = figure('Color', 'w');
set(hff, 'Toolbar', 'none', 'menubar', 'none', 'name', 'control constellition', 'NumberTitle', 'off');
hold on; axis equal; view(30, 0);
[Xe, Ye, Ze] = sphere(100);
image_file = 'EARTH.jpg';
cdata = imread(image_file);                        % Earth texture image
% create textured Earth surface
earth = surf(Re*Xe, Re*Ye, -Re*Ze, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
if ~isempty(gmst20)
    hgx = hgtransform;                              % parent transform (to rotate the Earth)
    set(earth, 'Parent', hgx);
end
set(earth, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
hold on; rotate3d on; view(0,0);

% -------------------------
% Pre-allocate arrays for multiple satellites
% -------------------------
R0 = zeros(length(tarray), 3);    % initial position vectors
V0 = R0;                           % initial velocity vectors
rm1 = cell(size(tarray));          % current position (per sat)
vm1 = rm1;                         % current velocity (per sat)
eull = rm1;                        % Euler angles (per sat)
xend = rm1;                        % final attitude state (per sat)
sat = rm1;                         % handle(s) for satellite pose visualization
vec = rm1;                         % handle(s) for pointing vector visualization

% -------------------------
% Initialize satellites: compute initial R/V and plot initial orbit
% -------------------------
for j = 1:length(tarray)
    % COE2RV converts classical orbital elements to position & velocity
    [R0(j,:), V0(j,:)] = COE2RV(sma, 0, tarray(j), RAarray(j), iarray(j), 0);
    % plotorbit2 sets up orbit plot, initial attitude, and pointing vector
    [rm1{j}, vm1{j}, eull{j}, xend{j}, sat{j}, vec{j}] = plotorbit2(R0(j,:), V0(j,:), date0);
end

% Set view and rotate Earth to epoch GMST
view(90, 0);
set(hgx, 'Matrix', makehgtform('zrotate', gmst2(date0-2400000.5)));

% Expand axis limits a bit for nicer framing
axis_data = get(gca);
xmin = axis_data.XLim(1)*1.2; xmax = axis_data.XLim(2)*1.2;
ymin = axis_data.YLim(1)*1.2; ymax = axis_data.YLim(2)*1.2;
zmin = axis_data.ZLim(1)*1.2; zmax = axis_data.ZLim(2)*1.2;
xlim([xmin xmax]); ylim([ymin ymax]); zlim([zmin zmax]);
axis('off')

% -------------------------
% Main loop: advance time, propagate, update attitude/visuals
% -------------------------
i = 1;
while(ishandle(hff))
    % Advance epoch by d1t seconds
    date0 = date0 + d1t/86400;
    theta = gmst2(date0-2400000.5);    % GMST angle for Earth rotation

    % Update each satellite state and visualization
    for j = 1:length(tarray)
        [eull{j}, xend{j}, rm1{j}, vm1{j}, sat{j}, vec{j}] = oporbit2(rm1{j}, vm1{j}, i, d1t, xend{j}, date0, sat{j}, vec{j});
    end

    % rotate Earth to match current GMST
    set(hgx, 'Matrix', makehgtform('zrotate', theta));
    drawnow;

    % show epoch as UI text (note: creates a new uicontrol each iter — small perf cost)
    uicontrol('Style','text', 'String', ['epoch : ', string(datetime(date0,'convertfrom','juliandate'))], ...
        'Units','normalized', 'Position',[0.007 0.005 0.5 0.06], 'BackgroundColor','none', 'ForegroundColor','red', 'FontSize',12, 'fontweight','bold');

    i = i + 1;
end

close all

% -------------------------
% Satellite orbit propagation (two-body integrator wrapper)
% -------------------------
function [rm, vm, State] = satOrbit2(R0, V0, t0, tf)
    options = odeset('reltol', 1.e-8, 'abstol', 1.e-8);
    X0 = [R0'; V0'];
    [~, State] = ode45(@(t,y) dEOM2(t,y), t0:1:tf, X0, options);
    rm = [State(end,1) State(end,2) State(end,3)];
    vm = [State(end,4) State(end,5) State(end,6)];
end

% -------------------------
% Convert lat/lon to ECEF (using GMST to rotate longitude)
% -------------------------
function R = Rmovee2(date, lat, thi)
    lon1g = gmst2(date-2400000.5)*(180/pi) + thi;   % Greenwich lon + local lon
    R = 6378 * [(cosd(lon1g).*cosd(lat))', (sind(lon1g).*cosd(lat))', sind(lat)'];
end

% -------------------------
% Two-body EOM
% -------------------------
function dydt = dEOM2(~, y)
    rv = y(1:3); rdot = y(4:6);
    R1 = norm(rv);
    mue = 398600.44189;
    rddot = -mue * rv / R1^3;
    dydt = [rdot; rddot];
end

% -------------------------
% GMST computation
% -------------------------
function gmst2ime = gmst2(Mjd_UT1)
    Secs = 86400; MJD_J2000 = 51544.5;
    Mjd_0 = floor(Mjd_UT1);
    UT1 = Secs * (Mjd_UT1 - Mjd_0);
    T_0 = (Mjd_0 - MJD_J2000)/36525;
    T = (Mjd_UT1 - MJD_J2000)/36525;
    gmst2 = 24110.54841 + 8640184.812866.*T_0 + 1.002737909350795.*UT1 + (0.093104-6.2e-6.*T).*T.*T;
    gmst2ime = 2*pi*Frac2(gmst2/Secs);
end

% -------------------------
% Attitude tracking helper (returns MRP error & angular velocity error)
% -------------------------
function [sigma_BR, w_BR_B] = tracking2(sigma_BN, w_BN_B, RN, w_RN_N)
    BN = ConvertAttitude(sigma_BN, 'MRP', 'DCM');
    w_RN_B = BN * w_RN_N;
    w_BR_B = w_BN_B - w_RN_B;
    BR = BN * RN';
    sigma_BR = ConvertAttitude(BR, 'DCM', 'MRP');
end

% -------------------------
% Pointing controller & attitude integrator (original point2)
% -------------------------
function [eull, xend, gvec] = point2(r0, v0, date, lat, thi, x0, t0, tf)
    dt = 1; mu = 398600; K = 0.005555555555555555; P = 0.16666666666666666;

    % time vector
    t = (t0:dt:tf)'; len = length(t);

    % initial attitude state (MRP + body rates)
    xlen = length(x0); x = zeros(xlen, len)'; x(1, :) = x0;

    % orbit state
    rm = r0; vm = v0;

    % compute ground target(s) ECEF positions for this epoch
    rmove = Rmovee2(date, lat, thi);

    % compute angle between satellite position and each candidate ground position
    angles = acos(sum(rm .* rmove, 2) ./ (sqrt(sum(rm .^ 2, 2)) .* sqrt(sum(rmove .^ 2, 2)))) * 180/pi;
    [angel, min_row] = min(angles);

    % choose pointing vector only if below threshold
    if (angel <= minAng)
        gvec = rmove(min_row, :);
    else
        gvec = [0 0 0];
    end

    % define local orbit frame: radial (er), along-track (eth), normal (eh)
    r = rm - gvec;
    er = r / norm(r);
    eh = cross(r, vm) / norm(cross(r, vm));
    eth = cross(eh, er);
    ON = [er; eth; eh];

    % reference attitude (flip sign depending on angle) - original logic retained
    if (angel <= 50)
        Ref = [1 0 0; 0 1 0; 0 0 1]' * [er; eth; eh];
    else
        Ref = [-1 0 0; 0 1 0; 0 0 -1]' * [er; eth; eh];
    end

    % reference angular velocity in reference frame
    Wref = ON' * [0, 0, sqrt(mu / (norm(rm))^3)]';

    % integrate attitude with PD control (MRP + rigid body dynamics)
    for k = 1:len-1
        if norm(x(k, 1:3)) > 1
            % MRP shadow set correction (original logic retained)
            x(k, 1:3) = -x(k, 1:3) / norm(x(k, 1:3))^2;
        end

        % compute attitude tracking error and control torque
        [sigma_BRn, w_BRn_B] = tracking2(x(k, 1:3)', x(k, 4:6)', Ref, Wref);
        u = -K * sigma_BRn - P * w_BRn_B;

        % RK4 integration step
        k1 = dt * xdot2(x(k, :), u);
        k2 = dt * xdot2(x(k, :) + k1 * 0.5 * dt, u);
        k3 = dt * xdot2(x(k, :) + k2 * 0.5 * dt, u);
        k4 = dt * xdot2(x(k, :) + k3 * 0.5 * dt, u);
        dx = dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        x(k+1, :) = x(k, :) + dx;

        if norm(x(k, 1:3)) > 1
            x(k, 1:3) = -x(k, 1:3) / norm(x(k, 1:3))^2;
        end
    end

    % return final attitude in Euler (321) degrees and final state
    xend = x(end, :);
    s = x(end, 1:3);
    eull = (ConvertAttitude(s', 'MRP', '321')') * 180/pi;
end

% -------------------------
% Rigid-body angular dynamics
% -------------------------
function dw = w_dot2(w, u)
    I = ([10,0,0; 0,5,0; 0,0,7.5]);
    WS = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    wBN_dot = (I) \ (-WS * I * w + u);
    dw = wBN_dot;
end

% -------------------------
% State derivative for attitude (MRP kinematics + angular accel)
% -------------------------
function xdot2 = xdot2(x, u)
    sigma = (x(1:3));
    omega = x(4:6);
    sigmadot = dMRP(sigma', omega');
    omegadot = w_dot2(omega', u);
    xdot2 = [sigmadot; omegadot]';
end

% -------------------------
% Plot orbit and initial satellite pose
% -------------------------
function [rm, vm, eull, xend, sat, vec] = plotorbit2(R0, V0, date0)
    % initial attitude state (MRP + angular rates)
    sigma_BN = ([0.3, -0.4, 0.5]');
    w_BN_B = ([1.0, 1.75, -2.2]') * pi/180;
    x0 = [sigma_BN', w_BN_B'];

    % compute initial pointing & attitude
    [eull, xend, gvec] = point2(R0, V0, date0, 30, 30, x0, 0, 1);

    % create satellite pose (quaternion initialized to zero Eulers)
    q = quaternion([0 0 0], "eulerd", "ZYX", "frame");
    sat = poseplot(q, R0, MeshFileName = meshF, ScaleFactor = 3, PatchFaceColor = 'b', PatchFaceAlpha = 1);

    % propagate orbit for a long time to draw trajectory (original logic)
    options = odeset('reltol', 1.e-8, 'abstol', 1.e-8);
    X0 = [R0'; V0'];
    n = sqrt(abs(mu/abs(sma^3)));
    T = abs(2*pi/n) * 100;
    [~, State] = ode45(@(t,y) dEOM2(t,y), 0:T, X0, options);
    plot3(State(:,1), State(:,2), State(:,3), '.b');

    rm = R0; vm = V0;
    vec = plot_3d_line2(gvec, rm);
end

% -------------------------
% Propagate one orbit step and update pose
% -------------------------
function [eull, xend, rm, vm, sat, vec] = oporbit2(rm, vm, i, dt, xend, date0, sat, vec)
    % propagate the orbit by dt seconds
    [rm, vm, ~] = satOrbit2(rm, vm, i, i + dt);

    % recompute pointing & attitude for the new state
    [eull, xend, gvec] = point2(rm, vm, date0, lat, long, xend, i, i + dt);

    % update pose visualization: delete old handles and create new ones
    zz = eull(1:3);
    q = quaternion([zz(1) zz(2) zz(3)], "eulerd", "ZYX", "frame");
    delete(vec);
    delete(sat);
    sat = poseplot(q, rm, MeshFileName = meshF, ScaleFactor = 2, PatchFaceColor = 'b', PatchFaceAlpha = 1);
    vec = plot_3d_line2(gvec, rm);
end

end
