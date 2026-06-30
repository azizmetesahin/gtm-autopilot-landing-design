%% ============================================================
% trajectory.mat içindeki trajectory'den
% WPs = [lon_deg, lat_deg, alt_ft, TAS_kts]
% waypoint matrisi oluştur
%
% Ayrıca altitude reference map oluşturur:
% AltMap = [station_m, alt_ref_ft]
% ============================================================

load('trajectory.mat');   % TRJ yapısını yükler

%% ------------------------------------------------------------
% PARAMETRELER
% -------------------------------------------------------------
WP_DIST_M      = 500;     % waypoint aralığı [m]
ALT_MAP_DIST_M = 20;      % altitude map aralığı [m] 

R_earth = 6371000;

BANK_DEG = 20;
g_ms2    = 9.81;

%% ------------------------------------------------------------
% Trajectory verileri 
% -------------------------------------------------------------
lat = TRJ.lat_deg(:);
lon = TRJ.lon_deg(:);
alt = TRJ.alt_ft(:);
vel = TRJ.vel_ms(:);

N = length(lat);

%% ------------------------------------------------------------
% Yol mesafesi - kümülatif
% -------------------------------------------------------------
ds = zeros(N,1);

for k = 2:N

    dlat = deg2rad(lat(k) - lat(k-1));
    dlon = deg2rad(lon(k) - lon(k-1));

    lat_avg = deg2rad(0.5*(lat(k) + lat(k-1)));

    dx = R_earth * dlon * cos(lat_avg);
    dy = R_earth * dlat;

    ds(k) = sqrt(dx^2 + dy^2);

end

s = cumsum(ds);

[s_unique, ia] = unique(s, 'stable');

lat_u = lat(ia);
lon_u = lon(ia);
alt_u = alt(ia);
vel_u = vel(ia);

%% ------------------------------------------------------------
% WP ÖRNEKLEME
% -------------------------------------------------------------
s_wp = (0:WP_DIST_M:s_unique(end))';

if s_wp(end) < s_unique(end)
    s_wp(end+1) = s_unique(end);
end

%% ------------------------------------------------------------
% WAYPOINT INTERPOLATION
% -------------------------------------------------------------
lat_wp = interp1(s_unique, lat_u, s_wp, 'linear');
lon_wp = interp1(s_unique, lon_u, s_wp, 'linear');
alt_wp = interp1(s_unique, alt_u, s_wp, 'linear');
vel_wp = interp1(s_unique, vel_u, s_wp, 'linear');

tas_kts = vel_wp * 1.94384;

%% ------------------------------------------------------------
% ALTITUDE REFERENCE MAP
%
% FORMAT:
% AltMap = [station_m, alt_ref_ft]
%
% station_m:
%   rotanın başlangıcından itibaren yatay kümülatif mesafe [m]
%
% alt_ref_ft:
%   o mesafedeki istenen altitude referansı [ft]
% -------------------------------------------------------------

s_alt = (0:ALT_MAP_DIST_M:s_unique(end))';

if s_alt(end) < s_unique(end)
    s_alt(end+1) = s_unique(end);
end

% Burada pchip kullanıyoruz çünkü altitude profilini daha yumuşak verir.
% Tracker içinde daha sonra lineer interpolation yapılacak.
alt_ref = interp1(s_unique, alt_u, s_alt, 'pchip');

AltMap = [s_alt, alt_ref];

%% ------------------------------------------------------------
% Psi Hesabı
% -------------------------------------------------------------
psi_deg = zeros(size(lat_wp));

for k = 1:length(lat_wp)-1

    lat1 = deg2rad(lat_wp(k));
    lat2 = deg2rad(lat_wp(k+1));

    dlon = deg2rad(lon_wp(k+1) - lon_wp(k));

    x = sin(dlon) * cos(lat2);

    y = cos(lat1)*sin(lat2) - ...
        sin(lat1)*cos(lat2)*cos(dlon);

    psi = atan2d(x,y);

    psi_deg(k) = mod(psi,360);

end

% son noktaya bir önceki heading'i ver
psi_deg(end) = psi_deg(end-1);

%% ------------------------------------------------------------
% WPs MATRİSİ
% [lon_deg, lat_deg, alt_ft, TAS_kts]
% -------------------------------------------------------------
WPs = [ ...
    lon_wp, ...
    lat_wp, ...
    alt_wp, ...
    tas_kts];

ApproachParams = TRJ.ApproachParams;

%% ------------------------------------------------------------
% save
% -------------------------------------------------------------
save('WPs.mat','WPs','psi_deg','AltMap','ApproachParams');

fprintf('Toplam waypoint sayısı: %d\n', size(WPs,1));
fprintf('Toplam altitude map noktası: %d\n', size(AltMap,1));
fprintf('Toplam rota uzunluğu: %.2f m\n', s_unique(end));

%% ------------------------------------------------------------
% PLOT
% -------------------------------------------------------------
figure('Color','w');

subplot(3,1,1)
plot(lon_wp,lat_wp,'b.-');
axis equal;
grid on;
xlabel('Lon');
ylabel('Lat');
title('Waypoint Ground Track');

subplot(3,1,2)
plot(psi_deg,'r','LineWidth',1.5);
grid on;
ylabel('\psi (deg)');
xlabel('WP Index');
title('Waypoint Heading');

subplot(3,1,3)
plot(s_unique, alt_u, 'k-', 'LineWidth', 1.0); hold on;
plot(AltMap(:,1), AltMap(:,2), 'r--', 'LineWidth', 1.5);
plot(s_wp, alt_wp, 'bo', 'MarkerSize', 4);
grid on;
xlabel('Station along trajectory [m]');
ylabel('Altitude [ft]');
legend('Original TRJ altitude','AltMap','WP altitude');
title('Altitude Reference Map');
