clc;

IC.lat_deg  =  40.88;
IC.lon_deg  =  29.30;
IC.alt_ft   =  1000;
IC.psi_deg  =  0;
IC.vel_ms   =  38.58;   % 75 kts

RWY.touchdown_lat = 40.902885;
RWY.touchdown_lon = 29.321008;
RWY.stop_lat      = 40.895348;
RWY.stop_lon      = 29.300500;
RWY.heading_deg   = 244.5;
RWY.alt_ft        = 318.24 + (318.9474942 - 318.24);

FP.cruise_alt_ft    = 1200;
FP.glide_angle_deg  =    2.5;
FP.flare_alt_agl_ft =   90;
FP.flare_tau_s      =    6.0;    
FP.cruise_vel_ms    =  38.58;
FP.approach_vel_ms  =  32;
FP.landing_vel_ms   =  28;
FP.stop_vel_ms      =   0.0;
FP.leg_straight_m   = 2500;
FP.climb_rate_fps   =   10;
FP.turn3_lat_offset_deg = 0.04;

BANK_DEG   = 15;  
g_ms2      = 9.81;
R_turn_m   = FP.cruise_vel_ms^2 / (g_ms2 * tand(BANK_DEG));

ALT_RATE_FPS = 10;
DT           = 0.5;
R_earth      = 6371000;

fprintf('Cruise:   %.2f m/s  (%.1f kts)\n', FP.cruise_vel_ms,   FP.cruise_vel_ms  *1.94384);
fprintf('Approach: %.2f m/s  (%.1f kts)\n', FP.approach_vel_ms, FP.approach_vel_ms*1.94384);
fprintf('Landing:  %.2f m/s  (%.1f kts)\n', FP.landing_vel_ms,  FP.landing_vel_ms *1.94384);
fprintf('Turn radius: %.1f m\n\n', R_turn_m);

%% ============================================================
%  Rota Noktaları
% ============================================================
inbound_hdg  = RWY.heading_deg;
outbound_hdg = mod(inbound_hdg + 180, 360);

% Glideslope geometrisi
alt_drop_m   = (FP.cruise_alt_ft - RWY.alt_ft) * 0.3048;
glide_dist_m = alt_drop_m / tand(FP.glide_angle_deg);

FAF_lat = RWY.touchdown_lat + cosd(outbound_hdg) * glide_dist_m / R_earth * (180/pi);
FAF_lon = RWY.touchdown_lon + sind(outbound_hdg) * glide_dist_m / ...
          (R_earth * cos(RWY.touchdown_lat*pi/180)) * (180/pi);

bearing_to_FAF = atan2d((FAF_lon - RWY.touchdown_lon)*cos(RWY.touchdown_lat*pi/180), ...
                         FAF_lat - RWY.touchdown_lat);

fprintf('TD->FAF bearing: %.2f°  (beklenen: %.1f°)\n', mod(bearing_to_FAF,360), outbound_hdg);
fprintf('Glide distance:  %.1f m\n\n', glide_dist_m);

lon0 = IC.lon_deg;
climb_dist_m = (FP.cruise_alt_ft - IC.alt_ft) / FP.climb_rate_fps * IC.vel_ms;
desc_dist_m  = climb_dist_m;

lat1 = IC.lat_deg + FP.leg_straight_m / R_earth * (180/pi);
lat2 = lat1 + climb_dist_m / R_earth * (180/pi);
lat3 = lat2 + FP.leg_straight_m / R_earth * (180/pi);
lat4 = lat3 + desc_dist_m  / R_earth * (180/pi);
lat5 = lat4 + FP.leg_straight_m / R_earth * (180/pi);

%% ============================================================
%  ARC geometri
% ============================================================
T3_dir = +1;
T3_ctr = geo_move([FAF_lat, FAF_lon], R_turn_m, mod(inbound_hdg + 90*T3_dir, 360), R_earth);
T3_p0  = geo_move(T3_ctr, R_turn_m, mod(180 - 90*T3_dir, 360), R_earth);
T3_p1  = [FAF_lat, FAF_lon];

T1_dir = +1;
T1_p0  = [lat5, lon0];
T1_ctr = geo_move(T1_p0, R_turn_m, mod(0  + 90*T1_dir, 360), R_earth);
T1_p1  = geo_move(T1_ctr, R_turn_m, mod(90 - 90*T1_dir, 360), R_earth);

L2_p0 = T1_p1;
R_turn_lon = R_turn_m / (R_earth * cos(T3_p0(1)*pi/180) * pi/180);
L2_p1 = [T1_p1(1), T3_p0(2) - R_turn_lon];

T2_dir = +1;
T2_p0  = L2_p1;
T2_ctr = geo_move(T2_p0, R_turn_m, mod(90  + 90*T2_dir, 360), R_earth);
T2_p1  = geo_move(T2_ctr, R_turn_m, mod(180 - 90*T2_dir, 360), R_earth);

L3_p0 = T2_p1;
EARLY_M = 20;
L3_p1   = geo_move(T3_p0, EARLY_M, 0, R_earth);
T3_p0   = L3_p1;

% L2 bacak uzunluğu kontrolü
L2_len_m = (L2_p1(2) - L2_p0(2)) * R_earth * cos(L2_p0(1)*pi/180) * pi/180;


%% ============================================================
%  ARC flare
% ============================================================
gamma_rad    = FP.glide_angle_deg * pi/180;
h_flare_m    = FP.flare_alt_agl_ft * 0.3048;

R_flare_m    = h_flare_m / (1 - cos(gamma_rad));
x_flare_m    = R_flare_m * sin(gamma_rad);

FLARE_lat = RWY.touchdown_lat + cosd(outbound_hdg) * x_flare_m / R_earth * (180/pi);
FLARE_lon = RWY.touchdown_lon + sind(outbound_hdg) * x_flare_m / ...
            (R_earth * cos(RWY.touchdown_lat*pi/180)) * (180/pi);
FLARE_alt_ft = RWY.alt_ft + FP.flare_alt_agl_ft;

FAF_alt_ft  = IC.alt_ft;
gs_dist_m   = (FAF_alt_ft - FLARE_alt_ft) * 0.3048 / tand(FP.glide_angle_deg);

d_flare_td = sqrt(((FLARE_lat-RWY.touchdown_lat)*pi/180*R_earth)^2 + ...
             ((FLARE_lon-RWY.touchdown_lon)*pi/180*R_earth*cos(FLARE_lat*pi/180))^2);

fprintf('gamma           : %.2f deg\n',  FP.glide_angle_deg);
fprintf('h_flare         : %.1f ft  (%.2f m)\n', FP.flare_alt_agl_ft, h_flare_m);
fprintf('R_flare         : %.1f m\n',   R_flare_m);
fprintf('x_flare (yatay) : %.1f m\n',   x_flare_m);
fprintf('FAF→FLARE dist  : %.1f m\n',   gs_dist_m);
fprintf('FLARE→TD dist   : %.1f m  (beklenen: %.1f m)\n\n', d_flare_td, x_flare_m);

%% ============================================================
%  Segment Listesi
% ============================================================
segs = {};

segs{end+1} = mkline([IC.lat_deg,IC.lon_deg],[lat1,lon0], IC.alt_ft,IC.alt_ft, IC.vel_ms,IC.vel_ms);
segs{end+1} = mkline([lat1,lon0],[lat2,lon0], IC.alt_ft,FP.cruise_alt_ft, IC.vel_ms,FP.cruise_vel_ms);
segs{end+1} = mkline([lat2,lon0],[lat3,lon0], FP.cruise_alt_ft,FP.cruise_alt_ft, FP.cruise_vel_ms,FP.cruise_vel_ms);
segs{end+1} = mkline([lat3,lon0],[lat4,lon0], FP.cruise_alt_ft,IC.alt_ft, FP.cruise_vel_ms,IC.vel_ms);
segs{end+1} = mkline([lat4,lon0],[lat5,lon0], IC.alt_ft,IC.alt_ft, IC.vel_ms,IC.vel_ms);
segs{end+1} = mkarc(T1_p0,T1_p1,T1_ctr,  0, 90, T1_dir, IC.alt_ft, IC.vel_ms);
segs{end+1} = mkline(L2_p0,L2_p1, IC.alt_ft,IC.alt_ft, IC.vel_ms,IC.vel_ms);
segs{end+1} = mkarc(T2_p0,T2_p1,T2_ctr, 90,180, T2_dir, IC.alt_ft, IC.vel_ms);
segs{end+1} = mkline(L3_p0,L3_p1, IC.alt_ft,IC.alt_ft, IC.vel_ms,FP.approach_vel_ms);
segs{end+1} = mkarc(T3_p0,T3_p1,T3_ctr,180,240, T3_dir, IC.alt_ft, FP.approach_vel_ms);

sg.type='GLIDE';
sg.p0=[FAF_lat,FAF_lon]; sg.p1=[FLARE_lat,FLARE_lon];
sg.h0=FAF_alt_ft; sg.h1=FLARE_alt_ft;
sg.v0=FP.approach_vel_ms; sg.v1=FP.approach_vel_ms;
sg.center=[0,0]; sg.hdg0=inbound_hdg; sg.hdg1=inbound_hdg; sg.dir=1;
sg.R_flare=R_flare_m; sg.gamma_rad=gamma_rad;
segs{end+1} = sg;

sf.type   = 'FLARE_ARC';
sf.p0     = [FLARE_lat, FLARE_lon];
sf.p1     = [RWY.touchdown_lat, RWY.touchdown_lon];
sf.h0     = FLARE_alt_ft;
sf.h1     = RWY.alt_ft;
sf.v0     = FP.approach_vel_ms;
sf.v1     = FP.landing_vel_ms;
sf.hdg0   = inbound_hdg;
sf.hdg1   = inbound_hdg;
sf.dir    = 1;
sf.center = [0,0];
sf.R_flare    = R_flare_m;
sf.gamma_rad  = gamma_rad;
segs{end+1} = sf;

% STOP
ss.type='LINE';
ss.p0=[RWY.touchdown_lat,RWY.touchdown_lon]; ss.p1=[RWY.stop_lat,RWY.stop_lon];
ss.h0=RWY.alt_ft; ss.h1=RWY.alt_ft;
ss.v0=FP.landing_vel_ms; ss.v1=0;
ss.center=[0,0]; ss.hdg0=inbound_hdg; ss.hdg1=inbound_hdg; ss.dir=1;
ss.R_flare=0; ss.gamma_rad=0;
segs{end+1} = ss;

%% ============================================================
%  Zaman vektörleri
% ============================================================
time_vec=[]; lat_vec=[]; lon_vec=[];
alt_vec=[];  psi_vec=[]; vel_vec=[];
t_now = 0;

for si = 1:numel(segs)
    seg = segs{si};

    if strcmp(seg.type,'ARC')
        dhdg    = mod((seg.hdg1 - seg.hdg0)*seg.dir + 360, 360);
        arc_len = R_turn_m * dhdg * pi/180;
        seg_dur = arc_len / seg.v0;

    elseif strcmp(seg.type,'FLARE_ARC')
        arc_len = seg.R_flare * seg.gamma_rad;
        avg_v   = (seg.v0 + seg.v1) / 2;
        seg_dur = arc_len / avg_v;
        fprintf('FLARE_ARC: R=%.1f m, arc_len=%.1f m, avg_v=%.2f m/s, dur=%.1f s\n',...
                seg.R_flare, arc_len, avg_v, seg_dur);

    else
        dlat_m  = (seg.p1(1)-seg.p0(1))*(pi/180)*R_earth;
        dlon_m  = (seg.p1(2)-seg.p0(2))*(pi/180)*R_earth*cos(seg.p0(1)*pi/180);
        seg_len = sqrt(dlat_m^2+dlon_m^2);
        avg_v   = max((seg.v0+seg.v1)/2, 0.5);
        seg_dur = seg_len / avg_v;
    end
    if seg_dur < DT; seg_dur = DT; end

    t_seg = (0:DT:seg_dur)';
    n     = numel(t_seg);
    tau   = t_seg / seg_dur;

    seg_lat = zeros(n,1); seg_lon = zeros(n,1);
    seg_alt = zeros(n,1); seg_psi = zeros(n,1); seg_vel = zeros(n,1);

    if strcmp(seg.type,'ARC')
        ctr      = seg.center;
        ang0_raw = seg.hdg0 - 90*seg.dir;
        dhdg     = mod((seg.hdg1 - seg.hdg0)*seg.dir + 360, 360);
        for k = 1:n
            ang_k      = ang0_raw + seg.dir * dhdg * tau(k);
            seg_lat(k) = ctr(1) + R_turn_m*cosd(ang_k)/R_earth*(180/pi);
            seg_lon(k) = ctr(2) + R_turn_m*sind(ang_k)/(R_earth*cos(ctr(1)*pi/180))*(180/pi);
        end
        seg_psi = seg.hdg0 + seg.dir * dhdg * tau;
        seg_alt = seg.h0 * ones(n,1);
        seg_vel = seg.v0 * ones(n,1);

    elseif strcmp(seg.type,'GLIDE')
        for k = 1:n
            seg_lat(k) = seg.p0(1) + tau(k)*(seg.p1(1)-seg.p0(1));
            seg_lon(k) = seg.p0(2) + tau(k)*(seg.p1(2)-seg.p0(2));
        end
        seg_alt = seg.h0 + tau*(seg.h1-seg.h0);
        seg_psi = inbound_hdg * ones(n,1);
        seg_vel = seg.v0 * ones(n,1);

    elseif strcmp(seg.type,'FLARE_ARC')
        Rf       = seg.R_flare;
        gam      = seg.gamma_rad;

        theta    = gam * (1 - tau);
        h_agl    = Rf * (1 - cos(theta));
        x_horiz  = Rf * (sin(gam) - sin(theta));

        for k = 1:n
            seg_lat(k) = seg.p0(1) + cosd(inbound_hdg)*x_horiz(k)/R_earth*(180/pi);
            seg_lon(k) = seg.p0(2) + sind(inbound_hdg)*x_horiz(k) / ...
                         (R_earth*cos(seg.p0(1)*pi/180))*(180/pi);
        end

        seg_alt = RWY.alt_ft + h_agl * 3.28084;
        seg_psi = inbound_hdg * ones(n,1);
        seg_vel = seg.v0 + tau*(seg.v1-seg.v0);

        fprintf('  FLARE_ARC : h_giriş=%.1f ft, h_çıkış=%.1f ft\n',...
                seg_alt(1), seg_alt(end));
        fprintf('  x_giriş=0 m, x_çıkış=%.1f m  (beklenen: %.1f m)\n',...
                x_horiz(end)*R_earth/(R_earth), x_flare_m);

    else  % LINE / STOP
        for k = 1:n
            seg_lat(k) = seg.p0(1) + tau(k)*(seg.p1(1)-seg.p0(1));
            seg_lon(k) = seg.p0(2) + tau(k)*(seg.p1(2)-seg.p0(2));
        end
        seg_alt = seg.h0 + tau*(seg.h1-seg.h0);
        seg_vel = seg.v0 + tau*(seg.v1-seg.v0);
        dlat_m = (seg.p1(1)-seg.p0(1))*(pi/180)*R_earth;
        dlon_m = (seg.p1(2)-seg.p0(2))*(pi/180)*R_earth*cos(seg.p0(1)*pi/180);
        if sqrt(dlat_m^2+dlon_m^2) > 1
            hdg = mod(atan2d(dlon_m,dlat_m), 360);
        else
            hdg = seg.hdg0;
        end
        seg_psi = hdg * ones(n,1);
    end

    if isempty(time_vec)
        time_vec = [time_vec; t_now+t_seg];
        lat_vec  = [lat_vec;  seg_lat]; lon_vec = [lon_vec; seg_lon];
        alt_vec  = [alt_vec;  seg_alt]; psi_vec = [psi_vec; seg_psi];
        vel_vec  = [vel_vec;  seg_vel];
    else
        time_vec = [time_vec; t_now+t_seg(2:end)];
        lat_vec  = [lat_vec;  seg_lat(2:end)]; lon_vec = [lon_vec; seg_lon(2:end)];
        alt_vec  = [alt_vec;  seg_alt(2:end)]; psi_vec = [psi_vec; seg_psi(2:end)];
        vel_vec  = [vel_vec;  seg_vel(2:end)];
    end
    t_now = t_now + seg_dur;

    fprintf('Seg %2d %-10s  dur=%6.1f s  psi=%.0f→%.0f  h=%.0f→%.0f\n',...
        si, seg.type, seg_dur, seg_psi(1), seg_psi(end), seg_alt(1), seg_alt(end));
end

%% ============================================================
%  Psi
% ============================================================
psi_norm = psi_vec;
for k = 2:numel(psi_norm)
    d = psi_norm(k) - psi_norm(k-1);
    if d >  180; psi_norm(k:end) = psi_norm(k:end) - 360; end
    if d < -180; psi_norm(k:end) = psi_norm(k:end) + 360; end
end
fprintf('\npsi_norm aralık: %.1f .. %.1f deg\n', min(psi_norm), max(psi_norm));

%% ============================================================
%  ALT
% ============================================================
BLEND_SIGMA_S = 6.0;
alt_rl = alt_vec;
for k = 2:numel(alt_rl)
    dtk       = time_vec(k) - time_vec(k-1);
    max_delta = ALT_RATE_FPS * dtk;
    delta     = alt_vec(k) - alt_rl(k-1);
    alt_rl(k) = alt_rl(k-1) + max(min(delta,max_delta), -max_delta);
end
sigma_n    = BLEND_SIGMA_S / DT;
half_k     = ceil(3*sigma_n);
kk         = (-half_k:half_k)';
gauss_k    = exp(-0.5*(kk/sigma_n).^2);
gauss_k    = gauss_k / sum(gauss_k);
alt_padded = [repmat(alt_rl(1),half_k,1); alt_rl; repmat(alt_rl(end),half_k,1)];
alt_smooth = conv(alt_padded, gauss_k, 'valid');
alt_final  = alt_smooth;

%% ============================================================
%  save
% ============================================================
TRJ.time       = time_vec;
TRJ.lat_deg    = lat_vec;
TRJ.lon_deg    = lon_vec;
TRJ.alt_ft     = alt_final;
TRJ.psi_deg    = psi_norm;
TRJ.vel_ms     = vel_vec;
TRJ.FP         = FP;
TRJ.RWY        = RWY;
TRJ.IC         = IC;
TRJ.total_time = time_vec(end);
TRJ.R_flare_m  = R_flare_m;

ApproachParams = [
    FAF_lon
    FAF_lat
    FAF_alt_ft

    FLARE_lon
    FLARE_lat
    FLARE_alt_ft

    RWY.touchdown_lon
    RWY.touchdown_lat
    RWY.alt_ft

    RWY.heading_deg
    FP.glide_angle_deg
    FP.flare_alt_agl_ft
];

TRJ.ApproachParams = ApproachParams;

save('trajectory.mat','TRJ','ApproachParams');
fprintf('\nToplam süre: %.1f s  |  Nokta sayısı: %d\n', TRJ.total_time, numel(TRJ.time));

%% ============================================================
%  plot
% ============================================================
figure('Color','w','Name','Ground Track');
plot(lon_vec, lat_vec,'b-','LineWidth',1.5, 'DisplayName', 'Trajectory'); hold on;
plot(RWY.touchdown_lon,RWY.touchdown_lat,'rv','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Touchdown');
plot(IC.lon_deg,IC.lat_deg,'gs','MarkerSize',12,'MarkerFaceColor','g','DisplayName','IC');
plot(FAF_lon,FAF_lat,'c^','MarkerSize',10,'MarkerFaceColor','c','DisplayName','FAF');
plot(FLARE_lon,FLARE_lat,'y^','MarkerSize',10,'MarkerFaceColor','y','DisplayName','FLARE start');
xlabel('Lon'); ylabel('Lat'); title('Ground Track'); legend show; grid on; axis equal; hold off;

figure('Color','w','Name','Reference Signals');
subplot(3,1,1);
plot(time_vec,alt_final,'b','LineWidth',1.5); hold on;
plot(time_vec,alt_vec,'b--','LineWidth',0.8); hold off;
ylabel('Alt (ft)'); grid on; title('Altitude');

subplot(3,1,2);
plot(time_vec,psi_norm,'r','LineWidth',1.5);
ylabel('Psi (deg)'); grid on; title('Heading');

subplot(3,1,3);
plot(time_vec,vel_vec,'g','LineWidth',1.5);
yline(FP.cruise_vel_ms,'g--','Cruise');
yline(FP.approach_vel_ms,'b--','Approach');
yline(FP.landing_vel_ms,'r--','Landing');
ylabel('Vel (m/s)'); xlabel('t (s)'); grid on; title('Velocity');

% Flare yakın çekim
figure('Color','w','Name','Arc Flare Profile');
flare_mask = alt_vec >= RWY.alt_ft & alt_vec <= (RWY.alt_ft + FP.flare_alt_agl_ft + 10) & ...
             time_vec > time_vec(round(end*0.7));
if any(flare_mask)
    t_fl  = time_vec(flare_mask) - time_vec(find(flare_mask,1));
    h_fl  = alt_vec(flare_mask) - RWY.alt_ft;
    v_fl  = vel_vec(flare_mask);
    subplot(2,1,1);
    plot(t_fl, h_fl*0.3048,'b','LineWidth',2);
    xlabel('t (s)'); ylabel('h AGL (m)'); title('Flare Altitude Profile (arc)'); grid on;
    subplot(2,1,2);
    plot(t_fl, v_fl,'r','LineWidth',2);
    xlabel('t (s)'); ylabel('V (m/s)'); title('Flare Speed Profile'); grid on;
end

%% ============================================================
% yardımcı fonksiyonlar
% ============================================================
function s = mkline(p0,p1,h0,h1,v0,v1)
    s.type='LINE'; s.p0=p0; s.p1=p1;
    s.h0=h0; s.h1=h1; s.v0=v0; s.v1=v1;
    s.center=[0,0]; s.hdg0=0; s.hdg1=0; s.dir=1;
    s.R_flare=0; s.gamma_rad=0;
end

function s = mkarc(p0,p1,ctr,hdg0,hdg1,dir,h,v)
    s.type='ARC'; s.p0=p0; s.p1=p1; s.center=ctr;
    s.hdg0=hdg0; s.hdg1=hdg1; s.dir=dir;
    s.h0=h; s.h1=h; s.v0=v; s.v1=v;
    s.R_flare=0; s.gamma_rad=0;
end

function p2 = geo_move(p1, dist_m, bearing_deg, Re)
    p2(1) = p1(1) + cosd(bearing_deg) * dist_m / Re * (180/pi);
    p2(2) = p1(2) + sind(bearing_deg) * dist_m / (Re*cos(p1(1)*pi/180)) * (180/pi);
end