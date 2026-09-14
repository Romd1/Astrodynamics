clear

Earth = AstroConstants();
a = Earth.r1_km;
b = Earth.r2_km;
e2 = 1 - (b / a)^2;  % == Earth.ecc^2
ellipsoid = wgs84Ellipsoid('kilometer'); % for distance() function

phi_geod = sort([0:15:90]);
%phi_geod = sort([0:5:90]);
lambda_ECEF = 0;

h0 = 0;
%h0 = 1000;
Rho = a + h0;

phi_geoc0 = 0 * phi_geod;
segment_dist = nan(size(phi_geod));
segment_dist1 = segment_dist;
segment_dist2 = segment_dist;

fprintf('h₀ = %g km\n', h0);
fprintf('  ϕ [⁰]  │   φ [⁰]  │  Δ = φ−ϕ  │ Δ [arcmin] │ Δ [arcsec] │    Δh⟂  │ N–S dist │  Δϕ [km]\n');
fprintf('─────────┼──────────┼───────────┼────────────┼────────────┼─────────┼──────────┼──────────\n');
for i = 1:length(phi_geod)

    [phi_geoc, h_geod] = geod_to_geoc(phi_geod(i), a, b, h0);
    phi_geoc0(i) = atand(tand(phi_geod(i)) * (b/a)^2); 

    [~, R_geoc] = geoc_to_pos(phi_geoc0(i),lambda_ECEF, 0, a,b);  
    h_geoc = Rho - R_geoc;
    
    delta_lat = phi_geod(i) - phi_geoc;  % Δlat = ϕ – φ [°]
    delta_h = h_geoc - h_geod;  % h_ϕ – h_φ [km] ==> cm-order

    % here computed radiuses R_e, R_N and R_M are independent on h (R_sub would be) 
    [~, ~, ~, R_N, R_M] = geod_to_pos(phi_geod(i),lambda_ECEF, 0, a,b); 
    
    % R_M = Meridional radius of curvature M(ϕ)  
    NS_dist_km = R_M * RAD(delta_lat);  % (arc length)
    if 0  % validation
        md = meridianarc(RAD(phi_geoc), RAD(phi_geod(i)), ellipsoid);  % [m]
        % Calculate geodesic distance (method "gc" is default for ellipsoid)
        [arclen, az] = distance(phi_geoc, 0, phi_geod(i), 0, ellipsoid);  % @ λ = 0°
        [ md/1000, arclen/1000, NS_dist_km]
    end

    if i > 1
        % (R_M is computed independent of h)
        [~, ~, ~, ~, R_M_mid] = geod_to_pos((phi_geod(i) + phi_geod(i-1))/2,lambda_ECEF, 0, a,b);   
        segment_dist(i) = R_M_mid * RAD(phi_geod(i) - phi_geod(i-1));  % (arc length)
        % Analytical solution
        segment_dist1(i) = a * (ellipticE(RAD(phi_geoc0(i)), e2) - ellipticE(RAD(phi_geoc0(i - 1)), e2));
        % Matlab function 
        segment_dist2(i) = distance(phi_geod(i), 0, phi_geod(i-1), 0, ellipsoid);
    end

    delta_lon = delta_lat;
    EW_dist_km = R_N * cosd(phi_geod(i)) * RAD(delta_lon); % ds = (R_N cos(ϕ)) dλ
    if 0  % validation
        [arclen, az] = distance(phi_geod(i), 0, phi_geod(i), delta_lon, ellipsoid);
        [ arclen/1000, EW_dist_km]
    end
    
    fprintf('%7.5g⁰ │ %7.4f⁰ │ %8.6f⁰ │ %8.3f′  │ %8.3f″  │ %5.2f m │ %5.2f km │ %.1f km\n', ...
        phi_geod(i), phi_geoc, delta_lat * [1, 60, 60^2], delta_h * 1000, NS_dist_km, segment_dist(i));
end
fprintf('─────────┴──────────┴───────────┴────────────┴────────────┴─────────┴──────────┴──────────\n');
fprintf('                                                         1/4 circumference: ∑ = %.1f km\n', ...
    sum(segment_dist(2:end)));

