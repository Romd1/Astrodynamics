function [phi_geoc, h_geod] = geod_to_geoc(phi_geod, a,b, h0)
% Retrieves geocentric latitude φ and geodetic altitude h_ϕ (above a) for the geodetic latitude ϕ
% h_geod ≠ R_e
% [phi_gc, h_gd] = geod_to_geoc(60.4350632585099, Earth.r1_km, 0.8 * Earth.r1_km) == [50, 831.969398327804]
% [phi_gc, h_gd] = geod_to_geoc([0, 90], a,b) == [0, a-b]
%{
% Test case:
    phi_geod = 50;
    h0 = 1000;
    [phi_geoc, h_geod] = geod_to_geoc(phi_geod, Earth.r1_km,Earth.r2_km, h0)
    % phi_geoc = 49.8364099440923°
    % h_geod = 1012.51026593283 km

    % → identical the Matlab Aerospace Toolbox function → using h_φ
    [phi_geoc, radii] = geod2geoc(phi_geod, h_geod, Earth.f,Earth.r1_km)
    % phi_geoc = 49.8364099440923°
    % radii - Earth.r1_km == h0

    [phi_geoc, radii] = geod2geoc(phi_geod, h0, Earth.f,Earth.r1_km)
    % phi_geoc = 49.8361320923381°
    % radii - Earth.r1_km = 987.489785146104 < h0 OK
%}
    if nargin < 4
        h0 = 0;
    end

    rho = a + h0;

    if 0  % Derivation #1 (for h0 = 0)
        %{
            Using identity for simplifying:
           ┌────────────────────────────┐ 
           │ 1 - e² sin(ϕ)² ≡ a² / R_N² │
           └────────────────────────────┘ 
        
            % Based on equation along the circle encompassing the ellipsoid
            x = (R_N + h_geod) * cosd(phi_geod)  % == x_sat 
            z = (R_N * (1 - e2) + h_geod) * sind(phi_geod)  % == z_sat
            r = hypot(x, z)  % == a
            [ r, norm(R_sat) , a ]
            
            h = h_geod  % (sought solution)
            a^2 == (R_N + h)^2 * cosd(phi_geod)^2 + (R_N * (1 - e2) + h)^2 * sind(phi_geod)^2
    
            % → isolate h
        %}
        e2 = 1 - (b/a)^2;
        R_N = a ./ sqrt(1 - e2 * sind(phi_geod).^2);

        % Quadratic form solution
        B = rho^2 ./ R_N;
        C = -rho^2 + R_N.^2 .* (2 * rho^2 ./ R_N.^2 + e2^2 * sind(phi_geod).^2 - 1);
        % (maybe this can still be simplified?)
        h_geod = sqrt(B.^2 - C) - B;

        % Once we know h_ϕ:
        [R_sat, phi_geoc] = geod_to_pos(phi_geod, 0, h_geod, a,b);
    
    else  % Derivation #2 (OK for any h0)

        % a) Computes sub-satellite point on the ellipsoid surface
        [R_sub, ~, ~, R_N] = geod_to_pos(phi_geod, 0, 0, a,b);  % (h = 0)
        x_sub = R_sub(1, :);
        z_sub = R_sub(3, :);

        % b) Computes normal line perpendicular to ellipsoid at sub-satellite point
        if 0
            m = b * -x_sub ./ (a^2 * sqrt(1 - (x_sub / a).^2));  % tangential slope to sub-sat point
            m = m .* sign(phi_geod);  % handles negative ϕ
            m = -1 ./ m;  % perpendicular
        else
            m = cotd(90 - phi_geod); % (direct solution)
        end
        b_ = z_sub - m .* x_sub;  % ordinate at origin ≠ ellipsoid parameter b
 
        % c) Intersecting normal vector at R_sub with spherical Earth @ ρ
        % m * xx + b_normal == sqrt(a_^2 - xx^2) @ x = x_sat
        % [ m^2*xx^2 + 2*b_normal*m*xx + b_normal^2 , a^2 - xx^2 ]
        % [ (m^2+1)*xx^2 + 2*b_normal*m*xx + b_normal^2 - a^2, 0 ]
        if 0
            A = m.^2 + 1;
            B = 2 * b_ .* m;
            C = b_.^2 - rho.^2;
            x_sat = (-B + sqrt(B.^2 - 4 * A .* C)) ./ (2*A);  % Quadratic equation solution
        else
            x_sat = (-b_ .* m + sqrt(rho.^2 * (m.^2 + 1) - b_.^2)) ./ (m.^2 + 1);  % other form
        end

        z_sat = m .* x_sat + b_;

        h_geod = x_sat ./ cosd(phi_geod) - R_N;  % h_ϕ (Vermeille's algoithm)
        phi_geoc = atan2d(z_sat, x_sat);  % φ

        % Handles NaNs
        k = ismember(phi_geod, [-90, 0, 90]);
        phi_geoc(k) = phi_geod(k);
        h_geod(phi_geod == 0) = h0;
        h_geod(abs(phi_geod) == 90) = a + h0 - b;
    end
    % could also be solved in R,ϕ :
    % y = R sin(ϕ) + z₀ ; y = sqrt((a+h0)² - R² cos(ϕ)²)
    % R^2 c1 + c2 == sqrt(c3 - R^2 c4)
    % (R^2 c1 + c2)^2 == c3 - R^2 c4
    % ... quadratic solve for R

    % ==> ! verify ± solution depending on quadrants...

end  % function