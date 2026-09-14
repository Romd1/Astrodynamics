function [R_sat, phi_geoc, R_e, R_N, R_M, x0,z0] = geod_to_pos(phi_geod, lambda, h_geod, a,b)
% Converts geodetic latitude ϕ [°] to position vector [x;y;z] and geocentric latitude, 
% at a given altitude h_ϕ above ellipsoid
% Vallado 2013, Eq.(3-7), r_δ & r_K 
% from Hedgley 1976, An Exact Transformation from Geocentric to Geodetic Coordinates for Nonzero Altitudes

% == Sph2Cart(lambda, phi_geoc, a + h_geoc) for b == a

    if nargin < 4
        % Earth = AstroConstants()
        a = 6378.137;  % = Earth.r1_km
        b = 6356.75231424518;  % = Earth.r2_km
        if nargin < 3
            h_geod = 0;
            if nargin < 2
                lambda = 0;
            end
        end
    end
    e2 = 1 - (b/a)^2;

    if lambda > 360
        error('>>> geod_to_pos: second argument θ out of bound?!')
    end

    % R_N(ϕ) Radius of curvature in Prime Vertical (terminated by minor axis)
    % or Radius of curvature in the meridian, (independent of h_geod)
    R_N = a ./ sqrt(1 + ((b/a)^2 - 1) * sind(phi_geod).^2);
    % == a ./ sqrt(1 - e2 * sind(phi_geod).^2)
    % == a^2 ./ sqrt(a^2 * cosd(phi_geod).^2 + b^2 * sind(phi_geod).^2)
    % R_N(ϕ=0)  ≡ a
    % R_N(ϕ=90) ≡ a^2/b ~ 2a - b

    % R_N(φ₀)
    % phi_geoc0 = atand(tand(phi_geod) * (b/a)^2);
    % num = a * sqrt(a^4 * sind(phi_geoc0)^2 + b^4 * cosd(phi_geoc0)^2);
    % den = b * sqrt(a^2 * sind(phi_geoc0)^2 + b^2 * cosd(phi_geoc0)^2);
    % R_N = num / den;
    % R_N = (a * sqrt(1 - (2*e2 - e2^2) * cosd(phi_geoc0)^2)) / ...
    %       (sqrt(1 - e2) * sqrt(1 - e2 * cosd(phi_geoc0)^2));
    
 % ┌─ (x,y,z) = f(ϕ) ───────────────────────────────────┐ 
    x = (R_N + h_geod) .* cosd(phi_geod) .* cosd(lambda);  % rδ @ h_ϕ  
    y = (R_N + h_geod) .* cosd(phi_geod) .* sind(lambda);
    z = ((b/a)^2 * R_N + h_geod) .* sind(phi_geod);        % rK @ h_ϕ
 % └────────────────────────────────────────────────────┘ 
    R_sat = [x; y; z];
    % Rho = norm(R_sat)

    % [~, phi_geoc] = Cart2Sph(R_sat);
    % phi_geoc = atan2d(z, hypot(x, y)) == asind(z / norm(R_sat))
 % ┌─ ϕ, h_ϕ = f(R(φ)) ───────────────────────────────────────────────────────────┐ 
    phi_geoc = atand(tand(phi_geod) .* (R_N * (b/a)^2 + h_geod) ./ (R_N + h_geod));
 % └──────────────────────────────────────────────────────────────────────────────┘ 
    % == atand(tand(phi_geod) * (b/a)^2) only when h_ϕ ≡ 0 → φ₀
    % h_geod == hypot(x, y) / cosd(phi_geod) - R_N

    z_sub = (b/a)^2 * R_N .* sind(phi_geod);

    r_delta = R_N .* cosd(phi_geod);
    r_K_ = R_N .* sind(phi_geod);  % extending down to z₀
    z0 = z_sub - r_K_;
    x0 = r_delta - z_sub .* r_delta ./ r_K_;

        r_K = R_N .* (1 - e2) .* sind(phi_geod);
        % r_sub = sqrt(r_delta^2 + r_K^2) == R_e

    % R_e(ϕ) = distance from the ellipsoid center to its surface measured along the geocentric radial direction
    R_e = a * sqrt(((1 - e2)^2 * sind(phi_geod).^2 + cosd(phi_geod).^2) ./ (1 - e2 * sind(phi_geod).^2));
    % == a * sqrt((1 - (2 * e2 - e2^2) * sind(phi_geod)^2) / (1 - e2 * sind(phi_geod)^2))
    % == R_N * sqrt(cosd(phi_geod)^2 + (b/a)^4 * sind(phi_geod).^2)
    % == sqrt((a^4 * cosd(phi_geod)^2 + b^4 * sind(phi_geod)^2) / (a^2 * cosd(phi_geod)^2 + b^2 * sind(phi_geod)^2))
    % == norm(R_sub)

    % R_e(φ) 
    % R_e = a ./ sqrt(1 + ((a/b)^2 - 1) .* sind(phi_geoc0).^2)
    %    == (a * sqrt(1 - e2)) / sqrt(1 - e2 * cosd(phi_geoc0)^2)
    % phi_geoc0 = atand(tand(phi_geod) * (b/a)^2);
    % R_e = (a * b) / sqrt(a^2 * sind(phi_geoc0)^2 + b^2 * cosd(phi_geoc0)^2);

    %--- Radius of curvature in Meridian (Meridian Radius of Curvature) 
    R_M = a * (1 - e2) ./ (1 - e2 * sind(phi_geod).^2).^(3/2);  % R_M(ϕ)
    % == R_N .* (1 - e2) ./ (1 - e2 * sind(phi_geod).^2)
    % R_M(ϕ=0)  ≡ a * (1 - e2)
    % R_M(ϕ=90) ≡ a^2/b

    if 0 
        % The following figure shows how much of a difference results from using an oblate spheroid, plotting the 
        % difference between geodetic and geocentric latitude as a function of geodetic latitude

        phi_geod = linspace(0, 90);
        lambda = 0;  % (N/A here)

        Figure(1, 'φ–ϕ difference');
        set(gca, 'Position', [0.1, 0.3, 0.85, 0.5]);
        %fontsize(gcf, scale=1.4)
        fontsize(gcf, 14, "points")

        phs = [];
        legends = {};
        for h_geod = 0:500:2000
            [R_sat, phi_geoc, R_e, R_N, R_M, x0,z0] = geod_to_pos(phi_geod, lambda, h_geod);

            if h_geod == 0
                phs(end+1) = plot(phi_geod, phi_geod - phi_geoc, '-', 'LineWidth',1.5);
            else
                phs(end+1) = plot(phi_geod, phi_geod - phi_geoc, '--', 'LineWidth',1.0);
            end
            legends{end+1} = sprintf('h_ϕ = %g km', h_geod);
        end
        phs(end+1) = plot(45 * [1,1], ylim, 'k:', 'Color',.6*[1,1,1], 'LineWidth',1.5);
        legends{end+1} = 'Maximum difference';
        xlabel('Geodetic Latitude $\varphi$ [$^\circ$]', 'Interpreter','latex', 'FontSize',16);
        ylabel('$\Delta = \phi - \varphi$ [$^\circ$]', 'Interpreter','latex', 'FontSize',16);
        sf = 1.7;
        legend(phs, legends, 'Location','best');
        Title(sprintf('Geocentric vs. Geodetic latitude'),sf);
        ylim([0, 0.2]);
    end

    if 0
        %--- Mean global radius (see AstroConstants.m)
        R1 = (2*a + b) / 3;  % Earth mean radius as defined in the 1984 World Geodetic System revision
        R2 = (a^2 * b)^(1/3);  % We).^2rtz 2001 9.1.5 Earth Oblateness
        % Radius of a sphere with the same volume as the ellipsoid
    
        %--- Local mean radius
        R_Gaussian = sqrt(R_N .* R_M);
        R_Eulerian = R_N .* R_M ./ (R_N * cosd(az)^2 + R_M * sind(az)^2);
    end

end  % function