function [phi_geod, h_geod, Lon,Lat,Rho] = geoc_to_geod(R_ECEF, a,b, method, tol)
% (ECEF) Geodetic coordinate conversion: (x,y,z) → (φ) ϕ,h_ϕ
% h_geod = deviation wrt to perfect sphere (0 km @ equator, 21.38 km at poles)
%
% Note: largest difference between φ and ϕ is at 45°
%
% See: Shu 2010, An iterative algorithm to compute geodetic coordinates
%      Table 1: Comparison of arithmetic operations (latitude and height)
% Escobal, Pedro Ramon: Methods of Orbit Determination. John Wiley & Sons, Inc., c.1965, p. 23.
%{
% Test case:
    phi_geoc = 50;
    rho = Earth.r1_km + 1000;
    R_ECEF = Sph2Cart(0, phi_geoc, rho);
    [phi_geod, h_geod, Lon,Lat,Rho] = geoc_to_geod(R_ECEF);
    % phi_geod = 50.1634243856166°
    % h_geod = 1012.57038273569 km
    % Rho == rho 

    % → idem to the Matlab Aerospace Toolbox function
    [phi_geod, h_geod] = geoc2geod(phi_geoc, Rho, Earth.f, Earth.r1_km)
%}
% See SOL_Geodetics.m for CelesTrak exercises

    if nargin < 5
        tol = 1e-14;  % tol must be < 1e-15 
        if nargin < 4
            method = 2;  % Exact non-iterative method (fast)
            if nargin < 2
                Earth = AstroConstants();  % Earth most of the time
                a = Earth.r1_km;
                b = Earth.r2_km;
            end
        end
    end

    e2 = 1 - (b/a)^2;   % First eccentricity squared 
    e_2 = (a/b)^2 - 1;  % Second eccentricity squared == e2 / (1 - e2)

    %--- 0) Computes Geocentric longitude λ & latitude φ and radius length ρ

    [Lon, Lat, Rho] = Cart2Sph(R_ECEF);  % Lat = phi_geoc

    %--- 1) Computes Geodetic latitude ϕ and altitude

    h_geod = [];  % (to be computed below)

    if method <= 1  &&  a ~= 6378.137
        error('>>> geoc_to_geod: Method 0 & 1 for Earth only: not applicable for normalized radiuses ==> Use method 2');
    end

    switch method

        case 0  % MATLAB Aerospace Toolbox
            [phi_geod, h_geod] = geoc2geod(Lat, Rho * 1000);
            h_geod = h_geod / 1000;

        case 1  % Long 1974, NASA TN D-7522 → FASTEST  Max err = 1.2×10⁻⁶ °

            [h_geod, phi_geod] = Geodetic_approx(Rho, Lat);  % (vectorial)

        case 2  % Heikkinen's Algorithm (1982)  Δ = 2.8e-14  *** BEST ***
            % Highly regarded closed-form (non-iterative) solution for converting Earth-Centered, Earth-Fixed (ECEF) 
            % Cartesian coordinates () into geodetic coordinates (latitude ϕ, longitude λ, and height h). 
            % It avoids the transcendental iterations found in methods like Bowring's by using a series of algebraic 
            % substitutions to solve the quartic equation involved in the transformation.
            % – Performance: No while loops or convergence checks make it ideal for high-throughput GPU or FPGA implementations.
            % – Precision: It provides sub-millimetre accuracy across the entire range of Earth's surface and near-Earth space.
            % – Stability: Unlike some iterative methods, it does not diverge near the Earth's poles or center.
        
            % r_𝛿: equatorial component along the semimajor axis (distance from Z-axis)
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  
            % r_K: vertical component parallel to the semiminor axis (distance from equatorial plane)
            z = R_ECEF(3,:);
        
            % 1. Initial calculations
            F = 54 * b^2 * z.^2;
            G = r.^2 + (1 - e2) * z.^2 - e2 * (a^2 - b^2);
            c = (e2^2 * F .* r.^2) ./ (G.^3);
            
            % 2. Intermediate variables
            s = (1 + c + sqrt(c.^2 + 2*c)).^(1/3);
            k = s + 1 + 1./s;
            P = F ./ (3 * k.^2 .* G.^2);
            Q = sqrt(1 + 2 * e2^2 * P);
            
            % 3. Calculate Surface Reference r0
            term1 = -(P .* e2 .* r) ./ (1 + Q);
            term2 = 0.5 * a^2 * (1 + 1./Q);
            term3 = (P .* (1 - e2) .* z.^2) ./ (Q .* (1 + Q));
            term4 = 0.5 * P .* r.^2;
            
            % Using max(0, ...) to prevent imaginary results from floating point noise
            r0 = term1 + sqrt(max(0, term2 - term3 - term4));
        
            % 4. Final Geodetic Coordinates
            U = sqrt((r - e2 * r0).^2 + z.^2);
            V = sqrt((r - e2 * r0).^2 + (1 - e2) * z.^2);
            
            z0 = (b^2 * z) ./ (a * V);
            
            h_geod = U .* (1 - b^2 ./ (a * V));  % h
            phi_geod = atan2d(z + e_2 * z0, r);
            % lambda = atan2d(y, x);

        case 3  % Zhu 1993: Direct exact transformation from ECEF coordinates to geodetic coordinates in closed form
            % Δ = 1.8e-13 

            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % equatorial component (distance from Z-axis)
            z = R_ECEF(3, :);

            l = e2 / 2;
            m = (r / a).^2;
            n = ((1 - e2) * z / b).^2;
            i = -(2 * l^2 + m + n) / 2;
            k = l^2 * (l^2 - m - n);
            q = (m + n - 4 * l^2).^3 / 216 + m .* n * l^2;
            D = sqrt((2 * q - m .* n * l^2) .* m .* n * l^2);
            beta = i / 3 - (q + D).^(1/3) - (q - D).^(1/3);
            t = sqrt(sqrt(beta.^2 - k) - (beta + i) / 2) - sign(m - n) .* sqrt((beta - i) / 2);
            w1 = r ./ (t + l);
            z1 = (1 - e2) * z ./ (t - l);

            phi_geod = atan2d(z1, (1 - e2) * w1);
            %Lon = 2 * atan2d(r - x, y);
            h_geod = sign(t - 1 + l) .* sqrt((r - w1).^2 + (z - z1).^2);

        case 4  % Borkowski 1987, Transformation of geocentric to geodetic coordinates without approximations
            % Δ = 6.0e-13

            if 0 % solution
                R_N = a / sqrt(1 - e2 * sind(phi_geod)^2);  
                h = h_geod;
                r = (R_N + h) * cosd(phi_geod); % == hypot(x_sat, y_sat)
                z = (R_N * (b / a)^2 + h) * sind(phi_geod);  % == z_sat
            end

            m = length(Lat);
            h_geod = zeros(1, m);
            phi_geod = zeros(1, m);
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % equatorial component (distance from Z-axis)

            for i = 1:m
                z = R_ECEF(3, i);

                % Longitude
                %lon = atan2(y, x);
                
                % Special polar case
                if r(i) < 1e-12
                    phi_geod(i) = sign(z) * 90;
                    h_geod(i) = abs(z) - b;
                else
                    % Borkowski auxiliary quantities
                    E = ((z + b) * b / a - a) ./ r(i);  % == (b * z - (a^2 - b^2)) ./ (a * r);
                    F = ((z - b) * b / a + a) ./ r(i);  % == (b * z + (a^2 - b^2)) ./ (a * r);
                    % Quartic equation: t^4 + 2*E * t^3 + 2*F * t - 1 == 0
                    
                    P = 4/3 * (E .* F + 1);
                    Q = 2 * (E.^2 - F.^2);
                    D = P.^3 + Q.^2;
                
                    % Solves cubic
                    if D >= 0
                        v = nthroot(sqrt(D) - Q, 3) - nthroot(sqrt(D) + Q, 3);
                    else
                        v = 2 * sqrt(-P) * cos(1/3 * acos(Q / (P * sqrt(-P))));
                    end
                
                    % Improve numerical stability
                    if abs(v^2) < abs(P)
                        v = -(v^3 + 2*Q) / (3*P);
                    end
                
                    G = 0.5 * (sqrt(E^2 + v) + E);
                    t = sqrt(G^2 + (F - v*G) ./ (2*G - E)) - G;
                
                    % Geodetic latitude
                    phi_geod(i) = atan2d(a * (1 - t.^2), (2 * b * t));
                
                    % Radius of curvature
                    sinLat = sind(phi_geod(i));
                    R_N = a / sqrt(1 - e2 * sinLat^2);
                
                    % Geodetic altitude
                    h_geod(i) = r(i) * cosd(phi_geod(i)) + (z + e2 * R_N .* sinLat) .* sinLat - R_N;
                end
            end


        case 10  % Bowring Iterative (convergent) solution 
            % as per CelesTrak method https://celestrak.org/columns/v02n03/

            m = length(Lat);
            h_geod = zeros(1, m);
            phi_geod = zeros(1, m);
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % r_δ equatorial component (distance from Z-axis)
            z = R_ECEF(3, :);

            for i = 1:m

                phi = Lat(i);  % from φ to ϕ
                
                if 0  % for loop
                    for iter = 1:10
                        R_N = a / sqrt(1 - e2 * sind(phi)^2);   % prime-vertical radius
                        phi_geod_ = atan2d(z(i) + e2 * R_N * sind(phi), r(i));
                        if abs(phi_geod_ - phi) < tol
                            break;
                        end
                        phi = phi_geod_;
                    end
                else  % while loop
                    phi_geod_ = Inf;
                    while abs(phi_geod_ - phi) > tol
                        phi_geod_ = phi;
                        R_N = a / sqrt(1 - e2 * sind(phi_geod_)^2);   % prime-vertical radius
                        phi = atan2d(z(i) + e2 * R_N * sind(phi_geod_), r(i));
                    end
                end

                phi_geod(i) = phi;
                h_geod(i) = r(i) / cosd(phi) - R_N;
            end

        case 11  % Fixed-point iteration of Bowring's formula to calculate the Geodetic Latitude
            % This is highly efficient, typically converging in 2-3 iterations
            % www.mathworks.com/help/aeroblks/geocentrictogeodeticlatitude.html

            m = length(Lat);
            phi_geod = zeros(1, m);
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % equatorial component (distance from Z axis)
            z = R_ECEF(3, :);

            for i = 1:m

                phi = Lat(i);  % [°] initial value
                phi_geod_ = phi;
                while 1  % typically converges in three iterations
                    beta = atand(tand(phi_geod_) * b/a); % [°]
                    phi_geod_ = atand((z(i) + b * e_2 * sind(beta)^3) / (r(i) - a * e2 * cosd(beta)^3));
                    if abs(phi_geod_ - phi) < tol
                        break;
                    end
                    phi = phi_geod_;
                end
                phi_geod(i) = phi;
            end

        case 12  % Clynch 2006, Geodetic Coordinate Conversions
            if 0
                % A) Latitude, Longitude, Height to ECEF xyz
                R_ECEF.x = (R_N + h) .* cosd(LatG_) .* cosd(Orbit.Lon);
                R_ECEF.y = (R_N + h) .* cosd(LatG_) .* sind(Orbit.Lon);
                R_ECEF.z = ((1 - e2) * R_N + h) .* sind(LatG_);

                R_ECEF.x = Orbit.ECEF.pos(1, :);
                R_ECEF.y = Orbit.ECEF.pos(2, :);
                R_ECEF.z = Orbit.ECEF.pos(3, :);

                % B) ECEF xyz to Latitude, Longitude, Height
                lon = atan2d(R_ECEF.y, R_ECEF.x);

                % The physical radius of the point and the radius in the x-y plane are used in an initial estimate of the altitude.
                r = sqrt(R_ECEF.x.^2 + R_ECEF.y.^2 + R_ECEF.z.^2);
                p = sqrt(R_ECEF.x.^2 + R_ECEF.y.^2);  % radius in the x-y plane

                % The geocentric latitude φ is computed exactly, and used as the initial value for the geodetic latitude ϕ in the loop.
                phi_geoc = atan2d(p, R_ECEF.z);
            end

            m = length(Lat);
            phi_geod = zeros(1, m);
            h_geod = zeros(1, m);
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % equatorial component
            z = R_ECEF(3, :); 

            for i = 1:m
                phi = Lat(i); % ϕ geodetic
                phi0 = phi;
                for j = 1:10  % typically converges in <= 5 iterations (for Earth geoid)
                    R_N = a / sqrt(1 - e2 * sind(phi).^2);
                    h = r(i) / cosd(phi) - R_N;  % diverges at the poles. There are two alternatives
                    phi = atand(z(i) ./ (r(i) .* (1 - e2 * R_N ./ (R_N + h))));
                    if abs(phi - phi0) < tol
                        break;
                    end
                    phi0 = phi;
                end
                phi_geod(i) = phi; % This is for positions even at earth satellite altitudes.
                h_geod(i) = h;
            end

        case 13  % Fukushima-style implementation

            m = length(Lat);
            phi_geod = zeros(1, m);
            h_geod = zeros(1, m);
            r = hypot(R_ECEF(1, :), R_ECEF(2, :));  % equatorial component
            z = R_ECEF(3, :);

            % Longitude
            %lon = atan2d(R_ECEF(2, :), R_ECEF(1, :));

            for i = 1:m
            
                % Polar singularity
                if r(i) < 1e-15
                    phi_geod(i) = sign(z(i)) * 90;
                    h_geod(i) = abs(z(i)) - b;
                    continue;
                end
                
                % Initial guess (Bowring-like)
                lat = atan2(z(i), r(i) * (1 - e2));
                
                % Fukushima/Halley iteration
                for k = 1:6  % 3 → 1e-9; 4 → 1e-11; 5 → 1e-13; 6 → 1e-14;
                
                    R_N = a / sqrt(1 - e2 * sin(lat)^2);
                    h = r(i) / cos(lat) - R_N;
                    lat_new = atan2(z(i), r(i) * (1 - e2 * R_N / (R_N + h)));
                
                    if abs(lat_new - lat) < tol
                        lat = lat_new;
                        break;
                    end
                    lat = lat_new;
                end
                
                R_N = a / sqrt(1 - e2 * sin(lat)^2);
                h_geod(i) = r(i) / cos(lat) - R_N;
                phi_geod(i) = DEG(lat);
            end

        case 14  % Nievergelt and Keeler [2000] (from Wertz 2001, p.863)
            %{
                (The publication includes references to a number of earlier, less satisfactory methods). 
                With a single iteration, this approach is good to 2x10⁻⁶° for the geodetic latitude and 1 mm in 
                geodetic altitude. Successive iterations can improve this, al­though that would rarely be needed. 
                The 4-step iterative approach is as follows:
            %}

            TODO
            

        case 20  % Series expansion

            TODO

%             m = length(Lat);
%             Lat_geod = zeros(1, m);
% 
%             % --- Constants for WGS84 ---
%             a = planet.r1_km;
%             b = planet.r2_km;
%             f = planet.f; % 1/298.257223563;  % flattening
%             e2 = 2*f - f^2;       % eccentricity squared
%             e = Earth.ecc;
%             [ e2 , e^2 ]  % ~= Earth.e
%             n = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2));  %  third flattening
% 
%             Lat(i) + DEG( e2 * (a / (a + h)) * sind(Lat(i)) * cosd(Lat(i)) )  % OK
 
% 
%             % f2 = 1/2*e^2 + 1/4*e^4 + 1/8*e^6+ 1/16*e^8 + 1/32*e^10;
%             % f2 = 2 * n - 2 * n^3 + 2 * n^5;
%             % 
%             % f4 = 1/8*e^4 + 1/8*e^6 + 3/32*e^8 + 1/16*e^10;
%             % f4 = 2*n^2 - 4*n^4;
%             % 
%             % Lat_geod = Lat(i) + DEG(f2 * sind(2 * Lat(i)) + f4 * sind(4 * Lat(i)))
%             % 
%             % atand(tand(Lat(i)) * (a/b)^2)
%             % 
%             % phi = RAD(Lat(i));
%             % Lat(i) - DEG( (1/2 * n - 13/16 * n^2) * sin(2 * phi))
% 
%             Lat(i)
%             Lat_geod(i) = geoc2geod(Lat(i), Rho(i) * 1000);  % 45.166064
%             h = Rho(i) - Earth.r1_km;
% 
%                 % h=h*1000;
%                 % a=a*1000;
% 
%             N = a / sqrt(1 - e2 * sind(Lat_geod(i))^2);  % radius of curvature in the prime vertical
%             [ tand(Lat(i)) , (N * (1 - e2) + h)/ (N + h) * tand(Lat_geod(i))]
% 
% 
%             for i = 1:m
%                 % first series approximation
%                 Lat_geod(i) = Lat(i) + DEG(a/(2*Rho(i)) * e2 * sind(2*Lat(i))); % <== OK
%                 Lat_geod(i) = Lat(i) + DEG( a*e2/(2*Rho(i)) * sind(2*Lat(i)) + ...
%                                            (a*e2/(2*Rho(i)))^2 * (sind(4*Lat(i)) - h/a * sind(2*Lat(i)))  ); % 
%             end
% 
%             a = Earth.r1_km;
%             r = Rho(i);
%             eps = e2 * (a/r);
%             theta = RAD(Lat(i));
%             phi_ref=
% 
%             %% Precompute trig
%             s = sin(theta);
%             c = cos(theta);
% 
% %% -----------------------------
% % Series expansion (correct form)
% %% -----------------------------
% 
% phi = theta ...
%     + eps * s * c ...
%     + (eps^2/2) * s * c * (1 - 2*s^2) ...
%     + (eps^3/6) * s * c * (1 - 6*s^2 + 6*s^4);
% phi_ref = geoc2geod(DEG(theta), r * 1000)
% %phi_ref = atand( tan(theta) / (1 - e2) );
% DEG(phi)-phi_ref

                % 
                % theta_deg = Lat(i);  % geocentric latitude [o]
                % %H = AltG;            % altitude [km]
                % H = Rho(min(i, length(Rho))) - a;
                % 
                % % Geocentric radius approximation
                % r = a + H;
                % a = planet.r1_km;  % [m]
                % r = Rho(min(i, length(Rho)));
                % 
                % theta = RAD(theta_deg);
                % 
                % % First-order correction coefficient
                % k1 = (a^2 * e2) / (2 * r^2);
                % 
                % % Second-order correction coefficient
                % k2 = (a^4 * e4) / (4*r^4);
                % 
                % % Series expansion
                % phi = theta + (k1 + k2)*sin(2*theta) + (a^4 * e4)/(8*r^4)*sin(4*theta);
                % 
                % Lat_geod(i) = DEG(phi);


%             % Convert to radians
%             phi_prime_rad = RAD(Lat(i));
% 
%             % --- Series Expansion Coefficients ---
%             % Expansion of (phi - phi_prime) in terms of n
%             c1 = 2*n - (2/3)*n^2 - 2*n^3;
%             c2 = (7/3)*n^2 - (8/5)*n^3;
%             c3 = (121/15)*n^3;
% 
%             % Calculate the difference using Fourier terms
%             % phi - phi_prime = c1*sin(2*phi') + c2*sin(4*phi') + c3*sin(6*phi')
%             delta_phi = c1*sin(2*phi_prime_rad) + ...
%                         c2*sin(4*phi_prime_rad) + ...
%                         c3*sin(6*phi_prime_rad);
% 
%                 % --- Corrected Series Coefficients (phi - phi_prime) ---
% % These coefficients expand the difference in terms of n
% % to O(n^3), suitable for sub-millimeter surface precision.
% c1 = 2*n + (2/3)*n^2 - 2*n^3;
% c2 = (7/3)*n^2 + (8/5)*n^3;
% c3 = (121/15)*n^3;
% 
% % Delta phi = phi - phi_prime
% delta_phi = c1*sin(2*phi_prime_rad) + ...
%             c2*sin(4*phi_prime_rad) + ...
%             c3*sin(6*phi_prime_rad);
% 
%   % --- Series Parameter ---
% % The expansion of atan(k*tan(x)) - x uses beta = (k-1)/(k+1)
% % For geocentric to geodetic: k = 1/(1-f)^2
% % This simplifies to beta = 2n / (1 + n^2)
% beta = (2*n) / (1 + n^2);
% 
% % Convert input to radians
% phi_prime_rad = deg2rad(phi_prime);
% 
% % --- Fourier Series Expansion ---
% % phi = phi' + beta*sin(2phi') + (1/2)*beta^2*sin(4phi') + (1/3)*beta^3*sin(6phi')
% delta_phi = beta * sin(2 * phi_prime_rad) + ...
%             (1/2) * (beta^2) * sin(4 * phi_prime_rad) + ...
%             (1/3) * (beta^3) * sin(6 * phi_prime_rad);
% 
% Lat(i) + DEG(delta_phi) % == 45.192423215


% %% WGS84 parameters
% a = 6378.137e3;              % semi-major axis [m]
% f = 1/298.257223563;         % flattening
% b = a*(1-f);                 % semi-minor axis [m]
% 
% e2 = 1 - (b^2/a^2);          % eccentricity squared
% 
% %% Example geocentric latitude
% theta_deg = 45;
% 
% theta = deg2rad(theta_deg);
% 
% %% ---------------------------------------------------------
% % Series expansion up to e^4
% %
% % phi = theta
% %     + (e^2/2 + e^4/4) sin(2 theta)
% %     + (e^4/8) sin(4 theta)
% %% ---------------------------------------------------------
% 
% phi_series = theta ...
%     + (e2/2 + e2^2/4)*sin(2*theta) ...
%     + (e2^2/8)*sin(4*theta);
% 
% phi_series_deg = rad2deg(phi_series);
% 
% %% Exact formula
% phi_exact = atan( tan(theta)/(1 - e2) );
% 
% phi_exact_deg = rad2deg(phi_exact);
% 
% %% Error
% error_arcsec = (phi_series_deg - phi_exact_deg)*3600;
                           
                % % Result in degrees
                % Lat_geod(i) = Lat(i) + DEG(delta_phi);
            % end


        case 990  % Solution through numerical solution of an algebraic expression
            % STUDY: doesn't converge for small φ and φ = 90°...

            m = length(Lat);
            h_geod = zeros(1, m);
            phi_geod = zeros(1, m);

            for i = 1:m
                % slope @ (x2,y2)
                % sx2 = -b * x2 / (a^2 * sqrt(1 - (x2/a)^2));
                % sx2 = -(xE - x2) / (zE - z2);
                % solve for x2: -(xE - x2) / (yE - y2) + b * x2 / (a^2 * sqrt(1 - (x2/a)^2)) = 0
                xE = hypot(R_ECEF(1, i), R_ECEF(2, i));
                zE = R_ECEF(3, i);
                x_ = 0.5 * R_ECEF(1, i);

                x2 = fzero(@(x) (zE - b * sqrt(1 - (x/a)^2)) / (x - xE) + (a^2 * sqrt(1 - (x/a)^2)) / (b * x), x_);
                % Solve[(d - b Sqrt[1 - (x/a)^2])/(x - c) == -(a^2 Sqrt[1 - (x/a)^2])/(b*x), x] → lenghty algebraic solution...
                z2 = b * sqrt(1 - (x2/a)^2);

                m = (zE - z2) / (xE - x2);
                b_ = z2 - m * x2;
                x0 = -b_ / m;
                phi_geod(i) = atan2d(zE, xE - x0);
                h_geod(i) = hypot(xE - x2, zE - z2);
            end

        case 991 % Shu 2010, An iterative algorithm to compute geodetic coordinates

            % Newton-Raphson method to solve quartic equation

            TO_BE_COMPLETED
            % p =
            % q =

            fk = (p * q)^2 - (r * q)^2 - (z * p)^2; % == 0
            fk_ = 2 * (b * p * q^2 + a * q * p^2 - a * q * r^2 - b * p * z^2);

            while delta < tol
                k_n1 = k_n - fk / fk_;
            end

            phi_geod = atan2d((a/b)^2 * z0 / r0);
            phi_geod = atan2d(a * p * z / (b * q * r));
            h_geod = k * sqrt((b * r / p)^2 + (a * z / q)^2);


        otherwise
            error('>>> geoc_to_geod: Unknown method #%d', method);

    end % switch

    R_N = a ./ sqrt(1 - e2 * sind(phi_geod).^2);  % R_N(ϕ)
    % R_N = ;  % R_N(φ)

    % Final assignment, if not computed above
    if isempty(h_geod)
        h_geod = hypot(R_ECEF(1, :), R_ECEF(2, :)) ./ cosd(phi_geod) - R_N;
    end

% STUDY:
% Shu 2020, An iterative algorithm to compute geodetic coordinates
% Zhang 2005, An alternative algebraic algorithm to convert Cartesian to geodetic coordinates.  (~)
% Hedgley 1976, An Exact Transformation from Geocentric to Geodetic Coordinates for Nonzero Altitudes.  (?)


%% Compare methods
if 0 

    Earth = AstroConstants();
    a = Earth.r1_km;
    b = Earth.r2_km;

    AltG = 1000;  % Satellite orbit
    %AltG=0;
    Rho = a + AltG;
    % no influence on φ–ϕ

    Lat = [0:0.5:90];  % small array for figure
    %Lat = [linspace(0,+90, 100000), 45];  % large array to compare evaluation speed
    Lon = linspace(0, 360, length(Lat));
    R_ECEF = Sph2Cart(Lon, Lat, Rho);

    tol = 1e-14; 
    tol = 1e-9; 
 
    Figure(0, 'Latitude ϕ geodetic study');

    phs = [];
    legends = {};
    for method = [0,1, 2:4, 10:13]
        tic
        [phi_geod, h_geod] = geoc_to_geod(R_ECEF, a, b, method, tol);
        elapsed = toc;

        [max_err, indx] = max(abs(Lat - phi_geod));
        fprintf('%2d: max φ–ϕ = %.7f° @ φ = %g°, ΔAlt = %.3e km (%.3f s)\n', method, ...
            max_err, Lat(indx), max(abs(AltG - h_geod)), elapsed);
       
        if method == 0  % uses MATLAB Aerospace Toolbox as reference
            phi_geod0 = phi_geod;
            Alt_geod0 = h_geod;
        elseif length(Lat) < 1000
            err = abs(phi_geod - phi_geod0);
            % err = abs(Alt_geod - Alt_geod0);
            phs(end+1) = plot(Lat, max(1e-15, err));
            legends{end+1} = sprintf('%d (%.2e)', method, max(err));
        end
    end

    X_label('Geocentric latitude λ [°]');
    Y_label('Geodetic latitude error Δϕ [°]');
    set(gca, 'YScale', 'log');
    plot(45 * [1,1], ylim, 'k:'); 
%plot(44.8821 * [1,1], ylim, 'r--'); % maximum error @ 44.8821°
    xlim([0,90]);
    ylim([1e-15, 2e-6]);  set(gca, 'YTick', [10.^(-15:-6)]);
    phs(end+1) = plot(xlim, tol * [1,1], 'm--');
    legends{end+1} = sprintf('tol = %g', tol);   
    plot(xlim, 1e-13 * [1,1], 'r--');
    set(gca, 'XTick', [0:15:90]);
    legend(phs, legends);
    Title(sprintf('Algorithms comparison of geodetic latitude computation: h = %g km', AltG));

    %{
        Comparative results:
         0: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.026 s)   Aerospace Toolbox for reference 
         1: max φ–ϕ = 0.1660667° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.003 s) ← Fastest NASA (BASELINE) (err < 1e-6°)
                                                                                 NON-ITERATIVE (exact):
         2: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.011 s)   Heikkinen's Algorithm (1982)
         3: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.011 s)   Zhu's Algorithm (1993)
         4: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.099 s)   Borkowski (1987)
                                                                                 ITERATIVE:
        10: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.052 s)   Bowring 
        11: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.109 s)   Bowring #2
        12: max φ–ϕ = 0.1660655° @ φ = 44.8821°, ΔAlt = 2.138e+01 km (0.057 s)   Clynch
        13: max φ–ϕ = 0.1660655° @ φ = 44.8816°, ΔAlt = 2.138e+01 km (0.048 s)   Fukushima
    %}

end % Test section

end % function