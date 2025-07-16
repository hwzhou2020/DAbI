function phi_aberration = add_aberration_zernike(x, y, zernike_coeffs, R)
% Generate a phase aberration map from Zernike coefficients
% ---------------------------------------------------------
% x, y            : meshgrid coordinates (same size)
% zernike_coeffs  : vector of coefficients ordered by Noll index (j = 1, 2, ...)
% R               : radius of the pupil (same unit as x and y)

    % Convert to polar coordinates
    rho   = sqrt(x.^2 + y.^2) / R;
    theta = atan2(y, x);

    % % Create a circular pupil mask
    % pupil_mask = rho <= 1;
    % rho(~pupil_mask)   = 0;
    % theta(~pupil_mask) = 0;

    % Initialize phase map
    phi_aberration = zeros(size(x));

    % Accumulate each Zernike term
    for j = 1:length(zernike_coeffs)
        if zernike_coeffs(j) == 0, continue; end
        [n, m] = noll_to_nm(j);  % Convert Noll index to (n, m)
        Z_nm = zernike_polynomial(n, m, rho, theta); % .* pupil_mask;
        phi_aberration = phi_aberration + zernike_coeffs(j) * Z_nm;
    end

    % Wrap phase to range [-pi, pi]
    phi_aberration = mod(phi_aberration + pi, 2*pi) - pi;
end
function [n, m] = noll_to_nm(j)
% Convert Noll index j to Zernike order (n, m)
% This follows Noll’s 1976 definition, assuming j >= 1

    % Find n such that j is within its group
    n = 0;
    while j > (n + 1)*(n + 2)/2
        n = n + 1;
    end

    % Zero-based index within group of same n
    p = j - n*(n + 1)/2 - 1;

    % For n = 0
    if n == 0
        m = 0;
        return;
    end

    % Even n: m = 0, ±2, ±4,...
    % Odd  n: m = ±1, ±3, ±5,...
    if mod(n, 2) == 0
        if p == 0
            m = 0;
        else
            k = ceil(p / 2);
            mag = 2 * k;
            m = mag * (-1)^(p + 1); % alternate sign
        end
    else
        mag = 2 * floor(p / 2) + 1;
        m = mag * (-1)^(p);
    end
end
function Z = zernike_polynomial(n, m, rho, theta)
% Compute the normalized Zernike polynomial Z_n^m(rho, theta)
% with orthonormality over the unit disk

    m_abs = abs(m);
    Rnm = zeros(size(rho));

    for s = 0:((n - m_abs)/2)
        c = (-1)^s * factorial(n - s) / ...
            (factorial(s) * factorial((n + m_abs)/2 - s) * factorial((n - m_abs)/2 - s));
        Rnm = Rnm + c * rho.^(n - 2*s);
    end

    % Compute normalization factor
    if m == 0
        norm_factor = sqrt(n + 1);
    else
        norm_factor = sqrt(2 * (n + 1));
    end

    % Apply angular term
    if m > 0
        Z = norm_factor * Rnm .* cos(m_abs * theta);
    elseif m < 0
        Z = norm_factor * Rnm .* sin(m_abs * theta);
    else
        Z = norm_factor * Rnm;
    end
end