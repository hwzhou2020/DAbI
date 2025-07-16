function [loc_p_subpixel] = subPixelFit(central_line, loc_p, fitRange, sub_res)
%   subpixel_fit refines a peak location to subpixel accuracy using spline fitting.
%
%   Inputs:
%       central_line - 1D signal (e.g., intensity profile along the central row)
%       loc_p        - Initial peak location (integer index)
%       fitRange     - Number of pixels to include in the fitting window
%       sub_res      - Subpixel resolution (e.g., 0.1 means 10× finer sampling)
%
%   Output:
%       loc_p_subpixel - Refined peak location with subpixel precision

    delta = floor(fitRange/2);
    left = loc_p - delta;
    right = loc_p + delta;

    if left < 1 || right > length(central_line)
        loc_p_subpixel = loc_p;
        return;
    end

    x = 1:fitRange;
    fitRange_line = central_line(left:right);
    fitted_line = fit(x', double(gather(fitRange_line)), 'smoothingspline');
    x_dense = 1:sub_res:fitRange;
    f_line = fitted_line(x_dense);

    [~, locs] = findpeaks(f_line, 'SortStr', 'descend');

    if isempty(locs)
        loc_p_subpixel = loc_p;
    else
        loc_p_subpixel = (locs(1) - 1) * sub_res + left;
    end
end