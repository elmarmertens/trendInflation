function h = plotynoncompact(dates,y, thisColor, linewidth)
% PLOTYNONCOMPACT ...
%   h = plotynoncompact(dates,y, thisColor)
%   ...

if nargin < 3
    thisColor = [0 0 0];
end
if nargin < 4
    linewidth = 1;
end

hanni = [];
if iscompact(y)
    hanni = plot(dates, y, '-', 'linewidth', linewidth, 'color', thisColor);
else
    theseDates = dates;
    while ~isempty(y)
        firstT = find(~isnan(y), 1, 'first');
        if ~isempty(firstT)
            lastT  = firstT;
            while (lastT < length(y)) && ~isnan(y(lastT+1))
                lastT = lastT + 1;
            end
            if lastT > firstT
                plot(theseDates(firstT:lastT), y(firstT:lastT), '-', 'color', thisColor)
            else
                hanni = plot(theseDates(firstT), y(firstT), 'x','color', thisColor);
                plot(theseDates(firstT), y(firstT), 'o','color', thisColor)
            end
            y      = y(lastT+1:end);
            theseDates = theseDates(lastT+1:end);
        else
            y = [];
        end
    end
end

if nargout > 0
    h = hanni;
end
