function fh = newfigure(s)
% NEWFIGURE
% creates a new figure window 
% implements also the bugfix for plotting patches under unix

%   Coded by  Elmar Mertens, em@elmarmertens.com

figure
clf reset

if nargout > 0
    fh = gcf;
end

% prior to rh5, we needed the following bugfix 
% ... to avoid crashes after plotting patches with large faces
% if isFRB && isunix
%    set(gcf, 'renderer', 'zbuffer')
% end

% set(gca, 'fontsize', getpref('embox', 'fontsize'))

h = rotate3d;
set(h, 'rotatestyle', 'box');
% set(gcf, 'Renderer', 'painters')

if nargin > 0
    set(gcf, 'name', s)
end