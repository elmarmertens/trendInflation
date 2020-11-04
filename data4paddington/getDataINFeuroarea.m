%% construct europrix data set for project "cambridge"
% Note: all simple rates of change are converted into log-differences

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% prepare variables

% quarters are time-stamped on the first day of the quarter (FRED convention)
% note: the following could also be automated, but might be good to set a
% few things manually (also forces us to check things when updating the data)


% % [ UPDATE ]
% dates    = genrMdates(1997,2020);
% dates    = dates(2:end-2); % from 1997M2 until 2020M10
% 
% % other
% T        = length(dates);
% Ny       = 1;     % actual inflation plus 3 SPF
% Ydata    = NaN(T,Ny); % this will be the output of this program


%% input files

COREimport     = importdata('HICPcoresdw.csv'); % headline and core
% GDPDimport     = importdata('GDPDsdw.csv');

Ylabel = {'HICP', 'HICPcore'};
Ny     = length(Ylabel);     % actual inflation plus 3 SPF
datalabel   = 'HICPeuroarea';


%% process HICP
sdwdates = datenum(COREimport.textdata(6:end,1), 'yyyymmm');

hicpdata = COREimport.data;
    

if length(sdwdates) ~= length(hicpdata)
    error('import dimensions do not match')
end

% flip order (SDW starts with youngest data)
sdwdates = flipud(sdwdates);
hicpdata = flipud(hicpdata);

%% Monthly data
hicpInflation = diff(log(hicpdata)) * 1200;
dates         = sdwdates(2:end);

% plot
figure
plot(dates, hicpInflation)
xtickdates(dates)


Ydata          = hicpInflation;
yNaNndx        = isnan(Ydata);
Ydata(yNaNndx) = 0;

mat2fortran(sprintf('%s.yData.txt', datalabel), Ydata)
logical2fortran(sprintf('%s.yNaN.txt', datalabel), yNaNndx)
mat2fortran(sprintf('%s.dates.txt', datalabel), dates)

filename = sprintf('%s.settings.txt', datalabel);
fid = fopen(filename, 'wt');
fprintf(fid, 'Ny = %d\n', size(Ydata,2));
fprintf(fid, 'T  = %d\n', size(Ydata,1));
fprintf(fid, 'YLABEL:\n');
for n = 1 : Ny
    fprintf(fid, '%s\n', Ylabel{n});
end
fclose(fid);
display(filename);
type(filename)
hrulefill

return

% 
% %% convert into quarterly data 
% if doQuarterly
%     m       = month(sdwdates);
%     qndx    = ismember(m, 3 : 3 : 12);
%     qlevels = hicpdata(qndx);
%     
%     hicpInflationQ = [NaN; diff(log(qlevels))] * 400;
%     % % check
%     % check = sumK(hicpInflation, 3) / 3;
%     % checkdiff(hicpInflationQ, check(qndx));
%     
%     qdates  = quarterlydates(sdwdates(qndx)); % beginning of quarter as in FRED
%     
% 
%     dates         = qdates(2:end);
%     hicpInflation = hicpInflationQ(2:end);
% else
%     hicpInflation = [NaN; diff(log(hicpdata))] * 1200;
%     dates         = sdwdates(2:end);
% end
% 
% %% plot
% figure
% plot(dates, hicpInflation)
% datetick('x')
% 
% 
% %% finish: store data
% Ydata = hicpInflation;
% 
% if doQuarterly
%     datalabel = strcat(datalabel, 'quarterly');
% end
% 
% % treat missing data
% yNaNndx        = isnan(Ydata);
% Ydata(yNaNndx) = 0;
% 
% % store data
% mat2fortran(sprintf('%s2020.yData.txt', datalabel), Ydata)
% logical2fortran(sprintf('%s2020.yNaN.txt', datalabel), yNaNndx)
% mat2fortran(sprintf('%s2020.dates.txt', datalabel), dates)
% 
% 
% %% finish
% dir *.txt
% 
% dockAllFigures
