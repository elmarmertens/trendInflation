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


%% input files

COREimport     = importdata('HICPcoresdw.csv'); % headline and core

Ylabel = {'HICP', 'HICPcore'};
Ny     = length(Ylabel);     
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
% figure
% plot(dates, hicpInflation)
% xtickdates(dates)


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

%% Add GDPD


Ylabel      = {'HICP', 'HICPcore', 'GDPD'};
Ny          = length(Ylabel);    
datalabel   = 'INFeuroarea';


GDPDimport  = importdata('GDPDsdw.csv');
sdwdates = datenum(GDPDimport.textdata(6:end,1), 'yyyyqq');
gdpddata = GDPDimport.data;
if length(sdwdates) ~= length(gdpddata)
    error('import dimensions do not match')
end
% flip order (SDW starts with youngest data)
sdwdates = flipud(sdwdates);
gdpddata = flipud(gdpddata);

gdpd     = diff(log(gdpddata)) * 400;
sdwdates = sdwdates(2:end);
% plot
figure
plot(sdwdates, gdpd)
xtickdates(sdwdates)



ud = union(dates, sdwdates);

dates2 = genrMdates(1995,2020,1);
dates2 = dates2(dates2 >= ud(1));
dates2 = dates2(dates2 <= ud(end));



T      = length(dates2);
Ydata2 = NaN(T, Ny);

ndxIn  = ismember(dates2, dates);
ndxOut = ismember(dates, dates2);
Ydata2(ndxIn,1:2) = hicpInflation(ndxOut,:);

ndxIn  = ismember(dates2, sdwdates);
ndxOut = ismember(sdwdates, dates2);
Ydata2(ndxIn,3) = gdpd(ndxOut);


% plot
figure
hold on
plot(dates2, Ydata2, '-x')
xtickdates(dates2)

dates          = dates2;
Ydata          = Ydata2;
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

