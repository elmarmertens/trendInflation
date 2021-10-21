%% construct europrix data set for project "paddington"
% Note: all simple rates of change are converted into log-differences

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% NOTE:
% make sure to prune SDWINFTRM from empty "-" fields in the SUPERCORE column (pre 2002)

%% import SDW file
SDWimport = importdata('SDWINFTRM.csv'); % headline and core
sdwdates  = datenum(SDWimport.textdata(6:end,1), 'yyyymmm');
sdwdates  = flipud(sdwdates);

% colheaders = SDWimport.textdata(2,:)

[y, m]      = datevec(sdwdates);
dates      = datenum(y,m,1); % FRED convention: dated at beginning of month
dates      = dates(2:end); % since data will be differenced
T          = length(dates);

%% collect HICP
Ylabel      = {'HICP', 'HICPcore'};
Ny          = length(Ylabel);     
datalabel   = 'HICP';

ndx = [12 13];    
hicpdata = SDWimport.data(:,ndx);
hicpdata = flipud(hicpdata);

hicpInflation = diff(log(hicpdata)) * 1200;

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


%% Add GDPD

GDPDimport  = importdata('SDWGDPD.csv');
sdwdates    = datenum(GDPDimport.textdata(6:end,1), 'yyyyqq');
gdpddata    = GDPDimport.data;
% flip order (SDW starts with youngest data)
sdwdates = flipud(sdwdates);
gdpddata = flipud(gdpddata);

gdpd     = diff(log(gdpddata)) * 400;
sdwdates = sdwdates(2:end);
% plot
figure
plot(sdwdates, gdpd)
xtickdates(sdwdates)

GDPD   = NaN(T,1);
[~, ndxIn, ndxOut] = intersect(dates, sdwdates);
GDPD(ndxIn,:)      = gdpd(ndxOut);


datalabel = 'SDWINF';
Ydata  = [Ydata, GDPD];
Ylabel = cat(2, Ylabel, 'GDPD');

Ny             = size(Ydata, 2);
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

%% add SUPERCORE
ndx         = 1;    
spxcoredata = SDWimport.data(:,ndx);
spxcoredata = flipud(spxcoredata);
SPXCORE     = diff(log(spxcoredata)) * 1200;

% use 12m avg to handle seasonality
SPXCORE12   = sumK(SPXCORE, 12) / 12;

figure
plot(dates, [SPXCORE SPXCORE12])
xtickdates(dates)

datalabel = 'SDWINFCORE';
Ydata     = [Ydata, SPXCORE12];

Ylabel    = cat(2, Ylabel, 'SPXCORE');

Ny             = size(Ydata, 2);
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

%% add TRIMMED
ndx      = 2 : 11;    
trimdata = SDWimport.data(:,ndx);
trimdata = flipud(trimdata);
TRM      = trimdata(2:end,:); % since other series were differenced
% transform into log rates
TRM      = log(1 + TRM / 100) * 100;

figure
plot(dates, TRM)
xtickdates(dates)

datalabel = 'SDWINFTRM';
Ydata     = [Ydata, TRM];

TRMlabel = cell(size(ndx));
trmcut   = 5 : 5 : 50;
for n = 1 : length(ndx)
    TRMlabel{n} = sprintf('TRM-%d', trmcut(n));
end
Ylabel    = cat(2, Ylabel, TRMlabel);

Ny             = size(Ydata, 2);
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

