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

%% IMPORT ALL INPUTFILES
SDWimport   = importdata('SDWINFTRM.csv'); % headline and core
GDPDimport  = importdata('SDWGDPD.csv');
SPFimport   = importdata('SDWINFSRV.csv');


%% import SDW file
sdwdates  = datenum(SDWimport.textdata(6:end,1), 'yyyymmm');
sdwdates  = flipud(sdwdates);

% colheaders = SDWimport.textdata(2,:)

[y, m]     = datevec(sdwdates);
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


% reset NaN (do not forget!)
Ydata(yNaNndx) = NaN;

%% Add GDPD

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

% save INF for later
Yinf     = Ydata;
INFlabel = Ylabel;

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


% reset NaN (do not forget!)
% Ydata(yNaNndx) = NaN;


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
Ydata     = [Yinf, SPXCORE12]; % using 12m to avoid seasonals

Ylabel    = cat(2, INFlabel, 'SPXCORE');

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


% reset NaN (do not forget!)
% Ydata(yNaNndx) = NaN;

%% add TRIMMED (but ex SUPERCORE)
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

TRMlabel = cell(size(ndx));
trmcut   = 5 : 5 : 50;
for n = 1 : length(ndx)
    TRMlabel{n} = sprintf('TRM-%d', trmcut(n));
end

ndx      = 3; % TRM15
TRM      = TRM(:,ndx);
TRMlabel = TRMlabel(ndx);

Ydata     = [Yinf, TRM];
Ylabel    = cat(2, INFlabel, TRMlabel);
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


datalabel = 'SDWINFTRMCORE';

Ydata     = [Yinf, TRM, SPXCORE12];
Ylabel    = cat(2, INFlabel, TRMlabel, 'SPXCORE');
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
% reset NaN (do not forget!)
% Ydata(yNaNndx) = NaN;

%% add SPF
spfdates     = datenum(SPFimport.textdata(6:end,1), 'yyyy-mm-dd');
spfdates     = flipud(spfdates);

% grab only longer-term SPF data
% col 1: HICP
% Col 2-4: 12m,24m,LTR
Nspf     = 3;
spfdata  = SPFimport.data(:,2:4);
spfdata  = flipud(spfdata);
spfdata  = log(1 + spfdata / 100) * 100;

Yspf     = NaN(T,Nspf); % target array 

for t = 1 : length(spfdates)
    
    % 12m: SDW stores survey response at target data
    i = 1;
    [y,m]        = datevec(spfdates(t));
    thisdate     = datenum(y-1,m,1); % always dated at first of month
    ndx          = dates == thisdate;
    Yspf(ndx, i) = spfdata(t,i);
    
    % 24m: SDW stores survey response at target data
    i = 2;
    [y,m]        = datevec(spfdates(t));
    thisdate     = datenum(y-2,m,1); % always dated at first of month
    ndx          = dates == thisdate;
    Yspf(ndx, i) = spfdata(t,i);
    
    % LTR: SDW stores survey response at last day of quarter SPF was conducted
    i = 3;
    % the 
    if ~isnan(spfdata(t,i))
        [y,m] = datevec(spfdates(t));
        % Figure out last month prior to survey quarter
        if m > 3
            thisdate = datenum(y,m-3,1); 
        elseif m == 3
            thisdate = datenum(y-1,12,1);
        else
            error('m lower than three not expected') % LTR data enterd for months 3,6,9 and 12
        end
        ndx              = dates == thisdate;
        Yspf(ndx, i) = spfdata(t,i);
    end
    
end

SPFlabel = {'SPF-12m', 'SPF-24m', 'SPF-LTR'};
figure
hold on
plot(dates, Yspf, '-x')
xtickdates(dates)
title('SPF')

datalabel = 'SDWINFTRMSRV';
Ydata     = [Yinf, TRM, Yspf];
Ylabel    = cat(2, INFlabel, TRMlabel, SPFlabel);


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


% reset NaN (do not forget!)
% Ydata(yNaNndx) = NaN;

datalabel = 'SDWINFSRV';
Ydata     = [Yinf, Yspf];
Ylabel    = cat(2, INFlabel, SPFlabel);
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