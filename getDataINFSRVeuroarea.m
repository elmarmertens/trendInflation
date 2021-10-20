%% construct europrix data set for project "paddington"
% Note: all simple rates of change are converted into log-differences

%% clear workspace
clear variables
clear global
close all
fclose all;
clc


%% process HICP
HICPimport     = importdata('HICPcoresdw.csv'); % headline and core

Ylabel = {'HICP', 'HICPcore'};
Ny     = length(Ylabel);
datalabel   = 'HICPeuroarea';

spfsdwdates = datenum(HICPimport.textdata(6:end,1), 'yyyymmm');

hicpdata = HICPimport.data;


if length(spfsdwdates) ~= length(hicpdata)
    error('import dimensions do not match')
end

% flip order (SDW starts with youngest data)
spfsdwdates = flipud(spfsdwdates);
hicpdata = flipud(hicpdata);

%% Monthly data
hicpInflation = diff(log(hicpdata)) * 1200;
dates         = spfsdwdates(2:end);

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
spfsdwdates = datenum(GDPDimport.textdata(6:end,1), 'yyyyqq');
gdpddata = GDPDimport.data;
if length(spfsdwdates) ~= length(gdpddata)
    error('import dimensions do not match')
end
% flip order (SDW starts with youngest data)
spfsdwdates = flipud(spfsdwdates);
gdpddata = flipud(gdpddata);

gdpd     = diff(log(gdpddata)) * 400;
spfsdwdates = spfsdwdates(2:end);
% plot
figure
plot(spfsdwdates, gdpd)
xtickdates(spfsdwdates)
title('GDPD')


ud = union(dates, spfsdwdates);

dates2 = genrMdates(1995,2021,1);
dates2 = dates2(dates2 >= ud(1));
dates2 = dates2(dates2 <= ud(end));



T      = length(dates2);
Ydata2 = NaN(T, Ny);

ndxIn  = ismember(dates2, dates);
ndxOut = ismember(dates, dates2);
Ydata2(ndxIn,1:2) = hicpInflation(ndxOut,:);

ndxIn  = ismember(dates2, spfsdwdates);
ndxOut = ismember(spfsdwdates, dates2);
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

%% Add SPF
close all


Ylabel      = {'HICP', 'HICPcore', 'GDPD', 'SPF12m', 'SPF24m', 'SPFltr'};
Ny          = length(Ylabel);
datalabel   = 'INFSRVeuroarea';

dates3 = dates2;
T      = length(dates3);
Ydata3 = NaN(T, Ny);

Ydata3(:,1:3) = Ydata2;

SPFimport   = importdata('SPFeuroprix.csv');

spfdates    = datenum(SPFimport.textdata(6:end,1), 'yyyy-mm-dd');
spfdata     = SPFimport.data(:,2:end); % first column is HICP

Nspf       = size(spfdata,2);

% flip time
% spfdates  = flipud(spfdates);
% spfdata   = flipud(spfdata);

% set spfdates to first of month
% [y,m,~]  = datevec(spfdates);
% spfdates = datenum(y,m,1);

% simple rates of change into log-differences
spfdata = log(1 + spfdata / 100) * 100;

% loop over spfdates, and place each obs as needed
% NOTE: SPFfile has uneven spacing of dates!

for t = 1 : length(spfdates)
    
    % 12m
    i = 1;
    [y,m] = datevec(spfdates(t));
    thisdate = datenum(y-1,m,1); % always dated at first of month
    ndx = dates3 == thisdate;
    Ydata3(ndx, 3+i) = spfdata(t,i);
    
    % 24m
    i = 2;
    [y,m] = datevec(spfdates(t));
    thisdate = datenum(y-2,m,1); % always dated at first of month
    ndx = dates3 == thisdate;
    Ydata3(ndx, 3+i) = spfdata(t,i);
    
    % LTR
    i = 3;
    if ~isnan(spfdata(t,i))
        [y,m] = datevec(spfdates(t));
        if m > 3
            thisdate = datenum(y,m-3,1); % always dated at first of month
        elseif m == 3
            thisdate = datenum(y-1,12,1);
        else
            error('m lower than three not expected')
        end
        ndx = dates3 == thisdate;
        Ydata3(ndx, 3+i) = spfdata(t,i);
    end
    
end

figure
hold on
plot(dates3, Ydata3, '-x')
xtickdates(dates2)

dates          = dates3;
Ydata          = Ydata3;
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
