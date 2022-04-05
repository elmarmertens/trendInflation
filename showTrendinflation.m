%% plot ucsv fortran results
path(pathdef)
clear var
close all
clc
wrap = pwd;

addpath matbox


% datalabel = 'INF';
% T = 741; % needs to be adapted to the length of the actual input data

datalabel = 'SDWINFTRM'; 
T = 302; % needs to be adapted to the length of the actual input data

mcmcModel = 'gapSVeqf'; % 'eqf' or 'cta'
% mcmcModel = 'gapCONST'; 

samplestamp = sprintf('T%d', T);
 
datadir    = pwd;
timestamp  = 'notrendslopes';



%% get data

y  = importdata(fullfile(datadir, sprintf('%s.yData.txt', datalabel)));
Ny = size(y,2);
Nstates = Ny * 2;
Nsv = Ny;


foo = importdata(fullfile(datadir, sprintf('%s.settings.txt', datalabel)));
Ylabel = foo(4:end,1);

% Ylabel = cellfun(@(x) sprintf('%s-%s', datalabel, x), Ylabel, 'UniformOutput', false) 
Ylabel = strcat(datalabel, '-', Ylabel);

filename = fullfile(datadir, sprintf('%s.yNaN.txt', datalabel));
if exist(filename, 'file')
    yNaN    = logical(importdata(filename));
    y(yNaN) = NaN;
end

ybar = mean(y,2, 'omitnan');

filename = sprintf('%s.label.txt', datalabel);


dates = importdata(fullfile(datadir, sprintf('%s.dates.txt', datalabel)));

if ~isempty(T)
    dates = dates(1:T);
    y     = y(1:T,:);
    ybar  = ybar(1:T);
    yNaN  = yNaN(1:T,:);
else
    T = length(y);
end

%% get parameters

if isempty(samplestamp)
    fileext = sprintf('%s.%s.dat', datalabel, mcmcModel);
else
    fileext = sprintf('%s.%s.%s.dat', datalabel, samplestamp, mcmcModel);
end
if ~isempty(timestamp)
    fileext = sprintf('%s.%s', timestamp,fileext);
end


type(fullfile(datadir, strcat('settings.', fileext)));


tau         = importdata(fullfile(datadir, sprintf('TAU1.%s', fileext)));

sigtrend    = importdata(fullfile(datadir, sprintf('SIGTREND.%s', fileext)));
siggap      = NaN(size(sigtrend,1), size(sigtrend,2), Ny);
for n = 1 : Ny
    siggap(:,:,n)  = importdata(fullfile(datadir, sprintf('SIGGAP%d.%s', n, fileext)));
end
% maxlambda   = importdata(fullfile(fortrandir, sprintf('MAXLAMBDA.%s', fileext)));
% hvarbar     = loaddat(fullfile(fortrandir, sprintf('HVARBAR.%s', fileext)));
% F           = loaddat(fullfile(fortrandir, sprintf('F.%s', fileext)));


Ndraws      = length(tau);

TAU = NaN(T, 12, Ny);
for i = 1 : Ny
    TAU(:,:,i)    = importdata(fullfile(datadir, sprintf('TAU%d.%s', i, fileext)));
end

STATES = NaN(T, 12, Nstates);
for i = 1 : Nstates
    STATES(:,:,i)    = importdata(fullfile(datadir, sprintf('STATES%d.%s', i, fileext)));
end

% % slopes
% shockslopes = loaddat(fullfile(fortrandir, sprintf('SHOCKSLOPES.%s', fileext)));
% trendslopes = loaddat(fullfile(fortrandir, sprintf('TRENDSLOPES.%s', fileext)));
% 
% % hinno
% hSigmavech = loaddat(fullfile(fortrandir, sprintf('HSIGMA.DRAWS.%s', fileext)));
% 
% % hgap0
% hgap0 = loaddat(fullfile(fortrandir, sprintf('HGAP0.%s', fileext)));

datalabel   = strrep(datalabel, '_', '');
trendname   = strcat('Trend',datalabel);
SV1name     = strcat('SVTREND',datalabel);
SV2name     = strcat('SVGAP',datalabel);


ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 8]; % 90 percent bands



%% plot
newfigure
plotCIlines(tau(:,ndxmean), tau(:,ndxtails), dates, [], [], true)
hold on
xtickdates(dates)
grid on
wrapcf(trendname, wrap)


ndx = dates >= datenum(2000,1,1);
newfigure
plotCIlines(tau(ndx,ndxmean), tau(ndx,ndxtails), dates(ndx), [], [], true)
hold on
theseLims = ylim;
if theseLims(1) > 0 
    ylim([0 theseLims(2)])
end
if theseLims(2) < 4 
    ylim([theseLims(1) 4])
end
xtickdates(dates(ndx))
grid on
wrapcf(strcat('RECENT', trendname), wrap)

%% Trend, data and shades
for i = 1 : Ny 
    newfigure 
    plotCIlines(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
    hold on
    plotynoncompact(dates, y(:,i), [1 0 0]);
    xtickdates(dates)
    grid on
    title(sprintf('%s Trend level for %s: %6.2f%% (90%% band: %4.2f%% - %4.2f%%)', Ylabel{i}, datestr(dates(end), 'mmm yyyy'), TAU(end,ndxmean,i), TAU(end,ndxtails,i)))
    wrapcf(sprintf('tau%s', Ylabel{i}), wrap)
end

for i = 1 : Ny
    newfigure 
    plotCIlines(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
    hold on
    
    plotynoncompact(dates, y(:,i), [1 0 0]);
    xtickdates(dates(dates >= datenum(1999,1,1)))
    
    grid on
    title(sprintf('%s Trend level for %s: %6.2f%% (90%% band: %4.2f%% - %4.2f%%)', Ylabel{i}, datestr(dates(end), 'mmm yyyy'), TAU(end,ndxmean,i), TAU(end,ndxtails,i)))
    wrapcf(sprintf('RECENTtau%s', Ylabel{i}), wrap)
end

%% Trend, 12m-data and shades
for i = 1 : Ny 
    newfigure 
    plotCIlines(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
    hold on
    y12m = sumK(y(:,i), 12) / 12;
    plotynoncompact(dates, y12m, [1 0 0]);
    xtickdates(dates)
    grid on
    title(sprintf('%s Trend level for %s: %6.2f%% (90%% band: %4.2f%% - %4.2f%%)', Ylabel{i}, datestr(dates(end), 'mmm yyyy'), TAU(end,ndxmean,i), TAU(end,ndxtails,i)))
    wrapcf(sprintf('tau12m%s', Ylabel{i}), wrap)
end

%% recent 12m data
for i = 1 : Ny 
    newfigure 
    plotCIlines(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
    hold on
    if iscompact(y(:,i))
        y12m = sumK(y(:,i), 12) / 12;
        plot(dates, y12m, 'r-', 'LineWidth', 2);
    else
        plot(dates, y(:,i), 'r-o', 'linewidth', 2);
    end

    xtickdates(dates(dates >= datenum(1999,1,1)))
    ylim([min([ylim, 0]) max([ylim, 4])])
    grid on
    title(sprintf('%s Trend level for %s: %6.2f%% (90%% band: %4.2f%% - %4.2f%%)', Ylabel{i}, datestr(dates(end), 'mmm yyyy'), TAU(end,ndxmean,i), TAU(end,ndxtails,i)))
    wrapcf(sprintf('RECENTtau12m%s', Ylabel{i}), wrap)
end

% for i = 1 : Ny
%     newfigure 
%     plotCIlines(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
%     hold on
%     
%     plotynoncompact(dates, y(:,i), [1 0 0]);
%     xtickdates(dates(dates >= datenum(1999,1,1)))
%     
%     grid on
%     wrapcf(sprintf('RECENTtau%s', Ylabel{i}), wrap)
% end

%% Trend SV
newfigure
plotCIlines(sigtrend(:,ndxmean), sigtrend(:,ndxtails), dates)
hold on
xtickdates(dates)
datetick('x', 'keeplimits')
wrapcf(SV1name, wrap)

%% SVgap and data

for i = 1 : Ny
    
    firstObs = find(~yNaN(:,i),1,'first');
    newfigure
    subplot(2,1,1)
    plotCIlines(STATES(:,ndxmean,Ny + i), STATES(:,ndxtails,Ny + i), dates)
    hold on
    shades(dates, 1 : T < firstObs)
    xtickdates(dates)
    grid on
    
    subplot(2,1,2)
    plotCIlines(siggap(:,ndxmean,i), siggap(:,ndxtails,i), dates)
    hold on
    xtickdates(dates)
    shades(dates, 1 : T < firstObs)
    datetick('x', 'keeplimits')
    wrapcf(sprintf('GAP%s', Ylabel{i}), wrap)
end

%% report little table
coreNdx = 1;
hrulefill
tableNdx = find(ismember(dates, [datenum(2007,12,1);datenum(2015,12,1);datenum(2016,12,1)';datenum(2017,3:3:12,1)';datenum(2018,3:3:12,1)']));
tableNdx = cat(1, tableNdx, (T-12:T)');
tableNdx = unique(tableNdx);
for i = 1 : length(tableNdx)
    
    t = tableNdx(i);

    
    fprintf('%s ', datestr(dates(t), 'YYYY-mmm'))
    fprintf('\t %6.4f', TAU(t,ndxmean,coreNdx))
    fprintf('\t (%6.4f -- %6.4f)', TAU(t,ndxtails,coreNdx))
    
    fprintf('\n')
end

%% store data table
% writedatatable(".", sprintf('Trend-%s', Ylabel{coreNdx}), ...
%     dates, TAU(:, [ndxmean ndxtails]), {'posterior mean', 'posterior 5\% quantile', 'posterior 95\% quantile'}, 'yyyymmm');

%% finish
dockAllFigures
