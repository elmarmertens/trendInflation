%% plot ucsv fortran results
initscript
if isunix
    initwrap
end

showSmoother = false;
showGains    = true;


datalabel = 'INFTRM';


mcmclabel = sprintf('notrendslopes.%s.T689', datalabel);

timestamp = [];
primaryNdx  = 1;

datadir   = pwd;

%% get parameters
if isempty(timestamp)
    filext = sprintf('particles.%s.gapSV.dat', datalabel);
else
    filext = sprintf('%s.particles.%s.gapSV.dat', timestamp, datalabel);
end

mcmcxt = sprintf('%s.gapSV.dat', mcmclabel);


%% get data

y     = importdata(fullfile(datadir, sprintf('YDATA.%s', filext)))';
yNaN  = logical(importdata(fullfile(datadir, sprintf('YNAN.%s', filext))))';
y(yNaN) = NaN;

foo = importdata(fullfile(datadir, sprintf('%s.settings.txt', datalabel)));
Ylabel = foo(4:end,1);

%% read results
T   = size(y,1);
Ny  = size(y,2);
Nstates = Ny * 2;
Nsv     = 1 + Ny;

dates = genrMdates(1960,2030);
dates = dates(1:T);

type(fullfile(datadir, strcat('settings.', filext)));

[TAU, GAP] = deal(NaN(T,12,Ny));
for s = 1 : Ny
    TAU(:,:,s)   = importdata(fullfile(datadir, sprintf('TAU%d.%s', s, filext)));
end
% for s = 1 : Ny
%     GAP(:,:,s)   = importdata(fullfile(datadir, sprintf('GAP%d.%s', s, filext)));
% end

[TAUhat, GAPhat] = deal(NaN(T,1,Ny));
for s = 1 : Ny
    TAUhat(:,:,s)   = importdata(fullfile(datadir, sprintf('TAUHAT%d.%s', s, filext)));
end

XBAR = importdata(fullfile(datadir, sprintf('XBAR.%s', filext)));

% for s = 1 : Ny
%     GAPhat(:,:,s)   = importdata(fullfile(datadir, sprintf('GAPHAT%d.%s', s, filext)));
% end

if showGains
    %     [GAINTAU, GAINGAP] = deal(NaN(T,12,Ny,Ny));
    %     for s = 1 : Ny
    %         for n = 1 : Ny
    %             GAINTAU(:,:,s,n)   = importdata(fullfile(datadir, sprintf('GAIN%dTAU%d.%s', n, s, filext)));
    %         end
    %         for n = 1 : Ny
    %             GAINGAP(:,:,s,n)   = importdata(fullfile(datadir, sprintf('GAIN%dGAP%d.%s', n, s, filext)));
    %         end
    %     end
    
    [GAINTAUhat, GAINGAPhat] = deal(NaN(T,Ny,Ny));
    for s = 1 : Ny
        GAINTAUhat(:,:,s)   = importdata(fullfile(datadir, sprintf('GAINTAUHAT%d.%s', s, filext)));
    end
    %     for s = 1 : Ny
    %         GAINGAPhat(:,:,s)   = importdata(fullfile(datadir, sprintf('GAINGAPHAT%d.%s', s, filext)));
    %     end
end

LOGLIKE         = importdata(fullfile(datadir, sprintf('LOGLIKE.%s', filext)));
ESS             = importdata(fullfile(datadir, sprintf('ESS.%s', filext)));

SV = NaN(T,12,Nsv);
for s = 1 : Nsv
    SV(:,:,s)   = importdata(fullfile(datadir, sprintf('SV%d.%s', s, filext)));
end

SVhat = NaN(T,1,Nsv);
for s = 1 : Nsv
    SVhat(:,:,s)   = importdata(fullfile(datadir, sprintf('SVHAT%d.%s', s, filext)));
end


if showSmoother
    [TAUsmoother, GAPsmoother] = deal(NaN(T,12,Ny)); %#ok<*UNRCH>
    for s = 1 : Ny
        TAUsmoother(:,:,s)   = importdata(fullfile(datadir, sprintf('smootherTAU%d.%s', s, filext)));
    end
    %     for s = 1 : Ny
    %         GAPsmoother(:,:,s)   = importdata(fullfile(datadir, sprintf('smootherGAP%d.%s', s, filext)));
    %     end
    
    SVsmoother = NaN(T,12,Nsv);
    for s = 1 : Nsv
        SVsmoother(:,:,s)   = importdata(fullfile(datadir, sprintf('smootherSV%d.%s', s, filext)));
    end
    
    SVmcmc = NaN(T, 12, Nsv);
    SVmcmc(:,:,1) = importdata(fullfile(datadir, sprintf('SIGTREND.%s', mcmcxt)));
    for n = 1 : Ny
        SVmcmc(:,:,1+n) = importdata(fullfile(datadir, sprintf('SIGGAP%d.%s', n, mcmcxt)));
    end
    
end

% read mcmc level
TAUmcmc = NaN(T,12,Ny);
for n = 1 : Ny
    TAUmcmc(:,:,n) = importdata(fullfile(datadir, sprintf('TAU%d.%s', n, mcmcxt)));
end
    
% parameters
hgap0 = importdata(fullfile(datadir, sprintf('HGAP0.%s', mcmcxt)));

ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 5 6 8]; % used for the smoother
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
horizons = [1, 4, 8, 12, 16];


%% filter: TAU and data
for i = 1 : Ny
    newfigure(sprintf('tau %d -- %s', i, Ylabel{i}))
    set(gcf, 'Renderer', 'painters')
    % plotCI(TAU(:,ndxmean,i), TAU(:,ndxtails,i), dates)
    plotCI(TAUhat(:,i), TAU(:,ndxtails,i), dates)
    plot(dates, TAU(:,ndxmean,i), 'r--')
    hold on
    
    plotynoncompact(dates, y(:,i), [1 0 0]);
    plothorzline(2, [], 'k-');
    ylim([-2 16])
    nberlines(dates)
    plotOrigin
    title(Ylabel{i})
    wrapcf(sprintf('tau%s', Ylabel{i}), wrap)
end

%% TAU: compare actual filter against "bar" filter
for i = 1 : Ny
    newfigure
    set(gcf, 'Renderer', 'painters')
    hold on
    plot(dates,TAUhat(:,i), 'b-', 'linewidth', 2)
    plot(dates,XBAR(i,:), 'r--', 'linewidth', 2)
    plothorzline(2, [], 'k-');
    ylim([-2 16])
    nberlines(dates)
    plotOrigin
    title(Ylabel{i})
    wrapcf(sprintf('taubar%s', Ylabel{i}), wrap)
end

%% TAU: compare particle filter against mcmc smoother
for i = 1 : Ny
    newfigure
    set(gcf, 'Renderer', 'painters')
    hold on
    plot(dates,TAUhat(:,i), 'b-', 'linewidth', 2)
    plot(dates,TAUmcmc(:,ndxmean,i), 'r--', 'linewidth', 2)
    plothorzline(2, [], 'k-');
    ylim([-2 16])
    nberlines(dates)
    plotOrigin
    title(Ylabel{i})
    wrapcf(sprintf('tauPFvsMCMC%s', Ylabel{i}), wrap)
end


%% filter: check gap plus trend
checkdiff(y,squeeze(TAU(:,1,:)+GAP(:,1,:))); % does not have to hold exactly since they are draws independently of each other
checkdiff(y,squeeze(TAUhat  + GAPhat)); % does not have to hold exactly since they are draws independently of each other

checkdiff(TAUhat,squeeze(TAU(:,1,:)));
checkdiff(GAPhat,squeeze(GAP(:,1,:)));

%% GAINS: TAU
if showGains
    % gains on primary trend
    primaryGAIN = GAINTAUhat(:,:,primaryNdx);
    if ~all(vec(primaryGAIN(yNaN) == 0))
        error('check gains')
    end
    primaryGAIN(yNaN) = NaN;
    
    %     for j = 1 : Ny
    %         newfigure(sprintf('primary gain of %s', Ylabel{j}))
    %         hold on
    %         set(gcf, 'Renderer', 'painters')
    %         % plot(dates, primaryGAIN(:,j), 'b-', 'linewidth', 2)
    %         bar(dates, primaryGAIN(:,j), 'Linestyle', '-', 'BarWidth', 1, ...
    %             'edgecolor', [1 0 0], 'facecolor', [1 0 0], 'linewidth', 1)
    %         ylim([-.3  1])
    %         set(gca, 'ytick', -1 : .25 : 2)
    %         xtickdates(dates)
    %         wrapcf(sprintf('primarygain%stau%s', Ylabel{j}), wrap)
    %     end
    
    
    %% monthly pictures
    %     for t = find(year(dates) == 2010, 1, 'first') : find(year(dates) == 2010, 1, 'last')
    %         figure
    %         bar(1 : Ny, primaryGAIN(t,:), 'Linestyle', '-', 'BarWidth', .8, ...
    %             'edgecolor', [1 0 0], 'facecolor', [1 0 0], 'linewidth', 1)
    %         title(datestr(dates(t), 'mmmm yyyy'))
    %         set(gca, 'xticklabel', Ylabel)
    %         %     ylim([-.1 .3])
    %         set(gca, 'ytick', -1 : .1 : 2)
    %         wrapcf(sprintf('GAIN%s', lower(datestr(dates(t), 'mmmyyyy'))), wrap)
    %     end
    
    

    %% stacked with negative values
    % close all
    %     for thisYear = 1990 : 2014
    %         ndxYear = find(year(dates) == thisYear, 1, 'first') : find(year(dates) == thisYear, 1, 'last');
    %         figure
    %
    %         if all(all(primaryGAIN(ndxYear,:) > 0))
    %             bar(dates(ndxYear), primaryGAIN(ndxYear,:), 'stacked', 'linewidth', 2)
    %         else
    %             hold on
    %
    %             % positive values
    %             theseValues = primaryGAIN(ndxYear,:);
    %             theseValues(theseValues > 0) = NaN;
    %             bar(dates(ndxYear), theseValues, 'stacked', 'linewidth', 2)
    %
    %             % negative values
    %             theseValues = primaryGAIN(ndxYear,:);
    %             theseValues(theseValues < 0) = NaN;
    %             bar(dates(ndxYear), theseValues, 'stacked', 'linewidth', 2)
    %
    %         end
    %         limmerick = ylim;
    %         ylim([limmerick(1) limmerick(2) * 1.7])
    %         set(gca, 'xtick', dates(ndxYear))
    %         xlim([dates(ndxYear(1) - 1)  dates(ndxYear(end) + 1)])
    %         datetick('x', 'mmmm', 'keeplimits', 'keepticks')
    %         legend(Ylabel{:}, 'location', 'northeast')
    %         title(sprintf('Year %d', thisYear))
    %         colormap bone
    %         wrapcf(sprintf('GAIN%d', thisYear), wrap)
    %     end
    
    %% one month all years
    
    monthLabels = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
    activeNdx   = any(~isnan(y),1);
    edgecolor = [1 1 1];
    for thisMonth = 1 :  12
        ndx = find(month(dates) == thisMonth & year(dates) >= 1965);
        figure
        thisGAIN = primaryGAIN(ndx,activeNdx);
        if all(all(primaryGAIN(ndx,:) > 0))
            bar(dates(ndx), thisGAIN, 'stacked', 'linewidth', 1, 'EdgeColor', edgecolor)
        else
            hold on
            
            % positive values
            theseValues = thisGAIN;
            theseValues(theseValues > 0) = NaN;
            bar(dates(ndx), theseValues, 'stacked', 'linewidth', 1, 'EdgeColor', edgecolor)
            
            % negative values
            theseValues = thisGAIN;
            theseValues(theseValues < 0) = NaN;
            bar(dates(ndx), theseValues, 'stacked', 'linewidth', 1, 'EdgeColor', edgecolor)
            
        end
        %     limmerick = ylim;
        %     ylim([limmerick(1) limmerick(2) * 1.2])
        set(gca, 'xtick', dates(1 : 60 : end))
        xlim([dates(ndx(1) - 12)  datenum(2016,1,1)])
        datetick('x', 'yyyy', 'keeplimits', 'keepticks')
        legend(Ylabel{activeNdx}, 'location', 'northeast')
        title(sprintf('%s', monthLabels{thisMonth}))
        colormap jet
        wrapcf(sprintf('GAINmonth%d', thisMonth), wrap)
        
        %% try surface
        %         %         thisGAIN = primaryGAIN(ndx,activeNdx);
        %         thisGAIN = primaryGAIN(:,activeNdx);
        %         thisGAIN(isnan(thisGAIN)) = 0;
        %         %         figure
        %         %         surface(thisGAIN)
        %         %         shading interp
        %
        %         figure
        %         bar3(thisGAIN, 1)
        %         legend(Ylabel{activeNdx})
        
    end
end

%% ESS
newfigure('ESS')
plot(dates,ESS)
ylim([0 1])
nbershades(dates)
wrapcf('ESS', wrap)

%% SV
for s = 1 : Nsv


    
    newfigure(sprintf('SV %d', s))
    set(gcf, 'Renderer', 'painters')
    hold on
    plotCI(SV(:,ndxmean,s), SV(:,ndxtails,s), dates)
    if showSmoother
        plotCIlines(SVsmoother(:,ndxmean,s), SVsmoother(:,ndxtails,s), dates, [], 'r')
    end
    %     plot(dates, SV(:,ndxmean,s), 'r-', 'linewidth', 2)
    xlim(dates([1 end]))
    ylimmer = ylim;
    ylim([0 ylimmer(2)])
    datetick('x', 'keeplimits')
    grid on
    wrapcf(sprintf('SV%d', s), wrap)
    
end

checkdiff(SVhat,squeeze(SV(:,1,:)));

%% LIKELIHOOD
newfigure('loglike')
plot(dates,cumsum(LOGLIKE) ./ (1:T)')
title(sprintf('LLF = %4.8f --- LLF / T = %4.8f', sum(LOGLIKE), mean(LOGLIKE)))
nbershades(dates)
wrapcf('loglike', wrap)


%% SMOOTHER X
if showSmoother
    % hanni = NaN(3,1);
    for i = 1 : Ny
        newfigure(sprintf('tau %d -- %s', i, Ylabel{i}))
        set(gcf, 'Renderer', 'painters')
        hold on
        %     hanni(1) = plotCIlines(TAUsmoother(:,ndxmean,i), TAUsmoother(:,ndxtails,i), dates, [], [1 0 0]);
        %     hanni(2) = plot(dates, TAU(:,ndxmean,i), 'k--', 'linewidth', 3);
        
        
        hanni(1) = plotCI(TAUhat(:,i), TAU(:,ndxtails,i), dates);
        hanni(2) = plotCIlines(TAUsmoother(:,ndxmean,i), TAUsmoother(:,ndxtails,i), dates, [], [1 0 0]);
        
        
        %     hanni(3) = plotCIlines(mcmc.YGAP(:,ndxmean), mcmc.YGAP(:,ndxtails), dates, [], [0 .5 0]);
        %     plot(dates, y(:,i), 'm-')
        %     if ~iscompact(y(:,i))
        %         plot(dates, y(:,i), 'mx')
        %         plot(dates, y(:,i), 'mo')
        %     end
        ylim([-2 16])
        plothorzline(2, [], 'k-');
        nberlines(dates)
        plotOrigin
        title(Ylabel{i})
        wrapcf(sprintf('SMOOTHERtau%s', Ylabel{i}), wrap)
    end
    
    %% compare smoother vs filter vs mcmc
    
    hanni = NaN(3,1);
    for i = 1 : Ny
        newfigure(sprintf('tau %d -- %s', i, Ylabel{i}))
        set(gcf, 'Renderer', 'painters')
        hold on
        hanni(1) = plot(dates, TAUsmoother(:,ndxmean,i), 'b-', 'linewidth', 3);
        hanni(2) = plot(dates, TAUhat(:,i), 'r--', 'linewidth', 3);
        hanni(3) = plot(dates, TAUmcmc(:,ndxmean,i), 'K-.', 'linewidth', 3);
        nberlines(dates)
        plotOrigin
        title(Ylabel{i})
        wrapcf(sprintf('COMPAREtau%s', Ylabel{i}), wrap)
    end
    
    %% SV comparison
    hanni = NaN(3,1);
    SVgap0     = mean(exp(hgap0 / 2));
    SVgap0tail = prctile(exp(hgap0 / 2), [25 75]);
    
    for i = 1 : Nsv
        newfigure(sprintf('SV %d ', i))
        set(gcf, 'Renderer', 'painters')
        hold on
        hanni(1) = plot(dates, SV(:,ndxmean,i), 'r--', 'linewidth', 3);
        hanni(2) = plot(dates, SVsmoother(:,ndxmean,i), 'b-', 'linewidth', 3);
        hanni(3) = plot(dates, SVmcmc(:,ndxmean,i), 'K-.', 'linewidth', 3);
        if i > 1
            plothorzline(SVgap0(i-1), [], 'k-.', 'linewidth', 2)
            plothorzline(SVgap0tail(1,i-1), [], 'k-.', 'linewidth', 1)
            plothorzline(SVgap0tail(2,i-1), [], 'k-.', 'linewidth', 1)
        end
        
        nberlines(dates)
        plotOrigin
        if i > 1
            shadeUntilFirstObs(dates, yNaN(:,i-1))
            title(Ylabel{i-1})
        end
        legend(hanni, 'Filter', 'Smoother', 'MCMC', 'location', 'best')
        wrapcf(sprintf('COMPAREsv%d', i), wrap)
    end
    
end



%% finish
finishwrap
dockAllFigures
finishscript
