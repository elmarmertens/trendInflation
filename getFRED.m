%% construct FRED-based data sets
% uses input files obtained from FRED and FRB-PHIL

clear
clc

samStart = datenum(1960,1,1);
samEnd   = datenum(2021,9,1);

dates    = genrMdates(1960,year(samEnd),1);
dates    = dates(dates <= samEnd);
T        = length(dates);

datalabel = 'INFTRM';

%% process datalabel
switch datalabel
    case 'INF'
        Ylabel = {'PCE', 'PCEcore', 'CPI', 'GDPD'};
        mndx1   = [1 2 3];
        mndx2   = [3 4 1];
        qndx1   = 4;
    case 'INFTRM'
        Ylabel  = {'PCE', 'PCEcore', 'CPI', 'GDPD', 'PCEtrim', 'CPItrim', 'CPImedian'};
        mndx1   = [3 7 1 2 5 6];
        mndx2   = 1:6;
        qndx1   = 4;
    case 'INFTRM2'
        Ylabel  = {'PCE', 'PCEcore', 'CPI', 'GDPD', 'CPItrim', 'CPImedian'};
        mndx1   = [3 7 1 5 6];
        mndx2   = [1 2 3 5 6];
        qndx1   = 4;
    otherwise
        error('datalabel <<%s>> not known', datalabel)
end

Ny = length(Ylabel);

%% load FRED data

qd = importdata('INFTRM_Quarterly.txt');
md = importdata('INFTRM_Monthly.txt');

qd.dates = datenum(qd.textdata(2:end,1));
md.dates = datenum(md.textdata(2:end,1));

%% transform all series into annualized log-changes
qd.data(2:end,:) = diff(log(qd.data)) * 400;
qd.data(1,:)     = NaN;


levelNdx = [1 3 4];
ppNdx    = [2 5 6];

md.data(2:end,levelNdx) = diff(log(md.data(:,levelNdx))) * 1200;
md.data(1,levelNdx)     = NaN;

md.data(:,ppNdx)        = log(1 + md.data(:,ppNdx) / 100) * 100;



%% prep destintion data and copy FRED data

data = NaN(T,Ny);

m2dates = ismember(md.dates, dates);
data(ismember(dates, md.dates), mndx1) = md.data(m2dates,mndx2);

q2dates = ismember(qd.dates, dates);
data(ismember(dates, qd.dates), qndx1) = qd.data(q2dates,:);
% push quarterly data two months later
data(3:end,qndx1) = data(1:end-2,qndx1);
data(1:2,qndx1)   = NaN;

ynan = isnan(data);
data(ynan) = 0;

mat2fortran(sprintf('%s.dates.txt', datalabel), dates);
logical2fortran(sprintf('%s.yNaN.txt', datalabel), ynan);
mat2fortran(sprintf('%s.yData.txt', datalabel), data);

filename = sprintf('%s.settings.txt', datalabel);
fid = fopen(filename, 'wt');
fprintf(fid, 'Ny = %d\n', size(data,2));
fprintf(fid, 'T  = %d\n', size(data,1));
fprintf(fid, 'YLABEL:\n');
for n = 1 : Ny
    fprintf(fid, '%s\n', Ylabel{n});
end
fclose(fid);
display(filename);
type(filename)
hrulefill
