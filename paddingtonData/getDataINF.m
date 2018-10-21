%% construct data set INFTRM for Paddington 2015
% uses input files obtained from FRED

clear
clc

samStart = datenum(1960,1,1);
samEnd   = datenum(2018,09,1);

dates    = genrMdates(1960,year(samEnd),1);
dates    = dates(dates <= samEnd);
T        = length(dates);


Ylabel = {'PCE', 'PCEcore', 'CPI', 'GDPD'};
Ny = length(Ylabel);

dataLabel = 'INF';

%% load FRED data

qd = importdata('INFTRM_Quarterly.txt');
md = importdata('INFTRM_Monthly.txt');

qd.dates = datenum(qd.textdata(2:end,1));
md.dates = datenum(md.textdata(2:end,1));

%% transforma ll series into annualized log-changes
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
mndx    = [3 4 1]; 
data(ismember(dates, md.dates), 1:3) = md.data(m2dates,mndx);

q2dates = ismember(qd.dates, dates);
qndx    = 4; 
data(ismember(dates, qd.dates), qndx) = qd.data(q2dates,:);
% push quarterly data two months later
data(3:end,qndx) = data(1:end-2,qndx);
data(1:2,qndx)   = NaN;

ynan = isnan(data);
data(ynan) = 0;
     
mat2fortran(sprintf('%s.dates.txt', dataLabel), dates);
logical2fortran(sprintf('%s.yNaN.txt', dataLabel), ynan);
mat2fortran(sprintf('%s.yData.txt', dataLabel), data);

filename = sprintf('%s.settings.txt', dataLabel);
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
