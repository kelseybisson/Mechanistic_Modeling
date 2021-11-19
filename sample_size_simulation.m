%            script to simulate what sample size is needed for model

%                   by: kelsey bisson, ucsb phd candidate
%                          -->  12 april 2017

% ________________________________________________________________________

clear all
%                                   steps

% 1.    randomly select lat /lons/ months to use in model given lat 
%       specified by satellite constraints. do this for different N sizes!! 
%               N = 5, 20, 100
% 2.    make it perfect data by running model with baseline parameters and
%       using the result as the data values. add 20 % noise to it. 

% 3.    then optimize the model with these different datasets and the 
%       baseline case as a priori values. choose base case,  BSA, BMA, BMSA

%% 0. select model run conditions 

howmany = 100;                       % how many bootstrap random samples?? 
nn = 5;                             % size of dataset
k = 1;                            % set noise level!!!

%% 0.5 set up bootstrap matrix and run time 

boot = nan(howmany,8); 
% cols are ll, (m-d/n), mean (m), mean (d) parameter #s
tic
for b = 1: howmany
%% 1. randomly select lat lons months & make perfect data


n = 4*nn;                            % overestimate data size first
lat = randi(65,[n 1]);
lon = randi(180,[n 1]);
sign1 = randi(2,[n 1]);              % 1 is positive, 2 is negative
lat(sign1==2) = -lat(sign1==2) ;     % correct lat according to sign 1
sign2 = randi(2,[n 1]);
lon(sign2==2) = -lon(sign2==2) ;     % correct lon according to sign 2
mo = randi(12,[n 1]);

%% 2. make data perfect & then add 20% noise to perfect data

% make data perfect 

load cafe.mat
load climAndTS_updated logChl mld t_an cbpm bbp z_eu vgpm Xi
clim  = {logChl mld t_an cbpm bbp z_eu vgpm CAFE Xi };

matches2 = matchmaker_kb(lat,lon,mo, clim);

% then pick first n that work (some random lat/lons/months wont work in
% model due to cloud cover /insufficient obs)

matches2(matches2(:,1)==0,:) =[];
par= matches2(1:nn,:);

% create poiiiifect data with baseline values

datas= ones(size(par,1),1);     % doesnt matter what data is for this purpose
o =1;               
[pos,l,m,d]=base_v18(log([0.1 0.1 0.3 0.1]),par,o,datas);
data =m;                            % the model outputs = perfect data

% now add 20% noise to it, where sample size is nn

v =  normrnd(k,k/10,1000,1);        % generate normal dist with mean 0.2 (20%)
r = randi(length(v),[nn 1]);         % randomly sample that normal distribution
sign = randi(2,[nn 1]);              % 1 is positive, 2 is negative
noise =v(r);
noise(sign==2) = -noise(sign==2);   % correct noise so some is negative
noisyd = data.* (1+noise);          % calculate new data with noise 
noisyd(noisyd<0) = 1;               % make negative data positive
m(m<0) = 1e-5;                      % make negative data positive

%% 3. Then optimize the model with these different data

%% COMPUTE SLICES FOR BASELINE MODEL 

o=1

for j =1:6
    
randn('state',sum(100*clock)); rand('state',sum(100*clock))
nsamples = 2000;

% fixed log-pdf target distribution
lpdf = @(p)base_v18(p,par,o,noisyd);

% initialize parameters
lb = [0; 0; 0; 0;];      % lower bounds
ub = [1; 1; 1; 1];  % upperbounds
x0 = random('unif',lb,ub);

% width of distribution
W = 0.5;

% gather parameter input samples
[sample,neval] = slicesample(log(x0),nsamples,'logpdf',lpdf,'width',W);

    for i = 1:nsamples
        p = sample(i,:);
    [pos(i),l,~,~]  = base_v18(p,par,o,noisyd);
    end

[R,neff,V,W,B] = psrf(pos');
R 

exp(sample);

if R > 0.998 && R < 1.04
    slice(:,:,j) = ans;
else
    slice(:,:,j) = NaN(2000,4);
end
end

if isreal(R) ==0
     slice(:,:,j) = -999.*ones(2000,4);
end

slices = [slice(1500:end,:,1); slice(1500:end,:,2); slice(1500:end,:,3);...
    slice(1500:end,:,4)];

a1= (nanmean(slices(:,1)));
a2= (nanmean(slices(:,2)));
a3= (nanmean(slices(:,3)));
a4= (nanmean(slices(:,4)));

p =[a1;a2;a3;a4];
[pos,l,m,d]=base_v18(log(p),par,o,noisyd);

boot(b,1) = l;
boot(b,2) = mean((m-d)/nn);
boot(b,3) = mean(m);
boot(b,4) = mean(d);
boot(b,5) = a1;
boot(b,6) = a2;
boot(b,7) = a3;
boot(b,8) = a4;

clear m d RMSE slices slice i b p

end
toc

disp(['ll =  ' num2str(nanmean(boot(:,1)))])
disp(['(m-d)/n =  ' num2str(nanmean(boot(:,2)))])
disp(['m =  ' num2str(nanmean(boot(:,3)))])
disp(['d =  ' num2str(nanmean(boot(:,4)))])



