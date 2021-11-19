% ========================================================================
%                           agg_v1.m

% code to predict aggregation globally as fn of concentration in waters

% critical Conc (ccr) defines max cell conc (threshold) in waters where agg
% occurs following coagulation theory. 
% OBJECTIVE: Compare Ccr to [chl] and predict where aggregation term is
% dominant. 

% all calc's done on 180 x 360 x 12 grid. 

% kelsey bisson,  7 july 2016
% phd student, university of california santa barbara

%       ---> updated 3 aug 2016

%                       Overview
% 1. calculate growth rate using satellite products 
% 2. calculate shear term wind velocities
% 3. define stickiness parameters
% 4. calculate critical concentration as f/n of satellite params
% 5. compare to satellite retrievals
% 6. and then do more stuff 
% ========================================================================

%% 1. Growth rate (u) = [2* chl/c * (1-e^-3Ig)] / [.022 + (.045-.022)e^-3Ig] 
 
load clims_4_agg.mat    
load kd_490.mat

carbon = 13000*(bbp - .00035); % from behrenfeld et al 2005
fraction = 10.^(logchl)./carbon; % chl/ c fraction

io = par;
ig3 = (3*io./24).*exp(-kd_490.*mld./2);  % moles photons m^-2 h^-1

u = 2*fraction .* (1-exp(-ig3))./ (.022 + (.045 - .022).*exp(-ig3));
u(u> 2) = NaN; u(u<0) = NaN; % divisions per day, remove unrealistic values

%                       ----------
%% 2.  calculate shear term wind velocities 
load wind.mat 
%                    variables 
a= 5.82e-6;     % kg m-3, constant term from MacKenzie & Leggett, '93
W = wind;       % windspeed, m s-1, 1988-2008 avg climatology
po = 1025;      % water density, kg m-3 
v = 1.0634e-6;  % m2 s-1, kinematic viscosity of seawater

E = a* (W.^3) ./mld;   % calcualate energy dissipation rate

g =( E./(1025*v)).^.5; % shear rate, in s^-1 
g = g* 86400;          % shear rate, in d^-1

%                       ----------
%% 3. Stickiness params
a1  = 0.1;             % low stickiness, 'other phyto's?'
a2  = 1;               % high stickiness, 'coccos' 'dias'
%                       ----------
%% 4. calculate critical concentration as f/n of satellite params

r = [.5e-6 20e-6 35e-6];     % .2 microns, .5 microns, 20 microns, 50 microns; 

for i=1:3; ccr_a1(:,:,:,i) = u./ (1.3 *a1 * g * 8 * (r(i)^3)); end
for i=1:3; ccr_a2(:,:,:,i) = u./ (1.3 *a2 * g * 8 * (r(i)^3)); end
% this represents cell counts.. 

Ccr08 = pi*u.*((g.*8).^-1); % critical concentration for all cell sizes considered, from 2008
% ^ this represents volumentric concentration, ppm

%                       ----------
%% 5. compare to satellite retrievals
load sizevars.mat

% create a function to evaluate the particle # within
% a certain size class. use kostadinov et al, 2010 algorith, to integrate N
% over the range of sizes considered

N0 = 10.^logNo;                      % convert log No to particle numbers 
fun = @(D) N0.* ((D/2e-6).^(-Xi));   % eqn 1 from Kostadinov et al 2010

q5_20 = integral(fun, .5e-6,20e-6,'ArrayValued',true);
q20_50 = integral(fun,20e-6,50e-6,'ArrayValued',true);

q5_20(isinf(q5_20))=NaN; q20_50(isinf(q20_50))=NaN;

%  calculate mean critical concentration for size class considered
cr5_20 = (ccr_a2(:,:,:,2) + ccr_a2(:,:,:,1)) ./3;
cr20_50= (ccr_a2(:,:,:,3) + ccr_a2(:,:,:,2)) ./2;

alpha = q5_20./ cr5_20; alpha(isinf(alpha))=NaN;
alpha2 = q20_50./cr20_50; alpha2(isinf(alpha2))=NaN;

an_q520= nanmean(q5_20,3); an_q2050= nanmean(q20_50,3);
an_cr5_20 = nanmean(cr5_20,3); an_cr2050 = nanmean(cr20_50,3);

an_alpha520 = an_q520./an_cr5_20;
an_alpha2050 = an_q2050./an_cr2050;

v_q520  = reshape(an_q520,[64800 1]); v_q2050  = reshape(an_q2050,[64800 1]); 
v_cr520 = reshape(an_cr5_20,[64800 1]); v_cr2050 = reshape(an_cr2050,[64800 1]); 

lat= 89.5:-1:-89.5; lon=-179.5:1:179.5; [Lon1,Lat1]= meshgrid(lon,lat);
v_lat = reshape(Lat1,[64800 1]);

load /Users/bisson/Documents/MATLAB/Homeless_CODE/ddd.mat

figure
subplot(221)
scatter(v_lat,v_q520./v_cr520,'k+'); axis([-70 70 0 2]);
title('C_{agg} vs latitude for 0.5-20um')
ylabel('Satellite # / # from theory')
xlabel('latitude')
set(gca,'fontsize',16);refline(0,1)

subplot(223)
imagesc(an_alpha520); colorbar;
caxis([0 .1]);
axis([0 360 10 165]);
title('C_{agg}')
colormap(ddd)
set(gca,'fontsize',16)

subplot(222)
scatter(v_lat,v_q2050./v_cr2050,'k+'); axis([-70 70 0 2]);
title('C_{agg} vs latitude for 20-50um ')
ylabel('Satellite # / # from theory')
xlabel('latitude')
set(gca,'fontsize',16); refline(0,1)

subplot(224)
imagesc(an_alpha2050); colorbar; 
caxis([0 .1]); 
axis([0 360 10 165]);
title('C_{agg}')
colormap(ddd)
set(gca,'fontsize',16)

% make plot of Ccr theory
r= linspace(0,50e-6,100); 
c = 1./(1.3*.1*8*(r.^3));

figure
plot(r,c,'k'); title('Critical concentration as f(r)');
xlabel('radius, m'); ylabel('Critical concentration, cells per m^3')
set(gca,'fontsize',16)

close figure 2

%% 6. Build index of aggregation on out

% what fraction of particles exceed the theory ? 

for i = 1:12;
high  = find(alpha(:,:,i)>1);
high2 = find(alpha2(:,:,i)>1);
end

% build ain index for aggregation fraction from the small portions! 
ind1 = reshape(alpha, [777600 1]);
ind1(ind1< .7) = .1;   
ind1(ind1> .7 & ind1< 1.3) = .2;
ind1(ind1>1.3) = .3; 

i_agg_small = reshape(ind1, [180 360 12]); 

% build ain index for aggregation fraction from the large portions! 
ind2 = reshape(alpha2, [777600 1]);
ind2(ind2< .7) = .1;   
ind2(ind2> .7 & ind1< 1.3) = .2;
ind2(ind2>1.3) = .3; 
i_agg_large = reshape(ind2, [180 360 12]); 

% ^ rationale is that theres a high and a low index -- low indices are .1
% and middle ones assigned a value of .2 and a high index is assigned .3

% Plot relationship between NPP and particle concentration 
% might expect highest NPP with higher particle concs, lower NPP at lower

load annualNPP.mat

conc_small = nanmean(q5_20,3);
conc_sm = reshape(conc_small, [64800 1]);
npp_v = reshape(avNPP, [64800 1]);

conc_big = nanmean(q20_50,3);
conc_b = reshape(conc_big, [64800 1]);

figure
subplot(121)
scatter(conc_sm,npp_v,'k','filled');
title('NPP vs particle concentration', 'fontsize',14)
xlabel('Small particle (5-20 um) concentration')
axis([0 2.5e12 0 2500]);
ylabel('NPP')

subplot(122)
scatter(conc_b,npp_v,'k','filled');
title('NPP vs particle concentration', 'fontsize',14)
xlabel('large particle (20-50 um) concentration')
ylabel('NPP')
axis([0 3.5e8 0 2500]);







