% EZflux_T_v1.m

% First written by Dave Siegel, modified by Kelsey Bisson, 2014

load climAndTS_updated

% The four adjustable model parameters...
falg = 0.1;       % Fraction of NPPM going into sinking aggregates
m_phyto = 0.1;    % Specific mortality by viral lysis and other biological/nongrazing losses
ffecM = 0.3;      % Fraction of GrzM going into flux
ffecN = 0.1;      % Fraction of GRZN going into flux

% Set up the interpolated time axes
% doy range (0 to 366 by 1's) is chosen so n=1 is doy=0 which also is doy=366 
% and n=2 is doy=1 and doy=367  - need to enforce the wrap around BC's
DelT = 1;
doy = 0:DelT:367;
ndays = length(doy);

% day4mon = central day of each month starting in Dec and going to Jan (14 months)
day4mon = [-45 -15 15 46 74 105 135 166 196 227 258 288 319 349 380 410]; 
monday = day4mon(3:14); % central doy for each 12 months
for i=1:12; idmon(i) = find(doy == monday(i)); end; idmon = idmon';

% Set parameters for Npts for obs and the lat/lon ranges
nobsmin = 8;
nlat = 180; nlon = 360; %90N-90S, 180W-180E

% Create 2D arrays for NPP, AlgEZ, FecEZ & TotEZ
avNPP = NaN(nlat,nlon); NPPmv = NaN(nlat,nlon,12); 
avAlgEZ = avNPP; avmonAlgEZ = NaN(nlat,nlon,12); 
avFecEZ = avNPP; avmonFecEZ = NaN(nlat,nlon,12); 
avTotEZ = avNPP; TotEZmv = NaN(nlat,nlon,12); 
fmicrom = avmonAlgEZ; 

ngood = zeros(nlat,nlon);

% Loop over time series...
for ilat = 1:nlat
for ilon = 1:nlon
    
% Test to see if there is a time series at this point
test1 = squeeze(10.^logChl(ilat,ilon,:))'; 
test2 = squeeze(mld(ilat,ilon,:))'; 
test3 = squeeze(t_an(ilat,ilon,1,:))'; 
test4 = squeeze(cbpm(ilat,ilon,:))'; 
test5 = squeeze(bbp(ilat,ilon,:))'; 
test6 = squeeze(z_eu(ilat,ilon,:))'; 

ngood(ilat,ilon) = length(find(test1 > 0 & test1 < 20 & test2 > 0 & ...
    ~isnan(test3) & test4 > 0 & test4 < 5000 & test5 > 0 & test5 < 0.01 ...
    & test6 > 0 & test6 < 250 ));

if (ngood(ilat,ilon) >= nobsmin); % Sets the mininum number of observations for each location

% Unpack the variables...
BBP = squeeze(bbp(ilat,ilon,:))';         BBPm = [BBP(11:12) BBP BBP(1:2)];
NPP = squeeze(cbpm(ilat,ilon,:))';        NPPm = [NPP(11:12) NPP NPP(1:2)];
XXi = squeeze(Xi(ilat,ilon,:))';          Xim =  [XXi(11:12) XXi XXi(1:2)];
zeu = squeeze(z_eu(ilat,ilon,:))';        zeum = [zeu(11:12) zeu zeu(1:2)];
zml = squeeze(mld(ilat,ilon,:))';         zmlm = [zml(11:12) zml zml(1:2)];
chl = squeeze(10.^logChl(ilat,ilon,:))';  chlm = [chl(11:12) chl chl(1:2)];

% Interpolate the seasonal time seris to a finer delT (= 1 day)
BBP = interp1(day4mon,BBPm,doy,'spline'); 
NPP = interp1(day4mon,NPPm,doy,'spline'); 
XXi = interp1(day4mon,Xim,doy,'spline'); 
zeu = interp1(day4mon,zeum,doy,'spline'); 
zml = interp1(day4mon,zmlm,doy,'spline'); 
chl = interp1(day4mon,chlm,doy,'spline'); 

% Calculate fmirco from Xi following Kostadinov et al 2010 using
% Dmax = 50um, Dmin = 0.5um & Micros from 20um to 50um
% Dspl = 20; Dmax = 50; Dmin = 0.5;
Dspl = 20; Dmax = 50; Dmin = 0.5;
if XXi == 4;
    fmicro = log(Dmax/Dspl)/log(Dmax/Dmin);
else
    fmicro = (Dmax.^(4-XXi) - Dspl.^(4-XXi)) ./ (Dmax.^(4-XXi) - Dmin.^(4-XXi));
end

% -------------------------------------------------------------------
% Calculate biomass fractions (assumes BBP method or Chl with a C/Chl=50)
PC = (BBP - 0.00035)*13000; % from Mike's GBC paper
% PhytoC = chl * 50;

% Algal flux (AlgEZ = falg(NPPm) * NPPm   
AlgEZ = falg .* NPP;

% Calculate dPM/dt & dPN/dt
dXXX = diff(PC)./DelT; 
dZZZ = diff(zml)./DelT;
    doymod = doy(1:length(dXXX))+DelT/2;
dPCdt = interp1(doymod,dXXX,doy,'linear','extrap'); 
dMLDdt = interp1(doymod,dZZZ,doy,'linear','extrap');
HdMLDdt = ones(size(dMLDdt)); 
iff = find(dMLDdt <= 0); 
if(~isempty(iff)); HdMLDdt(iff) = zeros(size(iff)); end
entrn = dMLDdt .* HdMLDdt;

% Estimate phytoplankton grazing rates in units of mgC/m3/day
Grz = NPP./zeu - dPCdt - AlgEZ./zeu - m_phyto*PC - (PC./zml).*entrn;
    iff = find(Grz < 0); if(~isempty(iff)); Grz(iff) = zeros(size(iff)); end

% Fecal flux (FecEZ = (ffecM*GrzM + ffecN*GrzN)*zeu)  
FecEZ = (ffecM*Grz) .* zeu;

if (max(Grz) > 0 )
    
% Stuff the stuff in variables...
avNPP(ilat,ilon) = mean(NPP); 
TotEZ = FecEZ + AlgEZ;
avTotEZ(ilat,ilon) = mean(TotEZ); 

end
end
end
end


















