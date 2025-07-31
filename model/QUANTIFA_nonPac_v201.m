%% Quantile Analysis of Temperature using Individual Foraminiferal Analyses
% Tropical Pacific Ocean Variant (version 2.0.1)
%
% ---- About ----     
%
% QUANTIFA is a user-friendly IFA proxy system model that combines routines
% for modeling the sensitivity of IFA populations to changes in climate 
% variability with tools for processing, plotting, and interpreting IFA-
% Mg/Ca data. Please see the README document on GitHub (https://github.com/
% rh-glaubke/QUANTIFA) for a quickstart guide on how to use this algorithm. 
%
% This is the second version of QUANTIFA (v2.0.1) and comes with all new 
% features, including:
%
% (1) a new "ENSO variability" metric that manipulates both ENSO amplitude 
% and frequency;
% (2) new weighting routines for depth habitat and seasonal bias inputs
% that allow users to (a) define a species' depth habitat as a range of
% depths, rather than a single depth horizon, and (b) weight the picking 
% algorithm toward a species' preferred season of productivity, rather than
% picking from those months exclusively;
% (3) a function that provides a rough quantitative measure of past climate
% variability by computing the barycenter (or "center of gravity") of the
% data-model consistency maps; and
% (4) other various performance and bug fixes.
% 
% For a more detailed description of these new features, please refer to
% the references below.
% 
% QUANTIFA is designed to work with potential temperature data from the  
% Ocean Reanalysis System 5 data assimilation (ORA-S5). This script variant
% is designed to work with data from the tropical Atlantic (TA_ORAS5.mat)
% or tropical Indian Ocean (TI_ORAS5.mat). Both of these datasets can be, 
% downloaded from the QUANTIFA repository on GitHub. For any information
% regarding the ORA-S5 dataset, please refer to Zuo et al. (2019) 
% (doi:10.5194/os-15-779-2019).
%
% ---- Author Info and Citation ----
% 
% Written by Ryan H. Glaubke (glaubke@arizona.edu) and Kaustubh Thirumalai 
% (kaustubh@arizona.edu)
% Latest Update: July 2025
% 
% Citations: 
% v2.0.1 — Glaubke et al. (2024). An inconsistent ENSO response to Northern
%          Hemisphere stadials over the Last Deglaciation. Geophysical
%          Research Letters. doi:10.1029/2023gl107634
%
% v1.0.1 — Glaubke et al. (2021). Discerning changes in high-frequency
%          climate variability using geochemical populations of individual 
%          foraminifera. Paleoceanography and Paleoclimatology, 36(2), 
%          e2020PA004065. doi:10.1029/2020PA004065
%
%% Input Window
% Enter data and relevant input parameters below, then execute algorithm.
% Information regarding each parameter is available below:
%
% Individual Foraminiferal Datasets:
%     X - Reference Population (e.g. core-top population)
%     Y - Comparison Population (e.g. down-core population)
%
%     NOTE: If you would like to run a detection sensitivity analysis,
%     please comment out both X and Y below. If you would like to compare
%     an IFA population against reanalysis data, please assign IFA data to
%     the Y variable and comment out X.
%
% Reanalysis Dataset:
%     The ORAS5 .mat file is a workspace file containing a gridded field of
%     potential temperature data ('votemper'; 1°x1° horizontal resolution) 
%     for 75 analyzed depth levels ('depths'). It is available for download 
%     here: https://github.com/rh-glaubke/QUANTIFA
%
% Location of Core:
%     lat - Latitude for temperature gridpoint (Format: XX.5)
%     lon - Longitude for temperature gridpoint (Format: XX)
%     dep - Mean depth and standard deviation for species of interest 
%           (Format: [XX,YY])
%
% Calibration Equation:
%     Select from a small bank of pre-loaded Mg/Ca calibration equations
%     (listed below), or program your own within the algorithm. Includes
%     the option to incorporate a dissolution correction term (D).
%
% Input Parameters:
%     A set of input conditions used to parameterize the modeled pseudoproxy
%     time series and initialize the picking scheme. A breif description of 
%     each term can be found below.
%
% Output Products:
%     QUANTIFA produces the following outputs:
%     1) A "conformity contour plot" showing the relative influence of
%        annual and interannual climate variability on quantiles in the 
%        body and tails of the paleotemperature distributions.
%     2) A Q-Q plot of the user's data (or the user's data vs. ORA-S5 data).
%     3) An array of false positive rates for each individual quantile 
%        (i.e. the probability each quantile commites a type I error).
%     4) A data-model consistency map that illustrates the proportion of
%        significant quantiles that align with each hypothetical climate
%        scenario. Darker colors = scenarios most consistent with the user's
%        data. Also includes barycenter coordinates.
%     
%     There are also a suite of supplemental figures, including:
%     1) The psuedoproxy time series representing modern variability 
%        plotted alongside time series of select altered climate scenarios.
%     2) Idealized quantile plots for select altered climate scenarios 
%        relative to modern variability.
% 
% FOR THE TROPICAL PACIFIC DATASET:
% Horizontal Domain:  [29.5°N-29.5°S, 120°E-70°W; Convention: -ve for S,W]
% Vertical Domain:    [0.5-5902 m; See 'depths']

clc;
clear;
close all;
tic

% ---- Individual Foraminiferal Datasets ----

%X = 'Insert X Population Here';
%Y = 'Insert Y Population Here';

% ---- Reanalysis Dataset ----

load('Scripts/TA_ORAS5.mat');   % Call up ORA-S5 dataset 
                                % (download link for .mat file above)

% ---- Core Location and Foraminiferal Ecology ----

lat = 11.5;                     % Coordinates of sediment core
lon = -67;                      % (example: VM12-107; Caribbean Sea)

dep = [0,0];                    % Foram depth habitat [m,std]
norm = 0;                       % 1 = implement as normal distribution
                                % 0 = implement as consistent depth window
                                
% ---- Calibration Equation ----

eqn = 1;                        % 1 = Anand et al. (2003) All species eqn
%D = 2.8;                       % 2 = Anand et al. (2003) G. ruber eqn
                                % 3 = Dekens et al. (2002) G. ruber eqn
                                %     (w/ dissolution correction term, D,
                                %     as core depth in km)
                                % 4 = Anand et al. (2003) G. sacculifer eqn
                                % 5 = Dekens et al. (2002) N. dutertrei eqn
                                %     (w/ dissolution correction term, D,
                                %     as core depth in km)
                                % 6 = Mohtadi et al. (2011) G. tumida eqn
                                % 7 = Anand et al. (2003) P. obliq eqn
                                
% ---- Input Parameters ----

% Model Initialization
run = 5000;                    % No. of Monte Carlo iterations
fp = 100;                      % No. of False Positive Exercises
num_p = 70;                    % No. of "forams" picked from time series
num_q = 35;                    % No. of quantiles calculated
tsl = 1000;                    % Length of synthetic time series (yrs)

% Uncertainties
cl = 90;                       % Confidence level (80, 90, 95, or 99)
anerr = 0.05;                  % Analytical error (mmol/mol)
calerr = 0;                    % Calibration error (°C) [OPTIONAL]

% Weighting Parameters
seas = [0,0];                  % Seasonal bias for picking scheme [m,std] 
                               % (J-D = 1-12; [0,0] = no bias)
fva = 1;                       % Relationship between ENSO freq and amp
                               % (1 = freq and amp are manipulated equally)
                               
% NOTE: You MUST manually input the years of known positive and negative
% phase events for the interannual climate mode you are interested in
% (i.e. the Atlantic Nino or the Indian Ocean Dipole). Also, you will want
% to change the starting month for the event year QUANTIFA averages over 
% (line XXX). As a placeholder, we have entered El Nino and La Nina years
% since 1958, and assigned the start month to April.

pyear = [1963;1965;1968;1972;1982;1986;1987;1991;1994;1997;2002;2009;2015];
nyear = [1970;1971;1973;1975;1984;1988;1995;1998;1999;2000;2007;2010;2011];

% ---- Output Products ----
% Note: All materials are saved to your working directory

saveoutput = 0;                % 1 = Save model results (.mat)
                               % 0 = Do not save model results
suppfig = 0;                   % 1 = Generate supplemental figures
                               % 0 = Do not generate supplemental figures
savefig = 0;                   % 1 = Save figures (PDF)
                               % 0 = Do not save figures
                             
                                
% ========================================================================
% ================= Modify Sections Below with Care! =====================
% ======================================================================== 


%% Input Parsing
% Checks input variables to ensure proper formatting.

disp('------------------------------------------------------------------------');
disp('Quantile Analysis of Temperature using Individual Foraminiferal Analyses');
disp('(QUANTIFA) v2.0.1 - Tropical Atlantic and Indian Variant');
disp(' ');
disp('Written by Ryan Glaubke and Kaustubh Thirumalai'); 
disp('Latest Update: July 2025');
disp('Latest Reference: Glaubke et al. (2024), Geophysical Research Letters');
disp('doi:10.1029/2023gl107634');
disp('------------------------------------------------------------------------');

% ---- Initiate Progress Bar ----

multiWaitbar('Close All');
multiWaitbar('Initializing QUANTIFA...', 0,'Color',[0.44,0.75,0.96]);

% ---- Check Inputs ----

% -- IFA Data --                    

if exist('X','var') == 0 && exist('Y','var') == 0
    disp(' ');
    disp('NOTE: No IFA data detected. Running sensitivity tests...');
    NoIFA = 1;
elseif exist('X','var') == 0 && exist('Y','var') == 1
    disp(' ');
    disp('NOTE: No IFA Data assigned to X. Using ORA-S5 data...');
    NoIFA = 0.5;
elseif exist('X','var') == 1 && exist('Y','var') == 0
    disp(' ');
    error('Error: No IFA Data assigned to Y. If comparing IFA Data to ORA-S5, please assign IFA Data to Y.');
else
    NoIFA = 0;
end

% -- Core Location and Depth --

% Check latitude formatting and domain
if rem(lat,1) == 0                  
    disp(' ');
    error('Error: Latitude is a whole number. Please input latitude as XX.5');
end

if (lat > 29.5) || (lat < -29.5)
    disp(' ');
    error('Error: Entered latitude is out of model domain [29.5°N - 29.5°S]');
end

% Check longitude formatting and domain
if rem(lon,1) ~= 0                 
    disp(' ');
    error('Error: Longitude not a whole number. Please revise.');
end

if exist('TA_ORAS5', 'var')         % Check model domain: longitude
    if lon < 0                          
        if (360+lon > 375) || (360+lon <= 259)
            disp(' ');
            error('Error: Entered longitude is out of model domain [100°W - 15°E].');
        end
    else
        if (lon > 15 && lon <= 259) 
            disp(' ');
            error('Error: Entered longitude is out of model domain [100°W - 15°W]');
        end
    end
elseif exist('TI_ORAS5', 'var')
    if (lon > 120) || (lon < 40)
        disp(' ');
        error('Error: Entered longitude is out of model domain [40°E - 120°E]');
    end
end

% Find nearest depth levels to inputted range
nrst_depth = [];                    
udep = dep(1) - dep(2);
ldep = dep(1) + dep(2);

[~,uidx] = min(abs(depths-udep)); [~,lidx] = min(abs(depths-ldep));
nrst_depth(1,:) = uidx; nrst_depth(2,:) = lidx;
dep_range = depths(nrst_depth(1):nrst_depth(2));

disp(' ');
disp(['NOTE: Entered depth range = ',num2str(udep),' - ',num2str(ldep),' m']);
disp(['      Nearest depth levels = ',num2str(dep_range(1)),' - ',num2str(dep_range(end)),' m']);
disp(' ');

% -- Mg/Ca-T Equation --

% Check for dissolution correction term
if eqn == 3 || eqn == 5                   
    if exist('D','var') == 0
        error('Error: no dissolution correction term found! Please assign the correction term native to your selected calibration equation under the variable "D".');
    end
end

% -- Input Parameters --

% run
if run > 1000
    disp('NOTE: At >1000 Monte Carlo iterations, QUANTIFA can take some');
    disp('      time depending on the performance of your machine.');
    disp(' ');
end

% seas (and initialize bias)
if length(seas) < 2
    error('Error: Unrecognized seasonal bias. Please enter a mean and standard deviation. If you prefer no seasonal bias, please enter [0,0].');
end

if seas(1) > 12
    error('Error: Unrecognized seasonal bias. Please enter 1 - 12 in the first position to represent January - December.');
end

if seas == [0,0]
    seas_wt = ones(tsl,1);                      
elseif seas(1) + seas(2) < 12
    mths = [1:12];
    seas_wt = normpdf(mths,seas(1),seas(2))';
    seas_wt = repmat(seas_wt,tsl,1);
elseif seas(1) + seas(2) > 12                   
    mths = [1:1:seas(1)+6];
    seas_wt = normpdf(mths,seas(1),seas(2))';
    for i = 1:length(find(mths>12))
        seas_wt(mths == i) = seas_wt(mths == 12+i);
    end
    seas_wt = seas_wt(1:12);
    seas_wt = repmat(seas_wt,tsl,1);
end
  
% calerr
if exist('calerr','var') == 0
    disp('NOTE: No calibration error detected. Uncertainty estimation will');
    disp('      not take calibration error into account.');
    disp(' ');
    calerr = 0;
end

% anerr
if anerr == 0
    error('Error: No analytical error detected. Please assign a value. If you are conducting a sensitivity analysis, a general value between 0.05 and 0.10 mmol/mol should be sufficient.');
end

% cl
if (cl ~= 80) && (cl ~= 90) && (cl ~= 95) && (cl ~= 99)
    error('Error: Unrecognized confidence level. Please enter either 80, 90, 95, or 99.');
end

%% Instrumental Data Extraction
% Extracts temperature time series data from ORA-S5 (Zuo et al. [2019])
% from specified lat/lon/dep parameters.

% ---- ORA-S5 lat/lon/depth indicies ----

[ilat,ilon,idepth] = oraco(lat,lon,nrst_depth,depths);

% ---- Time constraints ----

lyear = 2019;                   % Last year of time series
years = 61;                     % Length of time series: 1958(Jan)-2018(Dec)
syear = lyear - years;          % Starting year of time series
len = years*12;                 % No. of months
las = (lyear-1870)*12;          % Last month
fir = las - len + 1;            % First month
tmslice = tsl*12;               % Time slice in months

% ---- Declaration ----

T = zeros(len,1);               % Temperature data
M = zeros(len,1);               % Months within temperature time series

ys = syear;
c = 0;
for i=1:len                     % Converting the unit of time to years
    if (c<=11)
        M(i) = ys + (c/12);
        c = c+1;
    else
        c = 0;
        ys = ys+1;
        M(i) = ys + (c/12);
        c = c+1;
    end
end
% ---- Extract ORA-S5 Data ----

for i = 1:length(dep_range)
    
    depth = dep_range(i);

    % -- Index and Extract Temperature Time Series --

    for i = 1:len
        if exist('TA_ORAS5', 'var')
            [ilat,ilon,idepth] = ta_oraco(lat,lon,depth,depths);
            T(i) = TA_ORAS5.votemper(ilon,ilat,idepth,i);
        elseif exist('TI_ORAS5', 'var')
            [ilat,ilon,idepth] = ti_oraco(lat,lon,depth,depths);
            T(i) = TI_ORAS5.votemper(ilon,ilat,idepth,i);
        end
    end

    if isempty(t(isnan(t))) == 0    % Replace terminal NaNs by splicing
        le = length(t(isnan(t)));
        t(len-(le-1):len) = t(1:le);
    end

    % NOTE: only a few locations within the model domain have terminal NaNs
    % remedied by the code above. If using a different reanalysis dataset, 
    % this will need to be adjusted accordingly.
    
    T(:,i) = t;
    
end

% ---- Time Series Check ----

if nnz(~isnan(T)) == 0
   disp(' ');
   error('Error: Land Ho! It looks like theres no temperature data in your indexed time series... maybe youve struck land? Try another set of coordinates.');
end

%% Foraminiferal Mg/Ca Forward Model
% Transfers temperature record into Pseudo-Mg/Ca space.

prox = zeros(len,size(T,2));
switch eqn
    case 1  % Anand et al. (2003) all planktic species eqn
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.38*exp(0.09*T(i,t));
            end
        end        
    case 2  % Anand et al. (2003) G. ruber eqn
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.34*exp(0.102*T(i,t));
            end
        end
    case 3  % Dekens et al. (2002) G. ruber eqn (w/ dissolution correction)
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.38*exp(0.09*(T(i,t)-(0.61*D)-1.6));
            end
        end
    case 4  % Anand et al. (2003) G. sacculifer eqn
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 1.06*exp(0.048*T(i,t));
            end
        end   
    case 5  % Dekens et al. (2002) N. dutertrei eqn (w/ dissolution correction)
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.60*exp(0.08*(T(i,t)-(2.8*D)-5.4));
            end
        end 
    case 6  % Mohtadi et al. (2011) G. tumida eqn
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.410*exp(0.068*T(i,t));
            end
        end 
    case 7  % Anand et al. (2003) P. obliquiloculata eqn
        for t = 1:size(T,2)
            for i = 1:len
                prox(i,t) = 0.18*exp(0.12*T(i,t));
            end
        end       
end

% ---- Average Pseudoproxy Time Series ----

if norm == 1
    
    % -- Define weights --
    
    dep_weight = normpdf(dep_range,dep(1),dep(2));
    
    % -- Compute weighted average and depth habitat uncertainty --
    
    MgCa = sum(prox.*dep_weight,2)/sum(dep_weight);
    sdMgCa = std(prox,dep_weight/sum(dep_weight),2);
    
    derr = mean(sdMgCa);
    
elseif norm == 0

    % -- Compute unweighted average and depth habitat uncertainty --
    
    MgCa = mean(prox,2);
    sdMgCa = std(prox,[],2);
    
    derr = mean(sdMgCa);
    
end

%% Climatology
% Extracts El Niño, La Niña, and neutral year climatologies (beginning with 
% the month of April) from forward-modeled Mg/Ca time series.
% 
% If your time series is longer than 1958-present, you need to manually 
% add in El Niño years and La Niña years for the record. Alternatively, 
% you can code a function that defines ENSO events giving SST anomalies
% within the Nino3.4 region.

% ---- Define ENSO Event Years and Months ----

pen = (pyear-syear)*12 + 4;             % + 4 => Begin w/ April
nen = (nyear-syear)*12 + 4;

pin = (pyear - syear) + 1;              % Year for injecting ENSO event
nin = (nyear - syear) + 1;

beg_event = sort(vertcat(pen,nen),'ascend');     
event_mths = zeros(1,length(beg_event)*12);

for i = 1:length(beg_event)
    event_mths((i-1)*12 + 1) = beg_event(i);
    for j = 1:11
        next_num = beg_event(i) + j;
        event_mths((i-1)*12 + j + 1) = next_num;
    end
end

event_mths = event_mths';

% ---- Declaration ----

mp = zeros(12,length(pyear));
mn = zeros(12,length(nyear));

Cp = zeros(12,1);
Cpsd = zeros(12,1);

Cn = zeros(12,1);
Cnsd = zeros(12,1);

CC = zeros(12,1);
CCsd = zeros(12,1);

% ---- Yearly Data Extraction ----

% -- Positive Interannual Events --

for i = 1:length(pen)
    mp(:,i) = MgCa(pen(i):(pen(i)+11));
end
pt = reshape(mp,length(pyear)*12,1);

% -- Negative Interannual Events --

for i = 1:length(nen)
    mn(:,i) = MgCa(nen(i):(nen(i)+11));
end
nt = reshape(mn,length(nyear)*12,1);

% -- Mean ENSO Climatologies --

for month=1:12
    Cp(month) = nanmean(pt(month:12:end));
    Cn(month) = nanmean(nt(month:12:end));

    Cpsd(month) = std(pt(month:12:end));
    Cnsd(month) = std(nt(month:12:end));
end

% -- Mean Base Climatology (Seasonality) --

m_idx = (month:12:length(MgCa))';
ev_idx = intersect(m_idx,event_mths);
nt = MgCa; 
nt(event_mths,:) = [];

l = 1;
for month=4:15                           
    CC(l) = nanmean(nt(month:12:end));
    CCsd(l) = std(nt(month:12:end));
    l = l + 1;
end

% -- Interannual Deviations from Base Climatology --

dp = Cp - CC;
dn = Cn - CC;

% -- Interannual Relationship Check --

if (nanmean(Cp) < nanmean(Cn))        % eg. in W Atlantic, Atl Niño = cold
    w = 0;                            
else
    w = 1;
end

%% Synthetic Time Series: Modeled Modern Variability
% Constructs synthetic time series mimicking modern variability. Red and
% white noise is added and altered according to a while loop check with
% corresponding tolerances.

[rratai,rfw,rpwf,rwp,rrp] = powrat(MgCa);  % Function that extracts
                                           % spectral parameters (w/ noise)

% ---- Annual:Interannual Power Tolerance ----

aitol = 30;             % Percent tolerance of annual:interannual
mdif = 2*aitol;         % Initial condition
mag = 1;                % Magnitude of ENSO event to ENSO event variability

% ---- Noise Model Parameters ----

mnt = 0.75;             % Multiplication component
rnt = 0.25;             % Initial condition for red noise parameter
wnt = 0.25;             % Initial condition for white noise parameter

% ---- Noise Model Tolerance ----

wth = 30;               % White noise threshold
rth = 30;               % Red noise threshold
wnp = 2*wth;
rnp = 2*rth;

% ---- Modeled Time Series ----

xx = 0;
SCheck = zeros(6,5);

% ---- Establish Base Climatology (Seasonality) ----

mod_clim = zeros(12,years);         
for i = 1:years
    mod_clim(:,i) = CC(:);
end

while ((mdif > aitol) || (wnp > wth) || (rnp > rth))

    xx = xx + 1;

    % ---- Insert Positive Interannual Events ----

    for i = 1:length(pin)
        if (w == 0)
            mod_clim(:,pin(i)) = Cp - rand*Cpsd*mag;
        else
            mod_clim(:,pin(i)) = Cp + rand*Cpsd*mag;
        end
    end
    
    % ---- Insert Negative Interannual Events ----
    
    for i = 1:length(nin)
        if (w == 0)
            mod_clim(:,nin(i)) = Cn + rand*Cnsd*mag;
        else
            mod_clim(:,nin(i)) = Cn - rand*Cnsd*mag;
        end
    end

    elMgCa1 = reshape(mod_clim,len,1);

    % ---- Mimicking Reality ----

    fel = elMgCa1;

    % -- Noise Addition --

    for i=2:len
        fel(i) = mnt*fel(i) + rnt*fel(i-1) + wnt*(rand);
    end

    % ---- Spectral Parameters ----

    [modratai,mfw,mpwf,mwp,mrp] = powrat(fel);
    mdif = (abs(modratai-rratai)/rratai)*100;

    % ---- Noise Parameters ----

    rnp = abs((mrp-rrp)/rrp*100);
    wnp = abs((mwp-rwp)/rwp*100);

    % ---- Storage ----

    SCheck(xx,1) = mdif;
    SCheck(xx,2) = wnt;
    SCheck(xx,3) = wnp;
    SCheck(xx,4) = rnt;
    SCheck(xx,5) = rnp;

    % ---- Tweaking Noise ----

    if (rnp > rth && ((mrp-rrp)/rrp*100) < -50)   % Red noise check
        rnt = rnt + 0.05;
    elseif (rnp > rth && ((mrp-rrp)/rrp*100) > 50)
        rnt = rnt - 0.05;
    end

    if (wnp > wth && ((mwp-rwp)/rwp*100) < -50)   % White noise check
        wnt = wnt + 0.05;
    elseif (wnp > wth && ((mwp-rwp)/rwp*100) > 50)
        wnt = wnt - 0.05;
    end
end

% ---- Construct Full-length Modern Time Series ----

mod_series = zeros(tmslice,1);

if (tmslice ~= len)
    for i = 1:tmslice
        mod_series(i) = elMgCa1(mod(i,len)+1);    % Full Synthetic Series
    end
else
    mod_series = elMgCa1;
end

% ---- Noise Model Tolerance ----

wnt = SCheck(xx,2);
rnt = SCheck(xx,4);

for i = 2:tmslice
    mod_series(i) = mnt*mod_series(i) + rnt*mod_series(i-1) + wnt*(rand);
end

%% Synthetic Time Series: Modeled Past Variability
% Constructs 440 synthetic time series with altered annual/interannual
% variability. Seasonal Cycle amplitude and interannual amplitude+frequency
% vary from 100% dampened (or no signal) to 100% amplified (or double signal). 
% 21 hypothetical seasonal cycle scenarios x 21 hypothetical interannual 
% scenarios = 441 total scenarios (440 altered and 1 modern).

% ---- Declaration ----

alt_series = zeros(21,21,tmslice);      % Altered Time Series

% ---- Simulate Altered Time Series ----

r = 1;                          
for p = -100:10:100                         % Seasonal Cycle Amplitude
    c = 1;
    for w = -100:10:100                     % Interannual Variability

        % ---- Change Seasonal Cycle ----

        sc = (CC - mean(CC))*(1+(p/100));   % Change Seasonal Cycle
        nCC = mean(CC) + sc;                % Modified Seasonal Cycle

        % ---- Change ENSO Variability ----
        
        if w == 0                           % Change Interannual Variability:
            amp = (w/100)+1;                % Adjust Amplitude and Frequency
            fre = ((w/100) * fva)+1;        % according to prescribed weight
        else
            amp = (w/100)+1;
            fre = ((w/100) * fva)+1;
        end
        
        nCp = CC + amp.*dp;                % Modified Positive Event Amp
        nCn = CC + amp.*dn;                % Modified Negative Event Amp

        pfre = (length(pin)/years);        % Modified Positive Event Freq
        pnfre = pfre * fre;                % Modified Negative Event Freq
        pstep = years * pnfre;     

        nfre = (length(nin)/years);        % Modern La Niña Frequency
        nnfre = nfre * fre;                % Modified La Niña Frequency
        nstep = years * nnfre;

        npin = round(linspace(1,years-1,pstep));    % New Positive Events
        nnin = round(linspace(2,years,nstep));      % New Negative Events

        % ---- Insert Altered Base Climatology (Seasonality) ----

        alt_clim = zeros(12,years);
        for i = 1:years
            alt_clim(:,i) = (nCC(:));
        end

        % ---- Insert Altered ENSO Events ----

        for i = 1:length(npin)             % Insert Altered ENSO Events
            if (w == 0)
                alt_clim(:,npin(i)) = nCp - rand*Cpsd*mag;
                alt_clim(:,nnin(i)) = nCn + rand*Cnsd*mag;
            else
                alt_clim(:,npin(i)) = nCp + rand*Cpsd*mag;
                alt_clim(:,nnin(i)) = nCn - rand*Cnsd*mag;
            end
        end

        % ---- Construct Full-length Altered Time Series ----

        elMgCa2 = reshape(alt_clim,len,1);
        as = zeros(tmslice,1);

        if (tmslice ~= len)
           for i = 1:tmslice
               as(i) = elMgCa2(mod(i,len)+1);
           end
        else
           as = elMgCa2;
        end

        for i=2:tmslice
            as(i) = mnt*as(i) + rnt*as(i-1) + wnt*(rand);
        end
        alt_series(r,c,:) = as;
        c = c +1;
    end
    r = r + 1;
end

%% Quantile and Residual Computation: Modeled Data
% Creates idealized Q-Q patterns for each modeled climate scenario. First, 
% it iteratively subsamples the time series, translates to paleotemperature
% values, centers the data, then computes quantiles and residuals. These
% quantiles and residuals are then binned and averaged.

multiWaitbar('Initializing QUANTIFA...',0,...
             'Relabel','Picking Pseudo-IFA Data...');

% ---- For Modeled Modern Time Series ----

mod_qs = zeros(num_q,1);

% -- Initiate Picking Scheme --

tot = num_p * run;                              % Define total picks
picks = datasample(mod_series,tot, ...          % Pick Mg/Ca values
        'Weight',seas_wt);                      % Implement seasonal bias
picks = picks + anerr * randn(size(picks)) ...  % Add error terms
        + derr * randn(size(picks));
MC = reshape(picks,[num_p,run]);                % Pseudo-IFA datasets

% -- Compute Centered Quantiles --

switch eqn
    case 1
        mc_mod = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
    case 2
        mc_mod = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
    case 3
        mc_mod = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
            + calerr*randn(size(MC));
    case 4
        mc_mod = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
    case 5
        mc_mod = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
            + calerr*randn(size(MC));
    case 6
        mc_mod = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
    case 7
        mc_mod = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
end

for i = 1:run
    q = quantile(mc_mod(:,i)-mean(mc_mod(:,i)),num_q);
    mod_qs(:,i) = q;
end
m_mod_qs = mean(mod_qs,2);

% ---- For Altered Time Series ----

alt_qs = zeros(21,21,num_q,run);    % Quantiles of each altered series
yr_alt = zeros(21,21,num_q,run);    % Y Residuals b/w quantiles and 1:1 line

pb = 0;
for p = 1:21
    for w = 1:21
       
        as = squeeze(alt_series(p,w,:));
        
        % -- Initiate Picking Scheme --

        tot = num_p * run;                              % Define total picks
        picks = datasample(as,tot, ...                  % Pick Mg/Ca values
                'Weight',seas_wt);                      % Implement seasonal bias
        picks = picks + anerr * randn(size(picks)) ...  % Add error terms
                + derr * randn(size(picks));
        MC = reshape(picks,[num_p,run]);                % Pseudo-IFA datasets
        
        switch eqn
            case 1
                mc_alt = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
            case 2
                mc_alt = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
            case 3
                mc_alt = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
                    + calerr*randn(size(MC));
            case 4
                mc_alt = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
            case 5
                mc_alt = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
                    + calerr*randn(size(MC));
            case 6
                mc_alt = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
            case 7
                mc_alt = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
        end
        
        % -- Compute Centered Quantiles --
        
        for i = 1:run
            q = quantile(mc_alt(:,i)-mean(mc_alt(:,i)),num_q);
            alt_qs(p,w,:,i) = q;
        end
        
        % -- Calculate Residuals --

        for i = 1:num_q
            asq = squeeze(alt_qs(p,w,:,:));
            dist = asq(i,:) - mod_qs(i,:);
            yr_alt(p,w,i,:) = dist;
        end
        
        % -- Update Progress Bar --
        
        pb = pb + 1;
        if NoIFA == 0
            multiWaitbar('Picking Pseudo-IFA Data...',pb/882); % 0% --> 50%
        elseif NoIFA == 1
            multiWaitbar('Picking Pseudo-IFA Data...',pb/485); % 0% --> 90%
        end
        
    end
end

multiWaitbar('Picking Pseudo-IFA Data...',...   
             'Relabel','Generating Modeled Results...');         

if NoIFA == 0
    prog = pb/882;
elseif NoIFA == 1
    prog = pb/485;
end
    
% ---- Generate Bins ----

bins = quantile(m_mod_qs-mean(m_mod_qs),num_q)';                
d = diff(bins)/2;                               % Bins defined by quantiles
                                                % of MMV
edges = zeros(1,length(bins)+1);                % Bin edges are mid-way b/w
                                                % MMV quantiles
for i = 2:length(bins)
    edges(:,i) = bins(i-1,:) + d(i-1,:);
end
edges(:,1) = min(mod_qs(1,:));
edges(:,end) = max(mod_qs(num_q,:));

% ---- Bin Data ----

m_alt_qs = zeros(21,21,num_q);      % Mean quantiles for each comparison
m_resids = zeros(21,21,num_q);      % Corresponding mean residuals
sd_alt_qs = zeros(21,21,num_q);     % Standard deviation within each bin

pb = 0;
for p = 1:21                        
    for w = 1:21
        for i = 1:num_q
            
            % -- Identify Binned Quantiles --
            
            if i < num_q
                idx = find(mod_qs>=edges(:,i) & mod_qs<edges(:,i+1));
            else
                idx = find(mod_qs>=edges(:,i) & mod_qs<=edges(:,i+1));
            end
            
            % -- Calculate Mean Quantile --
            
            A = squeeze(alt_qs(p,w,:,:));
            alt_vals = A(idx);
            m = mean(alt_vals);
            m_alt_qs(p,w,i) = m;
            
            % -- Calculate SD --
            
            sd = std(alt_vals);
            sd_alt_qs(p,w,i) = sd;
            
            % -- Calculate Mean Residual --
            
            B = squeeze(yr_alt(p,w,:,:));
            dalt_vals = B(idx);
            m = mean(dalt_vals);
            m_resids(p,w,i) = m;
            
        end
        
        % -- Update Progress Bar --
        
        pb = pb + 1;
        multiWaitbar('Generating Modeled Results...',prog+pb/8820);
        
    end
end

prog = prog+pb/8820;

%% Modeling IFA Detection Sensitivity
% Tests the sensitivity of the body and the tails of paleotemperature 
% distributions (at the user-defined location) to changes in annual and 
% interannual climate variability.

% ---- Declaration ----

prb = round(num_q*0.16);
ct = (1:prb);                       % Cold Tail = Lower 16% of Distribution
wt = (num_q-prb:num_q);             % Warm Tail = Upper 16% of Distribution
bdy = (prb(end)+1:wt(1)-1);         % Body = Center 68% of Distribution

lval = zeros(21,21,num_q);          % Logical Values Matrix
ctp = zeros(21,21);                 % Cold Tail Sensitivity
bdyp = zeros(21,21);                % Body Sensitivity
wtp = zeros(21,21);                 % Warm Tail Sensitivity

% ---- Sensitivity Test ----

pb = 0;
for p = 1:21
    for w = 1:21
        for i = 1:num_q

            % -- Index Residuals --

            if i <= prb
                if (m_resids(p,w,i) >= sd_alt_qs(p,w,i)*-1) && (m_resids(p,w,i) <= sd_alt_qs(p,w,i))
                    lval(p,w,i) = 1;
                else
                    lval(p,w,i) = 0;
                end              
            elseif (i >= prb(end)+1) && (i <= wt(1)-1)
                if (m_resids(p,w,i) >= sd_alt_qs(p,w,i)*-1) && (m_resids(p,w,i) <= sd_alt_qs(p,w,i))
                    lval(p,w,i) = 1;
                else
                    lval(p,w,i) = 0;
                end
            elseif i >= num_q-prb
                if (m_resids(p,w,i) >= sd_alt_qs(p,w,i)*-1) && (m_resids(p,w,i) <= sd_alt_qs(p,w,i))
                    lval(p,w,i) = 1;
                else
                    lval(p,w,i) = 0;
                end                
            end
        end

        % -- Percent Conformity w/ Modern Time Series --

        for i = 1:num_q
            if i <= prb
                pctg = (sum(lval(p,w,1:prb))/length(ct))*100;
                ctp(p,w) = pctg;
            elseif (prb(end)+1 >= i) && (i <= wt(1)-1)
                pctg = (sum(lval(p,w,prb(end)+1:wt(1)-1))/length(bdy))*100;
                bdyp(p,w) = pctg;
            elseif i >= num_q-prb
                pctg = (sum(lval(p,w,num_q-prb:num_q))/length(wt))*100;
                wtp(p,w) = pctg;
            end
        end
        
        % -- Update Progress Bar --
        
        pb = pb + 1;
        multiWaitbar('Generating Modeled Results...',prog+pb/8820);
        
    end
end

prog = prog+pb/8820;

%% Output Product #1: Conformity Contour Plot
% Produces three contour plots of % conformity values for each altered
% climate scenario w/ respect to modeled modern time series. Tests the
% sensitivty of the tails (upper and lower 16%) and the central body
% (middle 68%) of the distribution to changes in ENSO and seasonality.

if NoIFA == 1
    multiWaitbar('Generating Modeled Results...',...
                 'Relabel','Finishing Up...',...
                 'Increment',0.40);
end

% ---- Generate Figure ----

f1 = figure(1);
clf;

pp = (-100:10:100)';
ww = (-100:10:100)';

for i = 1:3
    
    % -- Index Region --
    
    if i == 1
        pd = ctp;
    elseif i == 2
        pd = bdyp;
    elseif i == 3
        pd = wtp;
    end
    
    % -- Plot Contours --
    
    subplot(1,3,i)
    [cs,h]=contourf(ww,pp,pd,[0:5:100],'linestyle','none');
    hold on
    [cs2,h2]=contour(ww,pp,pd,[25,50,75,100],'linewidth',0.35);
    plot(0,0,'*y','MarkerSize',12);
    
    % -- Define Color Map --
    
    if i == 1
        title(gca,'Cold Tail Sensitivity (Bottom 16%)','fontsize',14);
        h2.Color = [0 0.2 0.4];
        clabel(cs2,h2,'FontName','Myriad Pro','fontsize',10,'color',[0 0.2 0.4],'LabelSpacing',250);
        
        co = [255 255 255;
            243 250 252;
            209 238 251;
            151 182 221;
            88 126 188;
            61 109 177;]./255;
        
        vs = [0; 10; 30; 70; 90; 100];
        cmap = interp1(vs./100,co,linspace(0,1,100));
        colormap(gca,cmap)
        colorbar
        
    elseif i == 2
        title(gca,'Body Sensitivity (Middle 68%)','fontsize',14);
        h2.Color = [0.3 0 0.4];
        clabel(cs2,h2,'FontName','Myriad Pro','fontsize',10,'color',[0.3 0 0.4],'LabelSpacing',250);
        
        co = [255 255 255;
            240 220 252;
            235 214 252;
            191 119 246;
            175 97 235;
            146 71 206;]./255;
        
        vs = [0; 10; 30; 70; 90; 100];
        cmap = interp1(vs./100,co,linspace(0,1,100));
        colormap(gca,cmap)
        colorbar
        
    elseif i == 3
        title(gca,'Warm Tail Sensitivity (Upper 16%)','fontsize',14);
        h2.Color = [0.6 0 0];
        clabel(cs2,h2,'FontName','Myriad Pro','fontsize',10,'color',[0.6 0 0],'LabelSpacing',250);
        
        co = [255 255 255;
            254 244 244;
            247 180 183
            239 151 158;
            238 81 85;
            235 59 63]./255;
        
        vs = [0; 10; 30; 70; 90; 100];
        cmap = interp1(vs./100,co,linspace(0,1,100));
        colormap(gca,cmap)
        colorbar
    end
    
    % -- Figure Properties --
    
    xlim([-100 100]);
    ylim([-100 100]);
    caxis([0 100]);
    xlabel('Change in Interannual Variability (%)','FontName','Myriad Pro','fontsize',14,'Fontweight','demi');
    ylabel('Change in Annual Variability (%)','FontName','Myriad Pro','fontsize',14,'Fontweight','demi');
    set(gca,'LineWidth',1.5,'fontsize',14,'Fontweight','demi');
    set(gca,'Xtick',-100:20:100)
    set(gca,'XtickLabel',[-100:20:100])
    
end
f1.GraphicsSmoothing = 'on';

% -- Save as PDF --

if savefig == 1
    f_filename = 'QUANTIFA_ConformityContourPlot';
    export_fig(f_filename, '-pdf', '-append', '-transparent', '-r600');
end

% ---- Terminate Script ----

if NoIFA == 1
    multiWaitbar('Finishing Up...','Close');
    disp('Run Complete!');
    toc
    return
end

%% False Positive Rate Estimation
% Computes the probability of committing a type I error for each individual
% quantile by repeatedly subsampling the modeled modern time series and
% comparing the two pseudo-IFA populations against one another via Q-Q
% analysis.

multiWaitbar('Generating Modeled Results...',...    
             'Relabel','Estimating False Positives...');

% ---- Declaration ----

fpr = zeros(num_q, fp);             

if exist('X','var') == 0            
    X = T;
    numX = num_p;
else                                   
    numX = length(X);
end
numY = length(Y);

% ---- For False Positive Exercise ----

pb = 0;
for kk = 1:fp
    
    fp_qx = zeros(num_q,run);          % Quantiles from Pseudo-IFA
    fp_qy = zeros(num_q,run);          % Quantiles from Pseudo-IFA

    fpm_qy = zeros(num_q,1);           % Mean Y Quantile
    fp_serr = zeros(num_q,1);          % Total Sampling Error (X and Y)
    
    % ---- Initiate Picking Scheme ----

    % -- For X Population --

    tot = numX * run;                               % Define total picks
    picks = datasample(mod_series,tot, ...          % Pick Mg/Ca values
            'Weight',seas_wt);                      % Implement seasonal bias
    picks = picks + anerr * randn(size(picks)) ...  % Add error terms
            + derr * randn(size(picks));
    MC = reshape(picks,[numX,run]);                % Pseudo-IFA datasets

    switch eqn
        case 1
            mc_x = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
        case 2
            mc_x = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
        case 3
            mc_x = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
                + calerr*randn(size(MC));
        case 4
            mc_x = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
        case 5
            mc_x = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
                + calerr*randn(size(MC));
        case 6
            mc_x = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
        case 7
            mc_x = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
    end

    for i = 1:run
        q = quantile(mc_x(:,i)-mean(mc_x(:,i)),num_q);
        fp_qx(:,i) = q;
    end

    % -- For Y Population --

    tot = numY * run;                               % Define total picks
    picks = datasample(mod_series,tot, ...          % Pick Mg/Ca values
            'Weight',seas_wt);                      % Implement seasonal bias
    picks = picks + anerr * randn(size(picks)) ...  % Add error terms
            + derr * randn(size(picks));
    MC = reshape(picks,[numY,run]);                % Pseudo-IFA datasets
    
    switch eqn
        case 1
            mc_y = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
        case 2
            mc_y = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
        case 3
            mc_y = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
                + calerr*randn(size(MC));
        case 4
            mc_y = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
        case 5
            mc_y = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
                + calerr*randn(size(MC));
        case 6
            mc_y = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
        case 7
            mc_y = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
    end

    for i = 1:run
        q = quantile(mc_y(:,i)-mean(mc_y(:,i)),num_q);
        fp_qy(:,i) = q;
    end

    % ---- Generate Bins ----

    b = mean(fp_qx,2);                % <-- Doesn't matter if I choose X or
                                      %     Y (both are subsamples of MMV)

    bins = quantile(b-mean(b),num_q)';                
    d = diff(bins)/2;                 % Bins defined by quantiles of MMV
    edges = zeros(1,length(bins)+1);  % New edges from resampled modern TS

    for i = 2:length(bins)
        edges(:,i) = bins(i-1,:) + d(i-1,:);
    end
    edges(:,1) = min(fp_qx(1,:,1));
    edges(:,end) = max(fp_qx(num_q,:,1));

    % ---- Bin Data ----

    for i = 1:num_q
        if i < num_q
            idx = find(fp_qx>=edges(:,i) & fp_qx<edges(:,i+1));
        else
            idx = find(fp_qx(:,:,1)>=edges(:,i) & fp_qx<=edges(:,i+1));
        end

        qy_vals = fp_qy(idx);
        m = mean(qy_vals);
        fpm_qy(i,1) = m;
        sig = std(qy_vals);
        fp_serr(i,1) = sig;

    end

    % ---- Calculate Residuals ----

    yr_data = zeros(num_q,run);     % Vert distance b/w quantiles & 1:1 line
    xr_data = zeros(num_q,run);     % Horiz distance b/w quantiles & 1:1 line

    for i = 1:num_q
        for ii = 1:run
            yr = fp_qy(i,ii) - fp_qx(i,ii);           
            xr = fp_qx(i,ii) - fp_qy(i,ii);           
            yr_data(i,ii) = yr;
            xr_data(i,ii) = xr;
        end
    end

    % ---- Find Significant Quantiles ----

    lval = zeros(num_q,run);        % 0 = not significant (at 90% CI)
                                    % 1 = significant (at 90% CI)

    % -- ID Significant Quantiles --

    for i = 1:num_q
        for ii = 1:run  
            if abs(yr_data(i,ii)) > fp_serr(i)*1.645
                if abs(xr_data(i)) > fp_serr(i)*1.645
                    lval(i,ii) = 1;
                else
                    lval(i,ii) = 0;
                end
            else
                lval(i,ii) = 0;
            end
        end
    end
    
    % -- Calculate False Positive Rate --

    for i = 1:num_q
        fpr(i,kk) = (sum(lval(i,:))/length(lval(i,:)))*100;            
    end
        
    % -- Update Progress Bar --
    
    pb = pb + 1;
    multiWaitbar('Estimating False Positives...','Increment',...
                 (0.35/fp));
   
end

% ---- Output False Positive Rates (Mean and SD) ----

fpr_m = mean(fpr,2);                % <--- mean FPR across fp
fpr_sd = std(fpr,[],2);             % <--- 1 sigma across fp

% NOTE: The way this bloc is constructed is very much like the bloc just 
% below it (constraining uncertainties for user data). Therefore, some of
% the operative variables here are overwritten as the code progresses. If
% you want to preserve this information, you should change the variables.
% In any case, these variables are erased in the cleaned workspace at the
% end of the model run. Bottom line: if all you want are the FPR's, then
% you're fine.

%% Uncertainty Estimation: Error on Quantiles from User Data
% Estimates total uncertainty for inputted IFA data in X and Y via a smooth
% bootstrap approach. Sampling uncertainty is estimated by repeatedly 
% resampling the modeled modern time series. Estimation takes analytical 
% and calibration certainty into account.

multiWaitbar('Estimating False Positives...',... 
             'Relabel','Finishing Up...',...
             'Increment',0.01);

% ---- Declaration ----

mc_qx = zeros(num_q,run,2);         % Quantiles from Pseudo-IFA
mc_qy = zeros(num_q,run,2);         % Quantiles from Pseudo-IFA

m_qx = zeros(num_q,1);              % Mean X Quantiles
m_qy = zeros(num_q,1);              % Mean Y Quantiles

serrX = zeros(num_q,1);             % One sigma total uncertainty for X
serrY = zeros(num_q,1);             % One sigma total uncertainty for Y

% ---- Picking Algorithm and Uncertainty Estimation ----

% -- For X Population --

tot = numX * run * 2;                           % Define total picks
picks = datasample(mod_series,tot, ...          % Pick Mg/Ca values
        'Weight',seas_wt);                      % Implement seasonal bias
picks = picks + anerr * randn(size(picks)) ...  % Add error terms
        + derr * randn(size(picks));
MC = reshape(picks,[numX,run,2]);              % Pseudo-IFA datasets

switch eqn
    case 1
        mc_x = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
    case 2
        mc_x = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
    case 3
        mc_x = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
            + calerr*randn(size(MC));
    case 4
        mc_x = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
    case 5
        mc_x = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
            + calerr*randn(size(MC));
    case 6
        mc_x = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
    case 7
        mc_x = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
end

for ii = 1:2
    for i = 1:run
        q = quantile(mc_x(:,i,ii)-mean(mc_x(:,i,ii)),num_q);
        mc_qx(:,i,ii) = q;
    end
end

% -- For Y Population --

tot = numY * run * 2;                           % Define total picks
picks = datasample(mod_series,tot, ...          % Pick Mg/Ca values
        'Weight',seas_wt);                      % Implement seasonal bias
picks = picks + anerr * randn(size(picks)) ...  % Add error terms
        + derr * randn(size(picks));
MC = reshape(picks,[numY,run,2]);              % Pseudo-IFA datasets

switch eqn
    case 1
        mc_y = (log(MC/0.38)/0.09) + calerr*randn(size(MC));
    case 2
        mc_y = (log(MC/0.34)/0.102) + calerr*randn(size(MC));
    case 3
        mc_y = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) ...
            + calerr*randn(size(MC));
    case 4
        mc_y = (log(MC/1.06)/0.048) + calerr*randn(size(MC));
    case 5
        mc_y = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
            + calerr*randn(size(MC));
    case 6
        mc_y = (log(MC/0.41)/0.068) + calerr*randn(size(MC));
    case 7
        mc_y = (log(MC/0.18)/0.12) + calerr*randn(size(MC));
end

for ii = 1:2
    for i = 1:run
        q = quantile(mc_y(:,i,ii)-mean(mc_y(:,i,ii)),num_q);
        mc_qy(:,i,ii) = q;
    end
end

% ---- Generate Bins ----

edges = zeros(1,length(bins)+1,2);  % New edges from resampled modern TS

for ii = 1:2
    for i = 2:length(bins)
        edges(:,i,ii) = bins(i-1,:) + d(i-1,:);
    end
    if ii == 1
        edges(:,1,ii) = min(mc_qx(1,:,1));
        edges(:,end,ii) = max(mc_qx(num_q,:,1));
    else
        edges(:,1,ii) = min(mc_qy(1,:,1));
        edges(:,end,ii) = max(mc_qy(num_q,:,1));
    end
end

% ---- Bin Data ----

for ii = 1:2                        % Calculate uncertainty for X data
    if ii == 1
        for i = 1:num_q
            if i < num_q
                idx = find(mc_qx(:,:,1)>=edges(:,i,ii) ...
                    & mc_qx(:,:,1)<edges(:,i+1,ii));
            else
                idx = find(mc_qx(:,:,1)>=edges(:,i,ii) ...
                    & mc_qx(:,:,1)<=edges(:,i+1,ii));
            end
            A = mc_qx(:,:,2);
            qx2_vals = A(idx);
            m = mean(qx2_vals);
            m_qx(i,1) = m;
            sig = std(qx2_vals);
            serrX(i,1) = sig;
        end
    else
        for i = 1:num_q             % Calculate uncertainty for Y data
            if i < num_q
                idx = find(mc_qy(:,:,1)>=edges(:,i,ii) ...
                    & mc_qy(:,:,1)<edges(:,i+1,ii));
            else
                idx = find(mc_qy(:,:,1)>=edges(:,i,ii) ...
                    & mc_qy(:,:,1)<=edges(:,i+1,ii));
            end
            A = mc_qy(:,:,2);
            qy2_vals = A(idx);
            m = mean(qy2_vals);
            m_qy(i,1) = m;
            sig = std(qy2_vals);
            serrY(i,1) = sig;
        end
    end
end

%% Quantile and Residual Computation: User Data
% Computes quantiles and residuals for the user's data, and computes the
% difference between these and the residuals of each hypothetical climate
% scenario (i.e. ∆R). The algorithm then ID's signficiant quantiles at the
% specified degree of confidence.

multiWaitbar('Finishing Up...','Increment',0.01);

% ---- Convert Mg/Ca-IFA to Temperature ----

if ~isequal(X,T)
    switch eqn
        case 1
            X = (log(X/0.38)/0.09);
        case 2
            X = (log(X/0.34)/0.102);
        case 3
            X = ((log(X/0.38)/0.09) + (0.61*D) + 1.6);
        case 4
            X = (log(X/1.06)/0.048);
        case 5
            X = ((log(X/0.6)/0.08) + (2.8*D) + 5.4);
        case 6
            X = (log(X/0.41)/0.068);
        case 7
            X = (log(X/0.18)/0.12);
    end
end

switch eqn
    case 1
        Y = (log(Y/0.38)/0.09);
    case 2
        Y = (log(Y/0.34)/0.102);
    case 3
        Y = ((log(Y/0.38)/0.09) + (0.61*D) + 1.6);
    case 4
        Y = (log(Y/1.06)/0.048);
    case 5
        Y = ((log(Y/0.6)/0.08) + (2.8*D) + 5.4);
    case 6
        Y = (log(Y/0.41)/0.068);
    case 7
        Y = (log(Y/0.18)/0.12);
end

% ---- Compute Centered Quantiles ----

qx = quantile(X - mean(X),num_q)';
qy = quantile(Y - mean(Y),num_q)';

% ---- Calculate Residuals ----

yr_data = zeros(num_q,1);         % Vert distance b/w quantiles & 1:1 line
xr_data = zeros(num_q,1);         % Horiz distance b/w quantiles & 1:1 line

for i = 1:num_q
    yr = qy(i) - qx(i);           
    xr = qx(i) - qy(i);           
    yr_data(i,1) = yr;
    xr_data(i,1) = xr;
end

% ---- Calculate Residual Differences ----

d_diff = zeros(21,21,num_q);

for p = 1:21
    for w = 1:21
        for i = 1:num_q
            dif = abs(yr_alt(p,w,i) - yr_data(i));
            d_diff(p,w,i) = dif;
        end
    end
end

% ---- ID Significant Quantiles ----

lval = zeros(num_q,1);          % Assign a number to each quantile ID'ing
                                % its level of confidence: 0 = not sig;
                                % 1 = 80%; 2 = 90%; 3 = 95%; 4 = 99%
                                
for i = 1:num_q
    if abs(yr_data(i)) < serrY(i)*1.282 || abs(xr_data(i)) < serrX(i)*1.282
        lval(i,:) = 0;
        
    elseif abs(yr_data(i)) > serrY(i)*1.282 && abs(yr_data(i)) < serrY(i)*1.645
        if abs(xr_data(i)) > serrX(i)*1.282
            lval(i,:) = 1;
        else
            lval(i,:) = 0;
        end
    
    elseif abs(yr_data(i)) > serrY(i)*1.645 && abs(yr_data(i)) < serrY(i)*2
        if abs(xr_data(i)) > serrX(i)*1.645
            lval(i,:) = 2;        
        elseif abs(xr_data(i)) > serrX(i)*1.282
            lval(i,:) = 1;
        else
            lval(i,:) = 0;
        end
        
    elseif abs(yr_data(i)) > serrY(i)*2 && abs(yr_data(i)) < serrY(i)*3
        if abs(xr_data(i)) > serrX(i)*2
            lval(i,:) = 3;
        elseif abs(xr_data(i)) > serrX(i)*1.645
            lval(i,:) = 2;
        elseif abs(xr_data(i)) > serrX(i)*1.282
            lval(i,:) = 1;
        else
            lval(i,:) = 0;
        end
        
    elseif abs(yr_data(i)) > serrY(i)*3
        if abs(xr_data(i)) > serrX(i)*3
            lval(i,:) = 4;
        elseif abs(xr_data(i)) > serrX(i)*2
            lval(i,:) = 3;
        elseif abs(xr_data(i)) > serrX(i)*1.645
            lval(i,:) = 2;
        elseif abs(xr_data(i)) > serrX(i)*1.282
            lval(i,:) = 1;
        else
            lval(i,:) = 0;
        end
    end
end

%% Output Product #2: Q-Q Plot of User Data
% Generates a Q-Q plot of the inputted IFA data w/ uncertainty bounds. If
% there is no X population, the Q-Q plot compares the Y IFA data to the
% extracted temperature time series from ORA-S5.

multiWaitbar('Finishing Up...','Increment',0.01);

% ---- Generate Figure ----

figure(2)
clf;
hold on

% -- Plot Confidence Bounds --

if cl == 80
    envX = serrX*1.282;
    envY = serrY*1.282;
elseif cl == 90
    envX = serrX*1.645;
    envY = serrY*1.645;
elseif cl == 95
    envX = serrX*2;
    envY = serrY*2;
elseif cl == 99
    envX = serrX*3;
    envY = serrY*3;  
end

boundedline(qx,qy,envY,'transparency',0.15,'cmap',[255 51 51]./255);
hold on
boundedline(qx,qy,envX,'orientation','horiz','transparency',0.15,'cmap',[255 51 51]./255);
errorbar(qx,qy,envY,'color',[1 0.3 0.3],'capsize',0,'linewidth',1.5);
errorbar(qx,qy,envX,'horiz','color',[1 0.3 0.3],'capsize',0,'linewidth',1.5);

% -- Plot User Data --

if cl == 80
    for i = 1:length(lval)
        if lval(i) >= 1
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','w','linewidth',1.5);
            hold on
        else
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','k','linewidth',1.5);
            hold on
        end
    end
elseif cl == 90
    for i = 1:length(lval)
        if lval(i) >= 2
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','w','linewidth',1.5);
            hold on
        else
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','k','linewidth',1.5);
            hold on
        end
    end
elseif cl == 95
    for i = 1:length(lval)
        if lval(i) >= 3
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','w','linewidth',1.5);
            hold on
        else
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','k','linewidth',1.5);
            hold on
        end
    end
elseif cl == 99
    for i = 1:length(lval)
        if lval(i) == 4
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','w','linewidth',1.5);
            hold on
        else
            plot(qx(i),qy(i),'kd-','MarkerSize',6,'MarkerFaceColor','k','linewidth',1.5);
            hold on
        end
    end
end

rl = refline(1,0);
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca,'LineWidth',1.5,'FontName','Myraid Pro','fontsize',14,'Fontweight','demi');
xlabel('X Population Normalized Quantiles (°C)','FontName','Myraid Pro','Fontsize',14,'Fontweight','demi');
ylabel('Y Population Normalized Quantiles (°C)','FontName','Myraid Pro','Fontsize',14,'Fontweight','demi');
    
% -- Plot Inset Figure w/ Kernel Density Functions --

[fx, xi] = ksdensity(X);
[fy, yi] = ksdensity(Y);
ax2 = axes('Position',[.2 .6 .25 .25]);
plot(xi,fx,'color',rgb('black'),'linewidth',1.8);
hold (ax2, 'on'); %[255 0 0]
plot(yi,fy,'color',[255 51 51]./255,'linewidth',1.5);
xlabel(ax2,'Temperature (°C)','FontName','Myraid Pro');
ylabel(ax2,'Probability Density','FontName','Myraid Pro');
set(ax2,'linewidth',1.5,'box','off','FontName','Myraid Pro','fontweight','demi','fontsize',11);
if isequal(X,T)
    legend('ORA-S5 Reanalysis Data',['IFA Population ({\itn} = ',num2str(length(Y)),')'],'box','off','location','northoutside','FontName','Myraid Pro','fontsize',12);
else
    legend(['X Population ({\itn} = ',num2str(length(X)),')'],['Y Population ({\itn} = ',num2str(length(Y)),')'],'box','off','location','northoutside','FontName','Myraid Pro','fontsize',12);
end

% -- Save as PDF --

if savefig == 1
    f_filename = 'QUANTIFA_QQplot';
    export_fig(f_filename, '-pdf', '-append', '-transparent', '-r600');
end
 
% ---- Terminate Script ----

if NoIFA == 0.5                 % Terminates QUANTIFA if performing IFA
    toc                         % comparison w/ reanalysis data
    return
end 

%% Output Product #3: Data-Model Consistency Map
% Generates a heat map showing the proportion of those signficiant quantiles 
% that exhibit good data-model fit for each scenario.

multiWaitbar('Finishing Up...','Increment',0.01);

% ---- Generate Heat Maps ----

% -- Define Confidence Level --

if cl == 80
    j = 1;
elseif cl == 90
    j = 2;
elseif cl == 95
    j = 3;
elseif cl == 99
    j = 4;
end

idx = find(lval>=j);

if isempty(idx)
    if j == 1
        warning('No quantiles signficant w/ 80% confidence.');
        disp(' ');
    elseif j == 2
        warning('No quantiles signficant w/ 90% confidence.');
        disp(' ');
    elseif j == 3
        warning('No quantiles signficant w/ 95% confidence.');
        disp(' ');
    elseif j == 4
        warning('No quantiles signficant w/ 99% confidence.');
        disp(' ');
    end  
    return
end

% -- Test Data-Model Fit --

lval2 = zeros(21,21,length(idx));

for p = 1:21
    for w = 1:21
        for ii = 1:length(idx)
            idr = d_diff(p,w,idx(ii));
            if idr <= serrY(ii)
                lval2(p,w,ii) = 1;
            else
                lval2(p,w,ii) = 0;
            end
        end
    end
end
s = sum(lval2,3);
pfit = (s./length(idx))*100;
pfit = round(pfit);

% -- Compute Barycenter --

[bX, bY] = barycenter2D(1:21,1:21,pfit);
bX = 10*(bX - 11); bY = 10*(bY - 11);   % Translates to native coordinates

% -- Generate Figure --

h = figure;
axis tight manual;
clf;

% -- Plot Heat Map --

hh = imagesc(ww,pp,pfit);
hold on
plot(bX,bY,'ko','markersize',30,'linewidth',2);
plot(bX,bY,'kx','markersize',30,'linewidth',2);

% -- Create Color Map --

crange = [0 100];
t = linspace(crange(1),crange(2),10000);
cmap = zeros(10000,3);

mask = t > 90;
cmap(mask,:) = repmat([190 0 0]./255,[sum(mask), 1]);
mask = t > 80 & t <= 90;
cmap(mask,:) = repmat([230 0 0]./255,[sum(mask), 1]);
mask = t > 70 & t <= 80;
cmap(mask,:) = repmat([255 0 0]./255,[sum(mask), 1]);
mask = t > 60 & t <= 70;
cmap(mask,:) = repmat([255 50 50]./255,[sum(mask), 1]);
mask = t > 50 & t <= 60;
cmap(mask,:) = repmat([255 90 90]./255,[sum(mask), 1]);
mask = t > 40 & t <= 50;
cmap(mask,:) = repmat([255 130 130]./255,[sum(mask), 1]);
mask = t > 30 & t <= 40;
cmap(mask,:) = repmat([255 160 160]./255,[sum(mask), 1]);
mask = t > 20 & t <= 30;
cmap(mask,:) = repmat([255 190 190]./255,[sum(mask), 1]);
mask = t > 10 & t <= 20;
cmap(mask,:) = repmat([255 235 235]./255,[sum(mask), 1]);    
mask = t <= 10;
cmap(mask,:) = repmat([255 255 255]./255,[sum(mask), 1]);

colormap(cmap);
caxis([0 100]);

% -- Figure Properties --

xlabel('Change in ENSO Variability (%)');
ylabel('Change in Annual Cycle (%)');
set(gca,'Ydir','normal')

if j == 1
    title({['Quantiles Significant at 80% Confidence'] ['Number of Quantiles = ' num2str(length(idx))]});
elseif j == 2
    title({['Quantiles Significant at 90% Confidence'] ['Number of Quantiles = ' num2str(length(idx))]});
elseif j == 3
    title({['Quantiles Significant at 95% Confidence'] ['Number of Quantiles = ' num2str(length(idx))]});
elseif j == 4
    title({['Quantiles Significant at 99% Confidence'] ['Number of Quantiles = ' num2str(length(idx))]});
end

% -- Save as PDF --

if savefig == 1
    f_filename = 'QUANTIFA_ConsistencyMap';
    export_fig(f_filename, '-pdf', '-append', '-transparent', '-r600');
end

%% Supplemental Figures

multiWaitbar('Finishing Up...','Increment',0.01);

if suppfig == 1

    % -- Supp Figure 1: Modeled Modern and Select Altered Time Series --
    
    figure()
    clf;
    
    yrs = [1:tmslice]./12;
    
    subplot(7,1,1)
    plot(yrs,mod_series,'k','LineWidth',1.5);
    title('Modeled Modern Variability','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,2)
    a = squeeze(alt_series(21,21,:));
    plot(yrs,a,'color',[0.6 0.06 0.06],'LineWidth',1.5);
    title('2x Seasonality / 2x ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,3)
    a = squeeze(alt_series(11,21,:));
    plot(yrs,a,'color',[1 0 0],'LineWidth',1.5);
    title('Normal Seasonality / 2x ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,4)
    a = squeeze(alt_series(1,21,:));
    plot(yrs,a,'color',[1 0.5 0.5],'LineWidth',1.5);
    title('No Seasonality / 2x ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,5)
    a = squeeze(alt_series(21,11,:));
    plot(yrs,a,'color',[0.4 0.53 1],'LineWidth',1.5);
    title('2x Seasonality / Normal ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,6)
    a = squeeze(alt_series(21,1,:));
    plot(yrs,a,'color',[0 0.22 1],'LineWidth',1.5);
    title('2x Seasonality / No ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    
    subplot(7,1,7)
    a = squeeze(alt_series(1,1,:));
    plot(yrs,a,'color',[0 0.01 0.7],'LineWidth',1.5);
    title('No Seasonality / No ENSO "Intensity"','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
    xlabel('Months','fontsize',12,'fontweight','bold');
    set(gcf,'position',[50,50,1200,1200]);
    
    % -- Supp Figure 2: Q-Q Plots of Altered Time Series --
    
    figure()
    clf;
    
    color = flipud(redblue(21,[-1,1],'white'));
    
    for i = 1:21
        scatter(bins,squeeze(m_alt_qs(11,i,:)),'o','filled',...
                'MarkerFaceColor',color(i,:),'MarkerEdgeColor','k');
        hold on
    end
    title('Q-Q Expressions w/ Varying ENSO "Intensity" (Seasonality = Modern)');
    legend('-100%','-90%','-80%','-70%','-60%','-50%','-40%','-30%',...
           '-20%','-10%','0%','10%','20%','30%','40%','50%','60%',...
           '70%','80%','90%','100%','location','northwest','box','off');
    
    rl = refline(1,0);
    rl.Color = 'k';
    rl.LineWidth = 1.5;
    rl.HandleVisibility = 'off';
    set(gca,'LineWidth',1.5,'FontName','Myraid Pro','fontsize',14,...
        'Fontweight','demi');
    xlabel('Modern Population Mean Centered Quantiles (°C)','FontName',...
           'Myraid Pro','Fontsize',14,'Fontweight','demi');
    ylabel('Altered Population Mean Centered Quantiles (°C)','FontName',...
           'Myraid Pro','Fontsize',14,'Fontweight','demi');
    
    figure()
    clf;
    
    for i = 1:21
        scatter(bins,squeeze(m_alt_qs(i,11,:)),'o','filled',...
                'MarkerFaceColor',color(i,:),'MarkerEdgeColor','k');
        hold on
    end
    
    title('Q-Q Expressions w/ Varying Seasonality (ENSO Intensity = Modern)');
    legend('-100%','-90%','-80%','-70%','-60%','-50%','-40%','-30%',...
           '-20%','-10%','0%','10%','20%','30%','40%','50%','60%',...
           '70%','80%','90%','100%','location','northwest','box','off');
    rl = refline(1,0);
    rl.Color = 'k';
    rl.LineWidth = 1.5;
    rl.HandleVisibility = 'off';
    set(gca,'LineWidth',1.5,'FontName','Myraid Pro','fontsize',14,...
        'Fontweight','demi');
    xlabel('Modern Population Mean Normalized Quantiles (°C)','FontName',...
           'Myraid Pro','Fontsize',14,'Fontweight','demi');
    ylabel('Altered Population Mean Normalized Quantiles (°C)','FontName',...
           'Myraid Pro','Fontsize',14,'Fontweight','demi');

end
%% Tidy Workspace
% Clean up workspace and save.

% NOTE: If you would rather have all of the variables saved to the workspace
% (this can be useful for debugging, if I've messed up somewhere...) just 
% comment out the bloc below.

clearvars -except alt_qs alt_series anerr bdyp bins calerr cl ctp D d_diff ...
    dep depths edges envX envY eqn fpr_m fpr_sd fx fy lat lon lval M m_qx ...
    m_qy MgCa mod_qs mod_series nrst_depth num_p num_q numX numY ORAS5 pfit...
    qx qy serrX serrY savefig T tmslice wtp X xi xr_data Y yi yr_alt yr_data

save('QUANTIFA_Output.mat')

multiWaitbar('Finishing Up...','Close');

disp('Run Complete!');
toc

%% Built-In Functions --> Homemade
function [ilat,ilon,idepth] = ta_oraco(lat,lon,nrst_depth,depths)
% Converts lat/lon/dep input into corresponding ORA-S5 matrix indicies
% (FOR TROPICAL ATLANTIC DATASET)
%
% - created by R. Glaubke (10/21/2019)
%
%--------------------------------------------------------------------------
% Convert latitude from inputted coordinate to index that matches input
ilat = lat + 30.5;                 % e.g. lat(0.5) = index 31

% Convert longitude from inputted coordinate to index that matches input
if (lon < 0)                       % -ve coordinates converted to °E
    ilon = (360 + lon) - 259;      % e.g. lon(-90) = 270 - 260 = index 10
elseif (lon > 260 && lon <= 360)
    ilon = lon - 259;
elseif (lon > 0 && lon <= 15)
    ilon = lon + 101;              % e.g. lon(15) = index 116 
end

% Generate interpolated function between depth levels and indicies
F = griddedInterpolant(depths,1:1:75);
idepth = F(nrst_depth);            % Pull out index of nearest depth level
end

function [ilat,ilon,idepth] = ti_oraco(lat,lon,nrst_depth,depths)
% Converts lat/lon/dep input into corresponding ORA-S5 matrix indicies
% (FOR TROPICAL INDIAN DATASET)
%
% - created by R. Glaubke (10/21/2019)
%
%--------------------------------------------------------------------------
% Convert latitude from inputted coordinate to index that matches input
ilat = lat + 30.5;                 % e.g. lat(0.5) = index 31

% Convert longitude from inputted coordinate to index that matches input
if (lon < 0)                       % -ve coordinates converted to °E
    ilon = (360 + lon) - 120;      % e.g. lon(-90) = 270 - 120 = index 150
else
    ilon = lon - 120;              % e.g. lon(165) = index 45
end

% Generate interpolated function between depth levels and indicies
F = griddedInterpolant(depths,1:1:75);
idepth = F(nrst_depth);            % Pull out index of nearest depth level
end

function [rata,fw,powf,wp,rp] = powrat(ts)
% This function gives you the ratio of the power of two frequencies (annual
% versus interannul) in an MTM Spectra of a (monthly) time series.
%   powrat(<time_series>) should give you an output of the ratio of annual
%   to interannual frequencies, frequency window & power
%
% - created by K. Thirumalai (03/25/2013)
%   kau@ig.utexas.edu
% - updated by R. Glaubke (11/19/2019)
%   rglaubke@marine.rutgers.edu

%--------------------------------------------------------------------------
dt = 1/12;
ts = ts((all((~isnan(ts)),2)),:);       % Work w/ non-NaN values -RG
[Pmtm,~,w]=pmtm(ts,3);
fw = w/(2*pi)/dt;
powf = fw.*Pmtm;
t = 1./fw;

aa = find(t>0.9 & t<1.11);
annu = sum(powf(aa(1):aa(end)));

ia = find(t>2 & t<7);
interann = sum(powf(ia(1):ia(end)));

wpf = find(t<0.4 & t>0.2);
wp = sum(powf(wpf(1):wpf(end)));

rpf = find(t>10 & t<35);
rp = sum(powf(rpf(1):rpf(end)));

rata = annu/interann;
end

%% Built-In Functions --> Outsourced (see descriptions for references)
function [cancel,figh] = multiWaitbar( label, varargin )
% multiWaitbar: add, remove or update an entry on the multi waitbar
%
% Citation:
% Ben Tordoff (2023). multiWaitbar (https://www.mathworks.com/matlabcentral/
% fileexchange/26589-multiwaitbar), MATLAB Central File Exchange. Retrieved 
% March 2, 2023.
%
%   multiWaitbar(LABEL,VALUE) adds a waitbar for the specified label, or
%   if it already exists updates the value. LABEL must be a string and
%   VALUE a number between zero and one or the string 'Close' to remove the
%   entry Setting value equal to 0 or 'Reset' will cause the progress bar
%   to reset and the time estimate to be re-initialized.
%
%   multiWaitbar(LABEL,COMMAND,VALUE,...)  or
%   multiWaitbar(LABEL,VALUE,COMMAND,VALUE,...)
%   passes one or more command/value pairs for changing the named waitbar
%   entry. Possible commands include:
%   'Value'       Set the value of the named waitbar entry. The
%                 corresponding value must be a number between 0 and 1.
%   'Increment'   Increment the value of the named waitbar entry. The
%                 corresponding value must be a number between 0 and 1.
%   'Color'       Change the color of the named waitbar entry. The
%                 value must be an RGB triple, e.g. [0.1 0.2 0.3], or a
%                 single-character color name, e.g. 'r', 'b', 'm'.
%   'Relabel'     Change the label of the named waitbar entry. The
%                 value must be the new name.
%   'Reset'       Set the named waitbar entry back to zero and reset its
%                 timer. No value need be specified.
%   'CanCancel'   [on|off] should a "cancel" button be shown for this bar
%                 (default 'off').
%   'CancelFcn'   Function to call in the event that the user cancels.
%   'ResetCancel' Reset the "cancelled" flag for an entry (ie. if you
%                 decide not to cancel).
%   'Close'       Remove the named waitbar entry.
%   'Busy'        Puts this waitbar in "busy mode" where a small bar
%                 bounces back and forth. Return to normal progress display
%                 using the 'Reset' command.
%
%   cancel = multiWaitbar(LABEL,VALUE) also returns whether the user has
%   clicked the "cancel" button for this entry (true or false). Two
%   mechanisms are provided for cancelling an entry if the 'CanCancel'
%   setting is 'on'. The first is just to check the return argument and if
%   it is true abort the task. The second is to set a 'CancelFcn' that is
%   called when the user clicks the cancel button, much as is done for
%   MATLAB's built-in WAITBAR. In either case, you can use the
%   'ResetCancel' command if you don't want to cancel after all. 
%
%   Author: Ben Tordoff
%   Copyright 2007-2020 The MathWorks, Inc.

persistent FIGH;
cancel = false;
% Check basic inputs
error( nargchk( 1, inf, nargin ) ); %#ok<NCHKN> - kept for backwards compatibility
if ~ischar( label )
    error( 'multiWaitbar:BadArg', 'LABEL must be the name of the progress entry (i.e. a string)' );
end
% Try to get hold of the figure
if isempty( FIGH ) || ~ishandle( FIGH )
    FIGH = findall( 0, 'Type', 'figure', 'Tag', 'multiWaitbar:Figure' );
    if isempty(FIGH)
        FIGH = iCreateFig();
    else
        FIGH = handle( FIGH(1) );
    end
end
% Check for close all and stop early
if any( strcmpi( label, {'CLOSEALL','CLOSE ALL'} ) )
    iDeleteFigure(FIGH);
    figh = [];
    return;
end
if nargout>1
    figh = FIGH;
end
% Make sure we're on-screen
if ~strcmpi( FIGH.Visible, 'on' )
    FIGH.Visible = 'on';
end
% Make sure the timer is still valid - it can be found and deleted
% externally.
if ~isvalid( getappdata( FIGH, 'BusyTimer' ) )
    setappdata( FIGH, 'BusyTimer', iCreateTimer(FIGH) );
end
% Get the list of entries and see if this one already exists
entries = getappdata( FIGH, 'ProgressEntries' );
if isempty(entries)
    idx = [];
else
    idx = find( strcmp( label, {entries.Label} ), 1, 'first' );
end
bgcol = getappdata( FIGH, 'DefaultProgressBarBackgroundColor' );
% If it doesn't exist, create it
needs_redraw = false;
entry_added = isempty(idx);
if entry_added
    % Create a new entry
    defbarcolor = getappdata( FIGH, 'DefaultProgressBarColor' );
    entries = iAddEntry( FIGH, entries, label, 0, defbarcolor, bgcol );
    idx = numel( entries );
end
% Check if the user requested a cancel
if nargout
    cancel = entries(idx).Cancel;
end
% Parse the inputs. We shortcut the most common case as an efficiency
force_update = false;
if nargin==2 && isnumeric( varargin{1} )
    entries(idx).LastValue = entries(idx).Value;
    entries(idx).Value = max( 0, min( 1, varargin{1} ) );
    entries(idx).Busy = false;
    needs_update = true;
else
    [params,values] = iParseInputs( varargin{:} );
    
    needs_update = false;
    for ii=1:numel( params )
        switch upper( params{ii} )
            case 'BUSY'
                entries(idx).Busy = true;
                needs_update = true;
                
            case 'VALUE'
                entries(idx).LastValue = entries(idx).Value;
                entries(idx).Value = max( 0, min( 1, values{ii} ) );
                entries(idx).Busy = false;
                needs_update = true;
                
            case {'INC','INCREMENT'}
                entries(idx).LastValue = entries(idx).Value;
                entries(idx).Value = max( 0, min( 1, entries(idx).Value + values{ii} ) );
                entries(idx).Busy = false;
                needs_update = true;
                
            case {'COLOR','COLOUR'}
                entries(idx).CData = iMakeColors( values{ii}, 16 );
                needs_update = true;
                force_update = true;
                
            case {'RELABEL', 'UPDATELABEL'}
                % Make sure we have a string as the value and that it
                % doesn't already appear
                if ~ischar( values{ii} )
                    error( 'multiWaitbar:BadString', 'Value for ''Relabel'' must be a string.' );
                end
                if ismember( values{ii}, {entries.Label} )
                    error( 'multiWaitbar:NameAlreadyExists', 'Cannot relabel an entry to a label that already exists.' );
                end
                entries(idx).Label = values{ii};
                if ~isempty(entries(idx).CancelButton)
                    set( entries(idx).CancelButton, 'Callback', @(src,evt) iCancelEntry(src, values{ii}) );
                end
                needs_update = true;
                force_update = true;
                
            case {'CANCANCEL'}
                if ~ischar( values{ii} ) || ~any( strcmpi( values{ii}, {'on','off'} ) )
                    error( 'multiWaitbar:BadString', 'Parameter ''CanCancel'' must be a ''on'' or ''off''.' );
                end
                entries(idx).CanCancel = strcmpi( values{ii}, 'on' );
                entries(idx).Cancel = false;
                needs_redraw = true;
                
            case {'RESETCANCEL'}
                entries(idx).Cancel = false;
                needs_redraw = true;
                
            case {'CANCELFCN'}
                if ~isa( values{ii}, 'function_handle' )
                    error( 'multiWaitbar:BadFunction', 'Parameter ''CancelFcn'' must be a valid function handle.' );
                end
                entries(idx).CancelFcn = values{ii};
                if ~entries(idx).CanCancel
                    entries(idx).CanCancel = true;
                end
                needs_redraw = true;
                
            case {'CLOSE','DONE'}
                if ~isempty(idx)
                    % Remove the selected entry
                    entries = iDeleteEntry( entries, idx );
                end
                if isempty( entries )
                    iDeleteFigure( FIGH );
                    % With the window closed, there's nothing else to do
                    return;
                else
                    needs_redraw = true;
                end
                % We can't continue after clearing the entry, so jump out
                break;
                
            otherwise
                error( 'multiWaitbar:BadArg', 'Unrecognized command: ''%s''', params{ii} );
                
        end
    end
end
% Now work out what to update/redraw
if needs_redraw
    setappdata( FIGH, 'ProgressEntries', entries );
    iRedraw( FIGH );
    % NB: Redraw includes updating all bars, so never need to do both
elseif needs_update
    [entries(idx),needs_redraw] = iUpdateEntry( entries(idx), force_update );
    setappdata( FIGH, 'ProgressEntries', entries );
    % NB: if anything was updated onscreen, "needs_redraw" is now true.
end
if entry_added || needs_redraw
    % If the shape or size has changed, do a full redraw, including events
    drawnow();
end
% If we have any "busy" entries, start the timer, otherwise stop it.
myTimer = getappdata( FIGH, 'BusyTimer' );
if any([entries.Busy])
    if strcmpi(myTimer.Running,'off')
        start(myTimer);
    end
else
    if strcmpi(myTimer.Running,'on')
        stop(myTimer);
    end
end
end % multiWaitbar
%-------------------------------------------------------------------------%
function [params, values] = iParseInputs( varargin )
% Parse the input arguments, extracting a list of commands and values
idx = 1;
params = {};
values = {};
if nargin==0
    return;
end
if isnumeric( varargin{1} )
    params{idx} = 'Value';
    values{idx} = varargin{1};
    idx = idx + 1;
end
while idx <= nargin
    param = varargin{idx};
    if ~ischar( param )
        error( 'multiWaitbar:BadSyntax', 'Additional properties must be supplied as property-value pairs' );
    end
    params{end+1,1} = param; %#ok<AGROW>
    values{end+1,1} = []; %#ok<AGROW>
    switch upper( param )
        case {'DONE','CLOSE','RESETCANCEL'}
            % No value needed, and stop
            break;
        case {'BUSY'}
            % No value needed but parsing should continue
            idx = idx + 1;
        case {'RESET','ZERO','SHOW'}
            % All equivalent to saying ('Value', 0)
            params{end} = 'Value';
            values{end} = 0;
            idx = idx + 1;
        otherwise
            if idx==nargin
                error( 'multiWaitbar:BadSyntax', 'Additional properties must be supplied as property-value pairs' );
            end
            values{end,1} = varargin{idx+1};
            idx = idx + 2;
    end
end
if isempty( params )
    error( 'multiWaitbar:BadSyntax', 'Must specify a value or a command' );
end
end % iParseInputs
%-------------------------------------------------------------------------%
function fobj = iCreateFig()
% Create the progress bar group window
bgcol = get(0,'DefaultUIControlBackgroundColor');
f = figure( ...
    'Name', 'QUANTIFA', ...
    'Tag', 'multiWaitbar:Figure', ...
    'Color', bgcol, ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'WindowStyle', 'normal', ... % We don't want to be docked!
    'HandleVisibility', 'off', ...
    'IntegerHandle', 'off', ...
    'Visible', 'off', ...
    'NumberTitle', 'off' );
% Resize and centre on the first screen
screenSize = get(0,'ScreenSize');
figSz = [360 42];
figPos = ceil((screenSize(1,3:4)-figSz)/2);
fobj = handle( f );
fobj.Position = [figPos, figSz];
setappdata( fobj, 'ProgressEntries', [] );
% Make sure we have the image
defbarcolor = [0.8 0.0 0.1];
barbgcol = uint8( 255*0.75*bgcol );
setappdata( fobj, 'DefaultProgressBarBackgroundColor', barbgcol );
setappdata( fobj, 'DefaultProgressBarColor', defbarcolor );
setappdata( fobj, 'DefaultProgressBarSize', [350 16] );
setappdata( fobj, 'MaxEntryRows', 10 );
% Create the timer to use for "Busy" mode, being sure to delete any
% existing ones
delete( timerfind('Tag', 'MultiWaitbarTimer') );
myTimer = iCreateTimer(f);
setappdata( fobj, 'BusyTimer', myTimer );
% Setup the resize function after we've finished setting up the figure to
% avoid excessive redraws
fobj.ResizeFcn = @iRedraw;
fobj.CloseRequestFcn = @iCloseFigure;
end % iCreateFig
%-------------------------------------------------------------------------%
function t = iCreateTimer(fig)
t = timer( ...
    'TimerFcn', @(src,evt) iTimerFcn(fig), ...
    'Period', 0.02, ...
    'ExecutionMode', 'FixedRate', ...
    'Tag', 'MultiWaitbarTimer' );
end
%-------------------------------------------------------------------------%
function cdata = iMakeColors( baseColor, height )
% Creates a shiny bar from a single base color
lightColor = [1 1 1];
badColorErrorID = 'multiWaitbar:BadColor';
badColorErrorMsg = 'Colors must be a three element vector [R G B] or a single character (''r'', ''g'' etc.)';
if ischar(baseColor)
    switch upper(baseColor)
        case 'K'
            baseColor = [0.1 0.1 0.1];
        case 'R'
            baseColor = [0.8 0 0];
        case 'G'
            baseColor = [0 0.6 0];
        case 'B'
            baseColor = [0 0 0.8];
        case 'C'
            baseColor = [0.2 0.8 0.9];
        case 'M'
            baseColor = [0.6 0 0.6];
        case 'Y'
            baseColor = [0.9 0.8 0.2];
        case 'W'
            baseColor = [0.9 0.9 0.9];
        otherwise
            error( badColorErrorID, badColorErrorMsg );
    end
else
    if numel(baseColor) ~= 3
        error( badColorErrorID, badColorErrorMsg );
    end
    if isa( baseColor, 'uint8' )
        baseColor = double( baseColor ) / 255;
    elseif isa( baseColor, 'double' )
        if any(baseColor>1) || any(baseColor<0)
            error( 'multiWaitbar:BadColorValue', 'Color values must be in the range 0 to 1 inclusive.' );
        end
    else
        error( badColorErrorID, badColorErrorMsg );
    end
end
% By this point we should have a double precision 3-element vector.
cols = repmat( baseColor, [height, 1] );
breaks = max( 1, round( height * [1 25 50 75 88 100] / 100 ) );
cols(breaks(1),:) = 0.6*baseColor;
cols(breaks(2),:) = lightColor - 0.4*(lightColor-baseColor);
cols(breaks(3),:) = baseColor;
cols(breaks(4),:) = min( baseColor*1.2, 1.0 );
cols(breaks(5),:) = min( baseColor*1.4, 0.95 ) + 0.05;
cols(breaks(6),:) = min( baseColor*1.6, 0.9 ) + 0.1;
y = 1:height;
cols(:,1) = max( 0, min( 1, interp1( breaks, cols(breaks,1), y, 'pchip' ) ) );
cols(:,2) = max( 0, min( 1, interp1( breaks, cols(breaks,2), y, 'pchip' ) ) );
cols(:,3) = max( 0, min( 1, interp1( breaks, cols(breaks,3), y, 'pchip' ) ) );
cdata = uint8( 255 * cat( 3, cols(:,1), cols(:,2), cols(:,3) ) );
end % iMakeColors
%-------------------------------------------------------------------------%
function cdata = iMakeBackground( baseColor, height )
% Creates a shaded background
if isa( baseColor, 'uint8' )
    baseColor = double( baseColor ) / 255;
end
ratio = 1 - exp( -0.5-2*(1:height)/height )';
cdata = uint8( 255 * cat( 3, baseColor(1)*ratio, baseColor(2)*ratio, baseColor(3)*ratio ) );
end % iMakeBackground
%-------------------------------------------------------------------------%
function entries = iAddEntry( parent, entries, label, value, color, bgcolor )
% Add a new entry to the progress bar
% Create bar coloring
psize = getappdata( parent, 'DefaultProgressBarSize' );
cdata = iMakeColors( color, 16 );
% Create background image
barcdata = iMakeBackground( bgcolor, psize(2) );
% Work out the size in advance
mypanel = uipanel( 'Parent', parent, 'Units', 'Pixels', 'BorderType', 'beveledout' );
labeltext = uicontrol( 'Style', 'Text', ...
    'String', label, ...
    'Parent', parent, ...
    'HorizontalAlignment', 'Left' );
etatext = uicontrol( 'Style', 'Text', ...
    'String', '', ...
    'Parent', parent, ...
    'HorizontalAlignment', 'Right' );
progresswidget = uicontrol( 'Style', 'Checkbox', ...
    'String', '', ...
    'Parent', parent, ...
    'Position', [5 5 psize], ...
    'CData', barcdata );
cancelwidget = uicontrol( 'Style', 'PushButton', ...
    'String', '', ...
    'FontWeight', 'Bold', ...
    'Parent', parent, ...
    'Position', [5 5 16 16], ...
    'CData', iMakeCross( 8 ), ...
    'Callback', @(src,evt) iCancelEntry( src, label ), ...
    'Visible', 'off' );
newentry = struct( ...
    'Label', label, ...
    'Value', value, ...
    'LastValue', inf, ...
    'Created', tic(), ...
    'LabelText', labeltext, ...
    'ETAText', etatext, ...
    'ETAString', '', ...
    'Progress', progresswidget, ...
    'ProgressSize', psize, ...
    'Panel', mypanel, ...
    'BarCData', barcdata, ...
    'CData', cdata, ...
    'BackgroundCData', barcdata, ...
    'CanCancel', false, ...
    'CancelFcn', [], ...
    'CancelButton', cancelwidget, ...
    'Cancel', false, ...
    'Busy', false );
if isempty( entries )
    entries = newentry;
else
    entries = [entries;newentry];
end
% Store in figure before the redraw
setappdata( parent, 'ProgressEntries', entries );
if strcmpi( get( parent, 'Visible' ), 'on' )
    iRedraw( parent, [] );
else
    set( parent, 'Visible', 'on' );
end
end % iAddEntry
%-------------------------------------------------------------------------%
function entries = iDeleteEntry( entries, idx )
delete( entries(idx).LabelText );
delete( entries(idx).ETAText );
delete( entries(idx).CancelButton );
delete( entries(idx).Progress );
delete( entries(idx).Panel );
entries(idx,:) = [];
end % iDeleteEntry
%-------------------------------------------------------------------------%
function entries = iCancelEntry( src, name )
figh = ancestor( src, 'figure' );
entries = getappdata( figh, 'ProgressEntries' );
if isempty(entries)
    % The entries have been lost - nothing can be done.
    return
end
idx = find( strcmp( name, {entries.Label} ), 1, 'first' );
% Set the cancel flag so that the user is told on next update
entries(idx).Cancel = true;
setappdata( figh, 'ProgressEntries', entries );
% If a user function is supplied, call it
if ~isempty( entries(idx).CancelFcn )
    feval( entries(idx).CancelFcn, name, 'Cancelled' );
end
end % iCancelEntry
%-------------------------------------------------------------------------%
function [entry,updated] = iUpdateEntry( entry, force )
% Update one progress bar
% Deal with busy entries separately
if entry.Busy
    entry = iUpdateBusyEntry(entry);
    updated = true;
    return;
end
% Some constants
marker_weight = 0.8;
% Check if the label needs updating
updated = force;
val = entry.Value;
lastval = entry.LastValue;
% Now update the bar
psize = entry.ProgressSize;
filled = max( 1, round( val*psize(1) ) );
lastfilled = max( 1, round( lastval*psize(1) ) );
% We do some careful checking so that we only redraw what we have to. This
% makes a small speed difference, but every little helps!
if force || (filled<lastfilled)
    % Create the bar background
    startIdx = 1;
    bgim = entry.BackgroundCData(:,ones( 1, ceil(psize(1)-filled) ),:);
    barim = iMakeBarImage(entry.CData, startIdx, filled);
    progresscdata = [barim,bgim];
    
    % Add light/shadow around the markers
    markers = round( (0.1:0.1:val)*psize(1) );
    markers(markers<=startIdx | markers>(filled-2)) = [];
    highlight = [marker_weight*entry.CData, 255 - marker_weight*(255-entry.CData)];
    for ii=1:numel( markers )
        progresscdata(:,markers(ii)+[-1,0],:) = highlight;
    end
    
    % Set the image into the checkbox
    entry.BarCData = progresscdata;
    set( entry.Progress, 'cdata', progresscdata );
    updated = true;
    
elseif filled > lastfilled
    % Just need to update the existing data
    progresscdata = entry.BarCData;
    startIdx = max(1,lastfilled-1);
    % Repmat is the obvious way to fill the bar, but BSXFUN is often
    % faster. Indexing is obscure but faster still.
    progresscdata(:,startIdx:filled,:) = iMakeBarImage(entry.CData, startIdx, filled);
    
    % Add light/shadow around the markers
    markers = round( (0.1:0.1:val)*psize(1) );
    markers(markers<startIdx | markers>(filled-2)) = [];
    highlight = [marker_weight*entry.CData, 255 - marker_weight*(255-entry.CData)];
    for ii=1:numel( markers )
        progresscdata(:,markers(ii)+[-1,0],:) = highlight;
    end
    
    entry.BarCData = progresscdata;
    set( entry.Progress, 'CData', progresscdata );
    updated = true;
end
% As an optimization, don't update any text if the bar didn't move and the
% percentage hasn't changed
decval = round( val*100 );
lastdecval = round( lastval*100 );
if ~updated && (decval == lastdecval)
    return
end
% Now work out the remaining time
minTime = 3; % secs
if val <= 0
    % Zero value, so clear the eta
    entry.Created = tic();
    elapsedtime = 0;
    etaString = '';
else
    elapsedtime = round(toc( entry.Created )); % in seconds
    
    % Only show the remaining time if we've had time to estimate
    if elapsedtime < minTime
        % Not enough time has passed since starting, so leave blank
        etaString = '';
    else
        % Calculate a rough ETA
        eta = elapsedtime * (1-val) / val;
        etaString = iGetTimeString( eta );
    end
end
if ~isequal( etaString, entry.ETAString )
    set( entry.ETAText, 'String', etaString );
    entry.ETAString = etaString;
    updated = true;
end
% Update the label too
if force || elapsedtime > minTime
    if force || (decval ~= lastdecval)
        labelstr = [entry.Label, sprintf( ' (%d%%)', decval )];
        set( entry.LabelText, 'String', labelstr );
        updated = true;
    end
end
end % iUpdateEntry
function eta = iGetTimeString( remainingtime )
if remainingtime > 172800 % 2 days
    eta = sprintf( '%d days', round(remainingtime/86400) );
else
    if remainingtime > 7200 % 2 hours
        eta = sprintf( '%d hours', round(remainingtime/3600) );
    else
        if remainingtime > 120 % 2 mins
            eta = sprintf( '%d mins', round(remainingtime/60) );
        else
            % Seconds
            remainingtime = round( remainingtime );
            if remainingtime > 1
                eta = sprintf( '%d secs', remainingtime );
            elseif remainingtime == 1
                eta = '1 sec';
            else
                eta = ''; % Nearly done (<1sec)
            end
        end
    end
end
end % iGetTimeString
%-------------------------------------------------------------------------%
function entry = iUpdateBusyEntry( entry )
% Update a "busy" progress bar
% Make sure the widget is still OK
if ~ishandle(entry.Progress)
    return
end
% Work out the new position. Since the bar is 0.1 long and needs to bounce,
% the position varies from 0 up to 0.9 then back down again. We achieve
% this with judicious use of "mod" with 1.8.
entry.Value = mod(entry.Value+0.01,1.8);
val = entry.Value;
if val>0.9
    % Moving backwards
    val = 1.8-val;
end
psize = entry.ProgressSize;
startIdx = max( 1, round( val*psize(1) ) );
endIdx = max( 1, round( (val+0.1)*psize(1) ) );
barLength = endIdx - startIdx + 1;
% Create the image
bgim = entry.BackgroundCData(:,ones( 1, psize(1) ),:);
barim = iMakeBarImage(entry.CData, 1, barLength);
bgim(:,startIdx:endIdx,:) = barim;
% Put it into the widget
entry.BarCData = bgim;
set( entry.Progress, 'CData', bgim );
end % iUpdateBusyEntry
%-------------------------------------------------------------------------%
function barim = iMakeBarImage(strip, startIdx, endIdx)
shadow1_weight = 0.4;
shadow2_weight = 0.7;
barLength = endIdx - startIdx + 1;
% Repmat is the obvious way to fill the bar, but BSXFUN is often
% faster. Indexing is obscure but faster still.
barim = strip(:,ones(1, barLength),:);
% Add highlight to the start of the bar
if startIdx <= 2 && barLength>=2
    barim(:,1,:) = 255 - shadow1_weight*(255-strip);
    barim(:,2,:) = 255 - shadow2_weight*(255-strip);
end
% Add shadow to the end of the bar
if endIdx>=4 && barLength>=2
    barim(:,end,:) = shadow1_weight*strip;
    barim(:,end-1,:) = shadow2_weight*strip;
end
end % iMakeBarImage
%-------------------------------------------------------------------------%
function iCloseFigure( fig, evt ) %#ok<INUSD>
% Closing the figure just makes it invisible
set( fig, 'Visible', 'off' );
end % iCloseFigure
%-------------------------------------------------------------------------%
function iDeleteFigure( fig )
% Actually destroy the figure
busyTimer = getappdata( fig, 'BusyTimer' );
if isvalid( busyTimer )
    stop( busyTimer );
end
delete( busyTimer );
delete( fig );
end % iDeleteFigure
%-------------------------------------------------------------------------%
function iRedraw( fig, evt ) %#ok<INUSD>
entries = getappdata( fig, 'ProgressEntries' );
fobj = handle( fig );
p = fobj.Position;
% p = get( fig, 'Position' );
border = 5;
textheight = 16;
barheight = 16;
panelheight = 10;
maxRows = getappdata( fig, 'MaxEntryRows' );
N = max( 1, numel( entries ) );
Nrows = min(maxRows, N);
Ncols = ceil(N ./ Nrows);
% Check the height is correct
heightperentry = textheight+barheight+panelheight;
requiredheight = 2*border + Nrows*heightperentry - panelheight;
if ~isequal( p(4), requiredheight )
    p(2) = p(2) + p(4) - requiredheight;
    p(4) = requiredheight;
    % In theory setting the position should re-trigger this callback, but
    % in practice it doesn't, probably because we aren't calling "drawnow".
    set( fig, 'Position', p )
end
width = floor((p(3) - Ncols*2*border) ./ Ncols);
setappdata( fig, 'DefaultProgressBarSize', [width barheight] );
for ii=1:numel( entries )
    col = ceil(ii./Nrows)-1;
    row = ii - col.*Nrows - 1;
    xpos = border + (width+2*border).*col;
    ypos = p(4) - border - heightperentry.*row;
    
    set( entries(ii).Panel, 'Position', [xpos-border ypos+panelheight/2-heightperentry width+2*border heightperentry] );
    set( entries(ii).LabelText, 'Position', [xpos ypos-textheight width*0.75 textheight] );
    set( entries(ii).ETAText, 'Position', [xpos+width*0.75 ypos-textheight width*0.25 textheight] );
    ypos = ypos - textheight;
    
    if entries(ii).CanCancel
        set( entries(ii).Progress, 'Position', [xpos ypos-barheight max(1,width-barheight+1) barheight] );
        entries(ii).ProgressSize = max(1,[width-barheight barheight]);
        set( entries(ii).CancelButton, 'Visible', 'on', 'Position', [p(3)-border-barheight ypos-barheight barheight barheight] );
    else
        set( entries(ii).Progress, 'Position', [xpos ypos-barheight width+1 barheight] );
        entries(ii).ProgressSize = [width barheight];
        set( entries(ii).CancelButton, 'Visible', 'off' );
    end
    entries(ii) = iUpdateEntry( entries(ii), true );
end
setappdata( fig, 'ProgressEntries', entries );
end % iRedraw
function cdata = iMakeCross( sz )
% Create a cross-shape icon of size sz*sz*3
cdata = diag(ones(sz,1),0) + diag(ones(sz-1,1),1) + diag(ones(sz-1,1),-1);
cdata = cdata + flip(cdata,2);
% Convert zeros to nans (transparent) and non-zeros to zero (black)
cdata(cdata == 0) = nan;
cdata(~isnan(cdata)) = 0;
% Convert to RGB
cdata = cat( 3, cdata, cdata, cdata );
end % iMakeCross
function iTimerFcn(fig)
% Timer callback for updating stuff every so often
entries = getappdata( fig, 'ProgressEntries' );
for ii=1:numel(entries)
    if entries(ii).Busy
        entries(ii) = iUpdateBusyEntry(entries(ii));
    end
end
setappdata( fig, 'ProgressEntries', entries );
end

function [barycenterX, barycenterY] = barycenter2D(x,y,C)
% BARYCENTER2D returns x and y barycenter values of density map C
%  Assumption: Input coordinated x and y are evenly spaced!
%
%  The same input as for imagesc(x,y,C) can be used.
%
%   INPUT (from Matlab(R))
%     x    ... vector of length n, x-coordinates of image
%     y    ... vector of length m, y-coordinates of image
%     C    ... matrix of size n x m, density value for each coordiante
%
%   OUTPUT
%     barycenterX  ... double, x-coordinate value of the barycenter
%     barycenterY  ... double, y-coordinate value of the barycenter
%
%   CREDIT
%   [1] Andrei Bobrov https://www.mathworks.com/matlabcentral/answers/
%                     363181-center-of-mass-and-total-mass-of-a-matrix
%   [2] Matt J(faster)https://www.mathworks.com/matlabcentral/answers/
%                     133395-compute-centroid-of-a-matrix
%
%  TEST Examples (uncomment the parts below and at the end and run cell)
%{
T{1} = [...
    4 1 1 1 1;...
    1 2 2 2 1;...
    1 2 3 7 1;...
    1 9 4 2 1;...
    1 2 5 2 1;...
    1 2 4 9 1;...
    1 7 3 2 1;...
    1 2 2 2 1;...
    1 1 1 1 4;...
    ];
T{2} = [...
    0 0 0 0 0;...
    1 0 0 0 0;...
    1 1 0 0 0;...
    1 1 1 0 0;...
    1 1 1 1 0;...
    ];
T{3} = [...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    1 1 0 0 0;...
    1 1 0 0 0;...
    ];
T{4} = [...
    5 0 0 0 3;...
    0 1 1 1 0;...
    0 1 0 1 0;...
    0 1 0 1 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 1 0 0;...
    0 0 1 0 0;...
    0 1 1 1 0;...
    3 0 0 0 5;...
    ];
T{5} = ones([80,140]);
T{6} = zeros(10);
C = T{1}; % select an example
x = 1:size(C,2);
y = 1:size(C,1);
x = x + 123; % shift x vector to show correct x transformation
y = y + 456; % shift y vector to show correct y transformation
%}
% (c) Assembled and tranformation by H. Penasso Feb/15/2020
%

% Method 1:
%{
    tic
    C = C/sum(C(:));
    [m,n]=size(C);
    [I,J]=ndgrid(1:m,1:n);
    barycenterIdxX = dot(J(:),C(:));
    barycenterIdxY = dot(I(:),C(:));
    toc
%}
% Method 2:
% {
    tot_mass = sum(C(:));
    [ii,jj] = ndgrid(1:size(C,1),1:size(C,2));
    barycenterIdxX = sum(jj(:).*C(:))/tot_mass;
    barycenterIdxY = sum(ii(:).*C(:))/tot_mass;
%}
% Align barycenter indices and coordinate values
% Index vectors corresponding to the previous results
    xIdx = 1:size(C,2); % x axis index vector
    yIdx = 1:size(C,1); % y axis index vector
% X: linear fit between x-index vector and x coordinate vector
    slopeX      = sum((xIdx-mean(xIdx)).*(x-mean(x)))./...
                  sum((xIdx-mean(xIdx)).^2); % slope
    interceptX  = mean(x)-slopeX*mean(xIdx); % intercept
  % evaluate barycenter x value using slope and intercept of fit in x:
    barycenterX = slopeX*barycenterIdxX+interceptX;
% Y: linear fit between y-index vector and y coordinate vector
    slopeY      = sum((yIdx-mean(yIdx)).*(y-mean(y)))./...
                  sum((yIdx-mean(yIdx)).^2); % slope
    interceptY  = mean(y)-slopeY*mean(yIdx); % intercept
  % evaluate barycenter y value using slope and intercept of fit in y:
    barycenterY = slopeY*barycenterIdxY+interceptY;
%{
% Plot for test example
    imagesc(x,y,C)
    hold on
    plot(barycenterX,barycenterY,'wx',...
         barycenterX,barycenterY,'wo','MarkerSize',20,'LineWidth',2)
%}
end

function y = redblue(varargin)
% generate a RED-BLUE colormap with zero as white.
% Syntax: y = redblue(n,clim,'black')
%   Typical usage: colormap(redblue(64))
%   Positive values are displayed as blue intensities and negative
%   values displayed as red intensities. Zero is white.
%   The clim values are used to find zero to set that to white.
% Arguments:
%   All arguments are optional and can be in any order.
%   n - number of color levels 
%       (optional with default as # of colors of current colormap)
%   clim - two element vector specifying the color limits
%       (optional with default as current axis color limits)
%   black - string ('k' or 'black') specifying zero as black.
%       (optional with default as zero is white)
% See: colormap
%
% Notes:
%   This creates a custom colormap for any image so that the value 
%   zero is either white or black. The colorbar scale will be skewed
%   toward red or blue depending on the caxis values of the image.
%
%   This version flattens the edges of the colormap to improve the 
%   visualization of the gradient. This is better for larger values of n.
%
%   Keep in mind that if the scale is very skewed, there will not
%   be much of a color gradient. The gradient can always be increased
%   by using your own clim values.
%   Example:
%   y = caxis;                      % e.g. y = [-11,-5] 
%   colormap(redblue(64)            % and not much gradient
%   colormap(redblue(64,[-11,0]))   % white is -5 with larger gradient
% created   3/24/2020   Mirko Hrovat    mihrovat@outlook.com
% modified  5/04/2020   Changed algorithm for colormap calculation, eliminating
%   an error which occurs for small "clim" differences.
ntotmax = 1e4;
n = [];     clim = [];      black = false;
for k = 1:nargin
    a = varargin{k};
    switch true
        case ischar(a)
            switch true
                case strcmpi(a,'k')
                    black = true;
                case strcmpi(a,'black')
                    black = true;
            end
        case isnumeric(a) && numel(a) == 1
            n = a;
        case isnumeric(a) && numel(a) == 2
            clim = a;
    end
end
if isempty(n)
    n = size(colormap,1);
end
if isempty(clim)
    clim = caxis;
end
cmin = min(clim);
cmax = max(clim);
% np controls the visual flattening of extreme regions.
% It is the number of color levels used in these extreme regions.
% 1/8 of the range (0 to 1) uses a different gradient slope. 
% Extreme regions are at max red, min red, max blue, min blue.
% np=n/16 provides a uniform gradient but is less appealing visually.
np = max(round(sqrt(n)/4),1);   % scales slowly with n               
switch true
    case cmin >= 0      % display just blue
        ntot = round(n*cmax/(cmax-cmin));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        n2 = ntot - 2*np;
        v = [linspace(1,.875,np),linspace(.875,.125,n2),...
            linspace(.125,0,np)];
        y = repmat(v',[1,3]);
        y(:,3) = 1;
        y(1:ntot-n,:) = [];
    case cmax <= 0      % display just red
        ntot = round(n*abs(cmin/(cmax-cmin)));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        n2 = ntot - 2*np;
        v = [linspace(0,.125,np),linspace(.125,.875,n2),...
            linspace(.875,1,np)];
        y = repmat(v',[1,3]);
        y(:,1) = 1;
        y(n+1:end,:) = [];
    otherwise           % display both red and blue
        ntot = round(n*2*max(cmax,abs(cmin))/(cmax-cmin));
        if ntot > ntotmax   % check if cmax is close to cmin
            ntot = ntotmax;
        end
        m = fix(ntot/2);
        n2 = round((ntot-4*np)/2);
        v = [linspace(-1,-.875,np),linspace(-.875,-.125,n2),...
            linspace(-.125,.125,2*np),linspace(.125,.875,n2),...
            linspace(.875,1,np)];
        y = 1-abs(repmat(v',[1,3]));
        y(end-m+1:end,3) = 1;   % set blue
        y(1:m,1) = 1;           % set red
        if cmax > abs(cmin)
            y(1:ntot-n,:) = [];
        else
            y(n+1:end,:) = [];
        end   
end
if black
    y2 = y;
    y(:,1) = (1-y2(:,2)).*(y2(:,1)==1);
    y(:,3) = (1-y2(:,2)).*(y2(:,3)==1);
    y(:,2) = 0;
end
end
