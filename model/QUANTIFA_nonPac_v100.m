%% Quantile Analysis of Temperature using Individual Foraminiferal Analyses
% Tropical Atlantic and Indian Ocean Variant
%
% ---- About ----
%
% QUANTIFA is a user-friendly IFA proxy system model that combines routines
% for modeling the sensitivity of IFA populations to changes in annual and 
% interannual climate variability with tools for processing, plotting, and 
% interpreting IFA-Mg/Ca data. Please see the README document for a 
% quickstart guide on how to use this algorithm. For a more detailed 
% description of the algorithm's statistical framework, please see our 
% associated publication cited below.
% 
% QUANTIFA is designed to work with potential temperature data from the  
% Ocean Reanalysis System 5 data assimilation (ORA-S5). This particular
% version of the script is designed to work with data from both the Tropical
% Atlantic (TA_ORAS5.mat) or Indian Oceans (TI_ORAS5.mat). Both of these
% datasets can be downloaded in the QUANTIFA repository on GitHub 
% (https://github.com/rh-glaubke/QUANTIFA). For any information regarding 
% the ORA-S5 dataset, please refer to Zuo et al. [2019] 
% (doi:10.5194/os-15-779-2019) or visit the Ocean Synthesis/ Reanalysis 
% Directory of the Integrated Climate Data Center: 
% https://icdc.cen.uni-hamburg.de/daten/reanalysis-ocean/easy-init-ocean/ecmwf-oras5.html
%
% ---- Author Info and Citation ----
% 
% Written by Ryan Glaubke (rglaubke@marine.rutgers.edu) and Kaustubh 
% Thirumalai (kaustubh@email.arizona.edu)
% Latest Update: January 2021
% 
% Citation: Glaubke et al. (2021). Discerning changes in high-frequency
%           climate variability using geochemical populations of individual 
%           foraminifera. Paleoceanography and Paleoclimatology, XXXXXX.
%           doi: 10.1029/2020PA004065

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
%     The ORA-S5 data structure (TA_ORAS5.mat for the Atlantic;
%     TI_ORAS5.mat for the Indian) contains a gridded field of potential 
%     temperature data ('votemper'; 1°x1° horizontal resolution) for 75 
%     analyzed depth levels ('depths'). It is available for download here: 
%     https://github.com/rh-glaubke/QUANTIFA
%
% Location of Core:
%     lat - Latitude for temperature gridpoint (Format: XX.5)
%     lon - Longitude for temperature gridpoint (Format: XX)
%     dep - Mean habitat depth for foram species of interest (Format: XX)
%
% Calibration Equation:
%     Select from a small bank of pre-loaded Mg/Ca calibration equations
%     (listed below), or program your own within the algorithm. Includes
%     the option to incorporate a dissolution correction term (D).
%
% Input Parameters:
%     A set of input conditions used to initialize the bootstrap picking
%     algorithm and parameterize the modeled pseudoproxy time series.
%     Breif description of each term is below.
%
% Output Products:
%     QUANTIFA produces the following outputs:
%     1) A "conformity contour plot" showing the relative influence of
%        annual and interannual climate variability on quantiles in the 
%        body and tails of the paleotemperature distributions.
%     2) A Q-Q plot of the user's data (or the user's data vs. ORA-S5 data).
%     3) An array, 'fp', which reports the false positive rate of each
%        computed quantile (i.e. probability each quantile commites a type
%        I error).
%     4) Data-model consistency maps that illustrate the proportion of
%        significant quantiles that align with each hypothetical climate
%        scenario. Darker colors = scenarios most consistent with the users
%        data.
%     
%     There are also a suite of supplemental figures, including:
%     1) The modeled modern time series (the modeled paleodata which 
%        represent modern climate variability), plotted alongside time 
%        series of select altered climate scenarios.
%     2) A Q-Q plot of the modeled modern time series and the select
%        altered time series.
% 
% FOR THE TROPICAL ATLANTIC DATASET:
% Horizontal Domain:  [29.5°N-29.5°S, 100°W-15°E; Convention: -ve for S,W]
% Vertical Domain:    [0.5-5902 m; See 'depths']
%
% FOR THE TROPICAL INDIAN DATASET:
% Horizontal Domain:  [29.5°N-29.5°S, 40°E-120°W; Convention: -ve for S,W]
% Vertical Domain:    [0.5-5902 m; See 'depths']

clc;
clear;
tic

% ---- Individual Foraminiferal Datasets ----

%X = 'Insert X Population Here';
%Y = 'Insert Y Population Here';

% ---- Reanalysis Dataset ----

load('TA_ORAS5.mat');           % Call up ORA-S5 dataset 
                                % (download link for .mat file above)

% ---- Location of Core ----

lat = 11.5;                     % Data point for VM12-107 (S Caribbean Sea)
lon = -67;
dep = 1;                       % NOTE: QUANTIFA selects the ORA-S5 depth
                                % level closest to the desired depth
                                % entered here. Upper ocean depth
                                % resolution is high, but decreases below
                                % ~1,000 m depth (see ORA-S5 documentation 
                                % for more information; link above).

                                
% ---- Calibration Equation ----

eqn = 1;                        % 1 = Anand et al. (2003) All species eqn
%D = 2.8;                        % 2 = Anand et al. (2003) G. ruber eqn
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

num_p = 70;                     % No. of "forams" picked from time series
num_q = 50;                     % No. of quantiles calculated
seas = (1:1:12);                % Seasonal weight for picking algorithm (e.g.
                                % (1:1:12) for all-year; [12 1 2] for DJF)
tsl = 60;                      % Length of synthetic time series (yrs)
run = 1000;                     % No. of Monte Carlo iterations
fp = 100;                       % No. of False Positive Exercises
cl = 90;                        % Confidence level (80, 90, 95, 99)
anerr = 0.1;                  % Analytical error (mmol/mol)
calerr = 1.2;                   % Calibration error (°C) [OPTIONAL]

% ---- Output Products ----

savefig = 0;                    % 0 = Do not save figures (PDF)
                                % 1 = Save figures (PDF)
                                % (Will save to working directory)

% ---- Interannual Event Years ----

% NOTE: You MUST manually input the years of known positive and negative
% phase events for the interannual climate mode you are interested in
% (i.e. the Atlantic Nino or the Indian Ocean Dipole). Also, you will want
% to change the starting month for the event year QUANTIFA averages over 
% (line 372). As a generic placeholder, we have entered El Nino and La Nina 
% years since 1958, and assigned the start month to April.

pyear = [1963;1965;1968;1972;1982;1986;1987;1991;1994;1997;2002;2009;2015];
nyear = [1970;1971;1973;1975;1984;1988;1995;1998;1999;2000;2007;2010;2011];
                                
                                                               
% ========================================================================
% ================= Modify Sections Below with Care! =====================
% ======================================================================== 


%% Input Parsing
% Checks input variables to ensure proper formatting.

disp('------------------------------------------------------------------------');
disp('Quantile Analysis of Temperature using Individual Foraminiferal Analyses');
disp('(QUANTIFA) - Tropical Atlantic and Indian Variant');
disp(' ');
disp('Written by Ryan Glaubke and Kaustubh Thirumalai'); 
disp('Last Updated: January 2021');
disp('Glaubke et al. (2021), Paleoceanography & Paleoclimatology');
disp('doi:10.1029/2020PA004065');
disp('------------------------------------------------------------------------');

% ---- Check Inputs ----

% -- IFA Data --                    % Check IFA Inputs

if exist('X','var') == 0 && exist('Y','var') == 0
    disp(' ');
    disp('NOTE: No IFA data detected. Running sensitivity tests for specified');
    disp('      location and depth...');
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

if rem(lat,1) == 0                  % Check proper formatting: latitude
    disp(' ');
    error('Error: Latitude is a whole number. Please input latitude as XX.5');
end

if (lat > 29.5) || (lat < -29.5)    % Check model domain: latitude
    disp(' ');
    error('Error: Entered latitude is out of model domain [29.5°N - 29.5°S]');
end

if rem(lon,1) ~= 0                  % Check proper formatting: longitude
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

[~,idx] = min(abs(depths-dep));     % Find nearest depth level to inputted depth
nrst_depth = depths(idx);
disp(' ');
disp(['NOTE: Entered depth = ',num2str(dep),' m']);
disp(['      Nearest depth level = ',num2str(nrst_depth),' m']);
disp(' ');

% -- Mg/Ca-T Equation --

if eqn == 3 || eqn == 5               % Check for dissolution correction term    
    if exist('D','var') == 0
        error('Error: no dissolution correction term found! Please assign the correction term native to your selected calibration equation under the variable "D".');
    end
end

% -- Input Parameters --

if exist('calerr','var') == 0
    disp('NOTE: No calibration error detected. Uncertainty estimation will');
    disp('      not take calibration error into account.');
    disp(' ');
    calerr = 0;
end

if anerr == 0
    error('Error: No analytical error detected. Please assign a value. If you are conducting a sensitivity analysis, a general value of 0.1 mmol/mol should be sufficient.');
end

if (cl ~= 80) && (cl ~= 90) && (cl ~= 95) && (cl ~= 99)
    error('Error: Unrecognized confidence level. Please enter either 80%, 90%, 95%, or 99%.');
end

%% Instrumental Data Extraction
% Extracts temperature time series data from ORA-S5 (Zuo et al. [2019])
% from specified lat/lon/dep parameters.

% ---- ORA-S5 lat/lon/depth indicies ----

if exist('TA_ORAS5', 'var')
    [ilat,ilon,idepth] = ta_oraco(lat,lon,nrst_depth,depths);
elseif exist('TI_ORAS5', 'var')
    [ilat,ilon,idepth] = ti_oraco(lat,lon,nrst_depth,depths);
end

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

% ---- Extract Temperature Time Series ----

for i=1:len
    if exist('TA_ORAS5', 'var')
        T(i) = TA_ORAS5.votemper(ilon,ilat,idepth,i);
    elseif exist('TI_ORAS5', 'var')
        T(i) = TI_ORAS5.votemper(ilon,ilat,idepth,i);
    end
end

if isempty(T(isnan(T))) == 0    % Replace NaNs at end by splicing time series
    le = length(T(isnan(T)));
    T(len-(le-1):len) = T(1:le);
end

% NOTE: only a few locations within the model domain have NaNs at the end
% of the time series which are remedied by the code above. If using a
% different reanalysis dataset, this will need to be adjusted accordingly.

% ---- Time Series Check ----

if nnz(~isnan(T)) == 0
   disp(' ');
   error('Error: Land Ho! It looks like there are no temperature data in your indexed time series... maybe youve struck land? Try another set of coordinates.');
end

%% Foraminiferal Mg/Ca Forward Model
% Transfers temperature record into Pseudo-Mg/Ca space.

MgCa = zeros(len,1);

switch eqn
    case 1  % Anand et al. (2003) all planktic species eqn
        for i = 1:len
            MgCa(i) = 0.38*exp(0.09*T(i));
        end
    case 2  % Anand et al. (2003) G. ruber eqn
        for i = 1:len
            MgCa(i) = 0.34*exp(0.102*T(i));
        end
    case 3  % Dekens et al. (2002) G. ruber eqn (w/ dissolution correction)
        for i = 1:len
            MgCa(i) = 0.38*exp(0.09*(T(i)-(0.61*D)-1.6));
        end 
    case 4  % Anand et al. (2003) G. sacculifer eqn
        for i = 1:len
            MgCa(i) = 1.06*exp(0.048*T(i));
        end        
    case 5  % Dekens et al. (2002) N. dutertrei eqn (w/ dissolution correction)
        for i = 1:len
            MgCa(i) = 0.60*exp(0.08*(T(i)-(2.8*D)-5.4));
        end
    case 6  % Mohtadi et al. (2011) G. tumida eqn
        for i = 1:len
            MgCa(i) = 0.410*exp(0.068*T(i));
        end
    case 7  % Anand et al. (2003) P. obliquiloculata eqn
        for i = 1:len
            MgCa(i) = 0.18*exp(0.12*T(i));
        end        
end

%% Climatology
% Extracts the seasonal and interannual climatology from forward-modeled 
% Mg/Ca time series.

% ---- ENSO Events in Record ----

pen = (pyear-syear)*12 + 4;            % + 4 => Begin w/ April
nen = (nyear-syear)*12 + 4;

pin = (pyear - syear) + 1;             % Year for injecting ENSO event
nin = (nyear - syear) + 1;

% ---- Declaration ----

mp = zeros(12,length(pyear));
mn = zeros(12,length(nyear));

Cp = zeros(12,1);
Cpsd = zeros(12,1);

Cn = zeros(12,1);
Cnsd = zeros(12,1);

CC = zeros(12,1);

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
lats = reshape(mn,length(nyear)*12,1);

% -- Mean Interanual Climatologies --

for month=1:12
    Cp(month) = nanmean(pt(month:12:end));
    Cn(month) = nanmean(lats(month:12:end));

    Cpsd(month) = std(pt(month:12:end));
    Cnsd(month) = std(lats(month:12:end));
end

% ---- Mean Base Climatology (Seasonality) ----

l = 1;
for month=4:15                          % 4-15 => April to March
    CC(l) = nanmean(MgCa(month:12:end));
    l = l + 1;
end

% ---- Interannual Deviations from Base Climatology ----

dp = Cp - CC;
dn = Cn - CC;

% ---- Interannual Relationship Check ----

if (nanmean(Cp) < nanmean(Cn))          % eg. in SW Pacific, ENSO = cold
    wp = 0;                             % and dry
else
    wp = 1;
end

%% Synthetic Time Series: Modeled Modern Variability
% Constructs synthetic time series mimicking modern variability. Read and
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
        if (wp == 0)
            mod_clim(:,pin(i)) = Cp - rand*Cpsd*mag;
        else
            mod_clim(:,pin(i)) = Cp + rand*Cpsd*mag;
        end
    end
    
    % ---- Insert Negative Interannual Events ----
    
    for i = 1:length(nin)
        if (wp == 0)
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
% amplitudes. Annual/interannual amplitudes vary from 100% dampened (or no
% signal) to 100% amplified (or double signal). 21 hypothetical seasonal
% cycle scenarios x 21 hypothetical interannual scenarios = 441 total 
% scenarios (440 altered and 1 modern).

% ---- Declaration ----

alt_series = zeros(21,21,tmslice);      % Altered Time Series

% ---- Simulate Altered Time Series ----

r = 1;                          
for p = -100:10:100                         % Seasonal Cycle Amplitude
    c = 1;
    for w = -100:10:100                     % Interannual Amplitude

        % ---- Change Seasonal Cycle ----

        sc = (CC - mean(CC))*(1+(p/100));   % Change seasonal cycle
        nCC = mean(CC) + sc;                % New seasonal cycle

        % ---- Change Interannual Amplitude ----

        amp = (w/100)+1;                    % Change in ENSO amplitude
        nCp = nCC + amp.*dp;                % New (altered) Positive Phase
        nCn = nCC + amp.*dn;                % New (altered) Negative Phase

        % ---- Insert Altered Base Climatology (Seasonality) ----

        alt_clim = zeros(12,years);
        for i = 1:years
            alt_clim(:,i) = (nCC(:));
        end

        % ---- Insert Altered El Nino Events ----

        for i = 1:length(pin)
            if (wp == 0)
                alt_clim(:,pin(i)) = nCp - rand*Cpsd*mag;
            else
                alt_clim(:,pin(i)) = nCp + rand*Cpsd*mag;
            end
        end

        % ---- Insert Altered La Nina Events ----

        for i = 1:length(nin)
            if (wp == 0)
               alt_clim(:,nin(i)) = nCn + rand*Cnsd*mag;
            else
               alt_clim(:,nin(i)) = nCn - rand*Cnsd*mag;
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

% ---- For Modeled Modern Time Series ----

mod_qs = zeros(num_q,1);

% -- Initiate Picking Loop --

r = rand(num_p,run);

tmslice_index = (1:tmslice)';           % Seasonal weighting
seas_wt = mod(tmslice_index,12);
seas_wt(seas_wt==0)=12;
seas_tmslice = tmslice_index(ismember(seas_wt,seas));

k = ceil(length(seas_tmslice)*r);
seas_mod_series = mod_series(ismember(seas_wt,seas));

MC = seas_mod_series(k) + anerr*randn;   % "Picked" Mg/Ca Values

% -- Compute Centered Quantiles --

switch eqn
    case 1
        mc_mod = (log(MC/0.38)/0.09) + calerr*randn;
    case 2
        mc_mod = (log(MC/0.34)/0.102) + calerr*randn;
    case 3
        mc_mod = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
    case 4
        mc_mod = (log(MC/1.06)/0.048) + calerr*randn;
    case 5
        mc_mod = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) + calerr*randn;
    case 6
        mc_mod = (log(MC/0.41)/0.068) + calerr*randn;
    case 7
        mc_mod = (log(MC/0.18)/0.12) + calerr*randn;
end

for i = 1:run
    q = quantile(mc_mod(:,i)-mean(mc_mod(:,i)),num_q);
    mod_qs(:,i) = q;
end
m_mod_qs = mean(mod_qs,2);

% ---- For Altered Time Series ----

alt_qs = zeros(21,21,num_q,run);    % Quantiles of each altered series
yr_alt = zeros(21,21,num_q,run);    % Y Residuals b/w quantiles and 1:1 line

% -- Initiate Picking Loop --

r = rand(21,21,num_p,run);

for p = 1:21
    for w = 1:21
       
        as = squeeze(alt_series(p,w,:));
        
        tmslice_index = (1:tmslice)';           % Seasonal weighting
        seas_wt = mod(tmslice_index,12);
        seas_wt(seas_wt==0)=12;
        seas_tmslice = tmslice_index(ismember(seas_wt,seas));
        
        k = ceil(length(seas_tmslice)*r(p,w,:,:));
        seas_alt_series = as(ismember(seas_wt,seas));

        MC = seas_alt_series(k) + anerr*randn;   % "Picked" Mg/Ca Values
        
        switch eqn
            case 1
                mc_alt = (log(MC/0.38)/0.09) + calerr*randn;
            case 2
                mc_alt = (log(MC/0.34)/0.102) + calerr*randn;
            case 3
                mc_alt = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
            case 4
                mc_alt = (log(MC/1.06)/0.048) + calerr*randn;
            case 5
                mc_alt = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) ...
                    + calerr*randn;
            case 6
                mc_alt = (log(MC/0.41)/0.068) + calerr*randn;
            case 7
                mc_alt = (log(MC/0.18)/0.12) + calerr*randn;
        end
        
        % -- Compute Centered Quantiles --
        
        for i = 1:run
            q = quantile(mc_alt(:,:,:,i)-mean(mc_alt(:,:,:,i)),num_q);
            alt_qs(p,w,:,i) = q;
        end
        
        % -- Calculate Residuals --

        for i = 1:num_q
            asq = squeeze(alt_qs(p,w,:,:));
            dist = asq(i,:) - mod_qs(i,:);
            yr_alt(p,w,i,:) = dist;
        end
        
    end
end

% ---- Generate Bins ----

bins = quantile(m_mod_qs-mean(m_mod_qs),num_q)';                
d = diff(bins)/2;                 % Bins defined by quantiles of MMV

edges = zeros(1,length(bins)+1);  % Bin edges are mid-way b/w MMV quantiles

for i = 2:length(bins)
    edges(:,i) = bins(i-1,:) + d(i-1,:);
end
edges(:,1) = min(mod_qs(1,:));
edges(:,end) = max(mod_qs(num_q,:));

% ---- Bin Data ----

m_alt_qs = zeros(21,21,num_q);      % Mean quantiles for each comparison
m_resids = zeros(21,21,num_q);      % Corresponding mean residuals
sd_alt_qs = zeros(21,21,num_q);     % Standard deviation within each bin

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
    end
end

%% Modeling IFA Detection Sensitivity
% Tests the sensitivity of the body and the tails of paleotemperature 
% distributions (at the user-defined location) to changes in the amplitude 
% of the Seasonal Cycle and interannual climate change.

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
    end
end

%% Output Product #1: Conformity Contour Plot
% Produces three contour plots of % conformity values for each altered
% climate scenario w/ respect to modeled modern time series. Tests the
% sensitivty of the tails (upper and lower 16%) and the central body
% (middle 68%) of the distribution.

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
    xlabel('Change in Interannual Climate Amplitude (%)','FontName','Myriad Pro','fontsize',14,'Fontweight','demi');
    ylabel('Change in Annual Climate Amplitude (%)','FontName','Myriad Pro','fontsize',14,'Fontweight','demi');
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

if NoIFA == 1                       % Terminates QUANTIFA if no IFA data
    toc                             % are present
    return
end

% ************************************************************************
% ************************************************************************
% ************************************************************************

%% False Positive Rate Estimation
% Computes the probability of committing a type I error for each individual
% quantile by repeatedly subsampling the modeled modern time series and
% comparing the two pseudo-IFA populations against one another via Q-Q
% analysis.

% ---- Declaration ----

fpr = zeros(num_q, fp);             % False Positive Rates

if exist('X','var') == 0            % Assign ORA-S5 Data to X if no IFA
    X = T;
    numX = num_p;
else                                   
    numX = length(X);
end
numY = length(Y);

for kk = 1:fp
    
    fp_qx = zeros(num_q,run);          % Quantiles from Pseudo-IFA
    fp_qy = zeros(num_q,run);          % Quantiles from Pseudo-IFA

    fpm_qy = zeros(num_q,1);           % Mean Y Quantile
    fp_serr = zeros(num_q,1);          % Total Sampling Error (X and Y)
    
    % ---- Picking Algorithm ----

    % -- For X Population --

    r = rand(numX,run);

    tmslice_index = (1:tmslice)';      % Seasonal Weighting Bloc
    seas_wt = mod(tmslice_index,12);                            % Find the modulus of 12 for each value to designate months
    seas_wt(seas_wt==0)=12;                                     % Assign 'December' as 12 and not 0.
    seas_tmslice = tmslice_index(ismember(seas_wt,seas));       % What are the index values required?

    k = ceil(length(seas_tmslice)*r);                           % Pick random values within the index
    seas_mod_series = mod_series(ismember(seas_wt,seas));       % What are the available timeseries values for those indices?

    MC = seas_mod_series(k) + anerr*randn;   % "Picked" Mg/Ca Values

    switch eqn
        case 1
            mc_x = (log(MC/0.38)/0.09) + calerr*randn;
        case 2
            mc_x = (log(MC/0.34)/0.102) + calerr*randn;
        case 3
            mc_x = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
        case 4
            mc_x = (log(MC/1.06)/0.048) + calerr*randn;
        case 5
            mc_x = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) + calerr*randn;
        case 6
            mc_x = (log(MC/0.41)/0.068) + calerr*randn;
        case 7
            mc_x = (log(MC/0.18)/0.12) + calerr*randn;
    end

    for i = 1:run
        q = quantile(mc_x(:,i)-mean(mc_x(:,i)),num_q);
        fp_qx(:,i) = q;
    end

    % -- For Y Population --

    r = rand(numY,run);

    tmslice_index = (1:tmslice)';      % Seasonal Weighting Bloc
    seas_wt = mod(tmslice_index,12);                            % Find the modulus of 12 for each value to designate months
    seas_wt(seas_wt==0)=12;                                     % Assign 'December' as 12 and not 0.
    seas_tmslice = tmslice_index(ismember(seas_wt,seas));       % What are the index values required?

    k = ceil(length(seas_tmslice)*r);                           % Pick random values within the index
    seas_mod_series = mod_series(ismember(seas_wt,seas));       % What are the available timeseries values for those indices?

    MC = seas_mod_series(k) + anerr*randn;   % "Picked" Mg/Ca Values

    switch eqn
        case 1
            mc_y = (log(MC/0.38)/0.09) + calerr*randn;
        case 2
            mc_y = (log(MC/0.34)/0.102) + calerr*randn;
        case 3
            mc_y = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
        case 4
            mc_y = (log(MC/1.06)/0.048) + calerr*randn;
        case 5
            mc_y = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) + calerr*randn;
        case 6
            mc_y = (log(MC/0.41)/0.068) + calerr*randn;
        case 7
            mc_y = (log(MC/0.18)/0.12) + calerr*randn;
    end

    for i = 1:run
        q = quantile(mc_y(:,i)-mean(mc_y(:,i)),num_q);
        fp_qy(:,i) = q;
    end

    % ---- Generate Bins ----

    b = mean(fp_qx,2); % <-- Doesn't matter if I choose X or Y
                       %     (both are subsamples of MMV)

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

% ---- Declaration ----

mc_qx = zeros(num_q,run,2);         % Quantiles from Pseudo-IFA
mc_qy = zeros(num_q,run,2);         % Quantiles from Pseudo-IFA

m_qx = zeros(num_q,1);              % Mean X Quantiles
m_qy = zeros(num_q,1);              % Mean Y Quantiles

serrX = zeros(num_q,1);             % One sigma total uncertainty for X
serrY = zeros(num_q,1);             % One sigma total uncertainty for Y

% ---- Picking Algorithm and Uncertainty Estimation ----

% -- For X Population --

r = rand(numX,run,2);

tmslice_index = (1:tmslice)';           % Seasonal weighting
seas_wt = mod(tmslice_index,12);
seas_wt(seas_wt==0)=12;
seas_tmslice = tmslice_index(ismember(seas_wt,seas));

k = ceil(length(seas_tmslice)*r);
seas_mod_series = mod_series(ismember(seas_wt,seas));

MC = seas_mod_series(k) + anerr*randn;   % "Picked" Mg/Ca Values

switch eqn
    case 1
        mc_x = (log(MC/0.38)/0.09) + calerr*randn;
    case 2
        mc_x = (log(MC/0.34)/0.102) + calerr*randn;
    case 3
        mc_x = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
    case 4
        mc_x = (log(MC/1.06)/0.048) + calerr*randn;
    case 5
        mc_x = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) + calerr*randn;
    case 6
        mc_x = (log(MC/0.41)/0.068) + calerr*randn;
    case 7
        mc_x = (log(MC/0.18)/0.12) + calerr*randn;
end

for ii = 1:2
    for i = 1:run
        q = quantile(mc_x(:,i,ii)-mean(mc_x(:,i,ii)),num_q);
        mc_qx(:,i,ii) = q;
    end
end

% -- For Y Population --

r = rand(numY,run,2);

tmslice_index = (1:tmslice)';           % Seasonal weighting
seas_wt = mod(tmslice_index,12);
seas_wt(seas_wt==0)=12;
seas_tmslice = tmslice_index(ismember(seas_wt,seas));

k = ceil(length(seas_tmslice)*r);
seas_mod_series = mod_series(ismember(seas_wt,seas));

MC = seas_mod_series(k) + anerr*randn;   % "Picked" Mg/Ca Values

switch eqn
    case 1
        mc_y = (log(MC/0.38)/0.09) + calerr*randn;
    case 2
        mc_y = (log(MC/0.34)/0.102) + calerr*randn;
    case 3
        mc_y = ((log(MC/0.38)/0.09) + (0.61*D) + 1.6) + calerr*randn;
    case 4
        mc_y = (log(MC/1.06)/0.048) + calerr*randn;
    case 5
        mc_y = ((log(MC/0.6)/0.08) + (2.8*D) + 5.4) + calerr*randn;
    case 6
        mc_y = (log(MC/0.41)/0.068) + calerr*randn;
    case 7
        mc_y = (log(MC/0.18)/0.12) + calerr*randn;
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

for ii = 1:2                        % Calculate uncertainy for X data
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
        for i = 1:num_q             % Calculate uncertainy for Y data
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

% -- Generate Figure --

h = figure;
axis tight manual;
clf;

% -- Plot Heat Map --

hh = heatmap(ww,flipud(pp),flipud(pfit));

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

hh.GridVisible = 'off';
% hh.ColorbarVisible = 'off';
hh.FontName = 'Myriad Pro';
hh.FontSize = 12;
% hh.CellLabelColor = 'none';

xlabel('Change in ENSO Amplitude (%)');
ylabel('Change in Seasonal Cycle (%)');

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
   
% ---- Supp Figure 1: Modeled Modern and Select Altered Time Series ----

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
    title('2x Seasonality / 2x ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
subplot(7,1,3)
    a = squeeze(alt_series(11,21,:));
    plot(yrs,a,'color',[1 0 0],'LineWidth',1.5);
    title('Normal Seasonality / 2x ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
subplot(7,1,4)
    a = squeeze(alt_series(1,21,:));
    plot(yrs,a,'color',[1 0.5 0.5],'LineWidth',1.5);
    title('No Seasonality / 2x ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
subplot(7,1,5)
    a = squeeze(alt_series(21,11,:));
    plot(yrs,a,'color',[0.4 0.53 1],'LineWidth',1.5);
    title('2x Seasonality / Normal ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
subplot(7,1,6)
    a = squeeze(alt_series(21,1,:));
    plot(yrs,a,'color',[0 0.22 1],'LineWidth',1.5);
    title('2x Seasonality / No ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');
subplot(7,1,7)
    a = squeeze(alt_series(1,1,:));
    plot(yrs,a,'color',[0 0.01 0.7],'LineWidth',1.5);
    title('No Seasonality / No ENSO','fontsize',14);
    ylabel('Mg/Ca (mmol/mol)','fontsize',12,'fontweight','bold');
    xlim([1 yrs(end)]);
    ylim([min(mod_series)-1 max(mod_series)+1]);
    set(gca,'fontweight','bold');    
xlabel('Months','fontsize',12,'fontweight','bold');
set(gcf,'position',[50,50,1200,1200]);
     
% ---- Supp Figure 2: Q-Q Plots of Select Altered Time Series ----

figure()
clf;

scatter(bins,squeeze(m_alt_qs(11,11,:)),32,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
scatter(bins,squeeze(m_alt_qs(21,21,:)),32,'o','filled','MarkerFaceColor',[0.6 0.06 0.06],'MarkerEdgeColor',[0.6 0.06 0.06]);
scatter(bins,squeeze(m_alt_qs(11,21,:)),32,'o','filled','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
scatter(bins,squeeze(m_alt_qs(1,21,:)),32,'o','filled','MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5]);
scatter(bins,squeeze(m_alt_qs(21,11,:)),32,'o','filled','filled','MarkerFaceColor',[0.4 0.53 1],'MarkerEdgeColor',[0.4 0.53 1]);
scatter(bins,squeeze(m_alt_qs(21,1,:)),32,'o','filled','MarkerFaceColor',[0 0.22 1],'MarkerEdgeColor',[0 0.22 1]);
scatter(bins,squeeze(m_alt_qs(1,1,:)),32,'o','filled','MarkerFaceColor',[0 0.01 0.7],'MarkerEdgeColor',[0 0.01 0.7]);
legend('Modeled Modern Variability','2x Seasonality / 2x ENSO','Normal Seasonality / 2x ENSO','No Seasonality / 2x ENSO','2x Seasonality / Normal ENSO','2x Seasonality / No ENSO','No Seasonality / No ENSO','location','northwest','box','off');

rl = refline(1,0);
rl.Color = 'k';
rl.LineWidth = 1.5;
rl.HandleVisibility = 'off';
set(gca,'LineWidth',1.5,'FontName','Myraid Pro','fontsize',14,'Fontweight','demi');
xlabel('Modern Population Mean Normalized Quantiles (°C)','FontName','Myraid Pro','Fontsize',14,'Fontweight','demi');
ylabel('Altered Population Mean Normalized Quantiles (°C)','FontName','Myraid Pro','Fontsize',14,'Fontweight','demi');

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

toc

%% Built-in Functions
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
