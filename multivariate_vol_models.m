%% **********************************************************
%  FM320 - Risk management and modelling
%  Michaelmas 2019 - summative 2
%  Date: 24 November 2019
%  **********************************************************

%% %% *** Part0: Preliminaries ***

% Clearing existing data from Matlab
clc;                   % Clears command window
clear all;             % Removes all variables from workspace
clf;                   % Clears current figure window
close all;             % Close all windows

%% *** Part 1: Loading Data from Excel File & Useful Information *** 

% Loading data
Data_Stocks            = xlsread('Homework 8 - Data.xlsx', ...
                         'Homework 8 Data', 'G12:AB4790');

AdjustedPrices         = Data_Stocks(:, 2:end);  
Dates                  = Data_Stocks(:, 1);

Tickers                = {'S&P500 ','AAPL','MSFT', 'AMZN', 'GOOG', 'FB', ...
                          'BRK.A', 'V', 'JPM', 'JNJ', 'WMT', 'PG', 'XOM', ...
                          'MA', 'T', 'BAC', 'HD', 'VZ', 'DIS', 'INTC', 'KO'};
                      
% Useful information & global variables
StartDate              = 1999;
EndDate                = 2018;
NDates                 = size(Dates,1);
NDatesForReturns       = NDates - 1;

NSecurities            = size(Data_Stocks,2)-1;    % including S&P500
NStocks                = NSecurities - 1;          % excluding S&P500 

NObsBurnIn             = 252;                      % number of burn in period
DaysPerYear            = 252;  

% Computing log returns
SimpleReturns          = AdjustedPrices(2:end, :) ./ AdjustedPrices(1:(end - 1), :) - 1;
LogReturns             = log(1 + SimpleReturns);
LogReturnsSq           = LogReturns .* LogReturns;
DatesForReturns        = Dates(2:end, :);

FigNo                  = 1;
%% *** Question (b) - Principal components analysis with stock data ***
Indices                = [10 11 12 15 17 18];

% Extracting rows
StockPCAData           = LogReturns(:, Indices);

% Performing principal components analysis (with 30-Year)
[Weights, ~, FactorVariance] = pca(StockPCAData);
% FactorVariance             = sum(FactorVariance);    % in windows version
% has to recover this line
TotalVar               = sum(FactorVariance);
FracVarExplained       = FactorVariance / TotalVar;

% Displaying principal component weights on each maturity
DisplayMat1            = [Tickers(1, Indices)' num2cell(Weights')];
DisplayMat2            = [num2str('Factor Variance') num2cell(FactorVariance)']; % in windows version recover num2cell(FactorVariance), without transpose
DisplayMat3            = [num2str('Fraction Total') num2cell(FracVarExplained)']; % without transpose
PCAResults             = [DisplayMat1; DisplayMat2; DisplayMat3];

%% %% *** Question (c): Choose stocks, estimate univariate GARCH as benchmark ***
% Compute univariate GARCH model to compare GARCH volatility with O-GARCH

% Choose JNJ, PG and HD to analyze DCC and OGARCH
Indices             = [10 12 17];
MVGARCHData         = LogReturns(:, Indices);
MVGARCHDataSq       = LogReturnsSq(:, Indices);

% Estimate GARCH for each series, save information for computing estimates
GARCHModel          = garch(1, 1);
GARCHInfo           = NaN(5, 3);

for i=1:3
   EstimatedModel   = estimate(GARCHModel, MVGARCHData(:, i));
   GARCHInfo(1, i)  = EstimatedModel.Constant;
   GARCHInfo(2, i)  = EstimatedModel.ARCH{1};
   GARCHInfo(3, i)  = EstimatedModel.GARCH{1};
   GARCHInfo(4, i)  = EstimatedModel.UnconditionalVariance;
   GARCHInfo(5, i)  = sqrt(DaysPerYear) * sqrt(GARCHInfo(4, i));
end

% Producing matrix to display results
Tickers2               = {'JNJ', 'PG', 'HD'}
RowHeaders             = {'Constant', 'ARCH(1)', 'GARCH(1)',...
                          'Variance', 'Unconditional Vol. (Ann.)'}
ColHeaders             = [{' '} Tickers2];
DisplayAux             = [RowHeaders' num2cell(GARCHInfo)];
GARCHDisplay           = [ColHeaders; DisplayAux];

% With the coefficient estimates, we can compute GARCH variances
GARCHVar            = NaN(NDatesForReturns, 3);
GARCHVar(1, :)      = GARCHInfo(4, :);            % first use unconditional variances

for i=2:NDatesForReturns
    GARCHVar(i, :)  = GARCHInfo(1, :) + ...
                      GARCHInfo(2, :) .* MVGARCHDataSq(i - 1, :) + ...
                      GARCHInfo(3, :) .* GARCHVar(i - 1, :);
end

% Convert to volatility
GARCHVol            = sqrt(GARCHVar);
GARCHVolJNJ         = GARCHVol(:, 1);
GARCHVolPG          = GARCHVol(:, 2);
GARCHVolHD          = GARCHVol(:, 3);

% Compute moving window correlations for comparison with models
CorrWindow          = 100;       % a window of lag 100
CorrLegend          = [num2str(CorrWindow) '-Day'];
MWCorr              = NaN(NDatesForReturns, 3);

for i=(CorrWindow+1):NDatesForReturns
    Indices         = (i-CorrWindow):(i-1);     
    Data            = MVGARCHData(Indices,:);
    Temp            = corrcoef(Data);
    MWCorr(i, 1)    = Temp(1, 2);    % the correlation btw stock 1 and stock 2 during the 100 datapoints (JNJ & PG)
    MWCorr(i, 2)    = Temp(1, 3);    % the correlation btw stock 1 and stock 3 (ie. JNJ and HD)
    MWCorr(i, 3)    = Temp(2, 3);    % the correlation btw stock 2 and stock 3 (ie. PG and HD)
end

%% *** Question (c): DCC with stock data ***
% Estimate DCC model
[DCCParms]          = dcc(MVGARCHData, [], 1, 0, 1);  % using GARCH (1,1) to estimate correlation matrix

% Construct DCC model matrices
DCCOmega            = [DCCParms(1) DCCParms(4) DCCParms(7)];   % GARCH(1,1) omega for volatility matrix
DCCAlpha            = [DCCParms(2) DCCParms(5) DCCParms(8)];   % ARCH term for volatility matrix
DCCBeta             = [DCCParms(3) DCCParms(6) DCCParms(9)];   % GARCH term for voalatility matrix
DCCRBar             = [1 DCCParms(10) DCCParms(11); DCCParms(10) 1 DCCParms(12); DCCParms(11) DCCParms(12) 1];  % its R bar in the equation
DCC_a               = DCCParms(13);   % apha   
DCC_b               = DCCParms(14);   % beta  

% Computing GARCH Volatilities for each asset, and standardized residuals
DCCVar              = NaN(NDatesForReturns, 3);
DCCVar(1, :)        = DCCOmega ./ (1 - DCCAlpha - DCCBeta);  % unconditional variance for GARCH model

for i=2:NDatesForReturns
    DCCVar(i, :)    = DCCOmega + ...
                      DCCAlpha .* MVGARCHDataSq(i - 1, :) + ...
                      DCCBeta .* DCCVar(i - 1, :);
end

DCCVol              = sqrt(DCCVar);
DCCEps              = MVGARCHData ./ DCCVol;  % Zts, the standardised residuals

% Compute value for DCC Qt and Rt, using recursive relation; 
% initialize with estimate for RBar
DCCQt               = NaN(3, 3, NDatesForReturns);  % Qt is the non-normalised Rt
DCCRt               = NaN(3, 3, NDatesForReturns);  % the normalised version (final outcome)
DCCQt(:, :, 1)      = DCCRBar;
DCCRt(:, :, 1)      = DCCRBar;

for i=2:NDatesForReturns
    DCCQt(:, :, i)  = (1 - DCC_a - DCC_b) * DCCRBar ...
                       + DCC_a * DCCEps(i-1, :)' * DCCEps(i-1, :) ...
                       + DCC_b * DCCQt(:, :, i-1);
    AuxMat          = [sqrt(DCCQt(1, 1, i))  sqrt(DCCQt(2, 2, i)) sqrt(DCCQt(3, 3, i))];   % the time-varying volatilities for each of the securities
    DCCRt(:, :, i)  = DCCQt(:, :, i) ./ (AuxMat' * AuxMat);   % normalising Qt to become Rt 
end

% Produce a proper correlation matrix
% Extracting volatility and correlation information
DCCVolJNJ           = DCCVol(:, 1);
DCCVolPG            = DCCVol(:, 2);
DCCVolHD            = DCCVol(:, 3);

DCCCorrJNJ_PG       = DCCRt(1, 2, :);
DCCCorrJNJ_PG       = reshape(DCCCorrJNJ_PG, NDatesForReturns, 1);

DCCCorrJNJ_HD       = DCCRt(1, 3, :);
DCCCorrJNJ_HD       = reshape(DCCCorrJNJ_HD, NDatesForReturns, 1);

DCCCorrPG_HD        = DCCRt(2, 3, :);
DCCCorrPG_HD        = reshape(DCCCorrPG_HD, NDatesForReturns, 1); 

% Produce charts comparing output of DCC with GARCH information
ChartDataJNJ        = [GARCHVolJNJ DCCVolJNJ];
ChartDataPG         = [GARCHVolPG  DCCVolPG];
ChartDataHD         = [GARCHVolHD  DCCVolHD];
ChartDataCorr       = [MWCorr DCCCorrJNJ_PG DCCCorrJNJ_HD DCCCorrPG_HD];

ChartDataJNJ        = ChartDataJNJ((NObsBurnIn+1):end, :);
ChartDataPG         = ChartDataPG((NObsBurnIn+1):end, :);
ChartDataHD         = ChartDataHD((NObsBurnIn+1):end, :);
ChartDataCorr       = ChartDataCorr((NObsBurnIn+1):end, :);

TimeLabels          = linspace(StartDate + 1, EndDate, NDatesForReturns - NObsBurnIn);

% Produce figure
figure(FigNo)
subplot(3, 1, 1);   % 3 graphs in total, and the first row and first column of graph
titleStr = 'Volatility Estimates for JNJ';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * ChartDataJNJ(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');    % plot GARCH estimate
plot(TimeLabels, sqrt(252) * ChartDataJNJ(:, 2), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');      % plot DCC volatility estimate
legend('GARCH', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);   % 3 graphs in total, and the second row and first column of graph
titleStr = 'Volatility Estimates for PG';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * ChartDataPG(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * ChartDataPG(:, 2), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);   % 3 graphs in total, and the third row and first column of graph
titleStr = 'Volatility Estimates for HD';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * ChartDataHD(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * ChartDataHD(:, 2), 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

FigNo = FigNo +1

% Produce figure for correlations estimates from GARCH and DCC
figure(FigNo)
subplot(3, 1, 1);
titleStr = 'Correlation Estimates for JNJ and PG';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, ChartDataCorr(:, 1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');    % correlation estimate using GARCH for JNJ&PG
plot(TimeLabels, ChartDataCorr(:, 4), 'LineStyle', '--', ...   % the fourth column bc MWCorr is composed of 3 columns
                 'LineWidth', 1, 'Color', 'red');      % correlation estimate using DCC for JNJ&PG
legend(CorrLegend, 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);
titleStr = 'Correlation Estimates for JNJ and HD';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, ChartDataCorr(:, 2), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');    % correlation estimate using GARCH for JNJ&HD
plot(TimeLabels, ChartDataCorr(:, 5), 'LineStyle', '--', ...   % the fifth column bc MWCorr is composed of 3 columns
                 'LineWidth', 1, 'Color', 'red');      % correlation estimate using DCC for JNJ&HD
legend(CorrLegend, 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);
titleStr = 'Correlation Estimates for PG and HD';
title(titleStr);
xlabel('Dates');
axis([(StartDate + 1) EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, ChartDataCorr(:, 3), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');    % correlation estimate using GARCH for PG&HD
plot(TimeLabels, ChartDataCorr(:, 6), 'LineStyle', '--', ...   % the sixth column bc MWCorr is composed of 3 columns
                 'LineWidth', 1, 'Color', 'red');      % correlation estimate using DCC for PG&HD
legend(CorrLegend, 'DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');
FigNo = FigNo + 1;

%% *** Question (c) - O-GARCH analysis for our model ***
% Stocks used: JNJ, PG and HD
Indices                = [10 12 17];
NStocks                = size(Indices, 2);  

% We first redo the PCA analysis for the three stocks only
% Extracting only rows without missing data
StockPCAData           = LogReturns(:, Indices);
DataIsMissing          = isnan(StockPCAData);
DataIsMissing          = sum(DataIsMissing, 2);
RowsToUse              = find(DataIsMissing == 0);
NDates                 = size(RowsToUse, 1);
StockPCAData           = StockPCAData(RowsToUse, :);

% Performing principal components analysis
[Weights, ~, FactorVariance] = pca(StockPCAData);

% FactorVariance        = sum(FactorVariance);   % remove the % in windows
TotalVar                = sum(FactorVariance);
FracVarExplained        = FactorVariance / TotalVar;

% Displaying principal component weights on each maturity
DisplayMat1            = [Tickers(1, Indices)' num2cell(Weights')];
DisplayMat2            = [num2str('Factor Variance') num2cell(FactorVariance)'];   % no transpose in Windows
DisplayMat3            = [num2str('Fraction Total') num2cell(FracVarExplained)'];  % no transpose in Windows
PCAResults2            = [DisplayMat1; DisplayMat2; DisplayMat3];

% The analysis seems to indicate one common factors, capturing market
% exposure, which is a commonality between JNJ PG and HD (retail and manufacturers)
NumFactors             = 1;
[Parms, HT, W, PC]     = o_mvgarch(StockPCAData, NumFactors, 1, 0, 1);
% HT: historical estimate for variance and covariances
% PC: the number of arch terms and leverage and garch terms i would like to
% have
% HT is three dimensional bc i have 3 diff securities 

% Producing displays comparing volatilities for individual securities from 
% GARCH and from O-GARCH
GARCHInfo1             = NaN(6, NStocks);
GARCHVol1              = NaN(NDates, NStocks);
figure(FigNo);

for i=1:NStocks
    % First, extract returns and OGARCH volatilities
    StockRet           = LogReturns(RowsToUse, Indices(1, i));
    StockRetSq         = StockRet .* StockRet;
    VolOGARCH          = sqrt(squeeze(HT(i, i, :)));  % squeeze convert it to a 2-dimensional matrix
                                                      % HT (i,i,:) extracts the volatility
                                                      % of security i
                                      
    % Next, estimate GARCH
    [Parameters]       = tarch(StockRet, 1 , 0, 1); 
    
    Omega              = Parameters(1, 1);
    Alpha              = Parameters(2, 1);
    Beta               = Parameters(3, 1);

    GARCHInfo1(1, i)   = Omega;
    GARCHInfo1(2, i)   = Alpha;
    GARCHInfo1(3, i)   = Beta;
    GARCHInfo1(4, i)   = Alpha + Beta;
    GARCHInfo1(5, i)   = Omega / (1 - Alpha - Beta);
    GARCHInfo1(6, i)   = sqrt(DaysPerYear) * sqrt(GARCHInfo1(5, i));
    
    % Compute GARCH conditional volatilities
    GARCHVar           = NaN(NDates, 1);
    GARCHVar(1, 1)     = GARCHInfo1(5, 1);

    for j=2:NDates
        GARCHVar(j, 1) = Omega + Alpha * StockRetSq(j - 1, 1) + ...
                           Beta * GARCHVar(j - 1, 1);
    end

    GARCHVol1(:, i)    = sqrt(GARCHVar);
    
    % Graph the volatility estimates, together with returns (dates omitted)
    subplot(3, 1, i);
    ChartData          = [StockRet 2*GARCHVol1(:, i) -2*VolOGARCH];
    AuxVar1            = DatesForReturns(RowsToUse(1, 1), 1);
    StartDate          = (AuxVar1 - mod(AuxVar1, 10000))/10000;    
    TimeLabels         = linspace(StartDate, EndDate, NDates);
    YTickPoints        = linspace(-0.25, 0.25, 11);

    % Charting the series
    figure(FigNo);
    titleStr = strcat({'Daily Returns and +/- 2 Std. Dev. GARCH and O-GARCH for '}, Tickers(1, Indices(i)));
    title(titleStr);
    xlabel('Dates');
    axis([StartDate EndDate -0.25 0.25]);
    grid on;
    hold on;
    plot(TimeLabels, ChartData(:, 1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'blue');
    plot(TimeLabels, ChartData(:, 2), 'LineStyle', '--', ...
                     'LineWidth', 1, 'Color', 'red'); 
    plot(TimeLabels, ChartData(:, 3), 'LineStyle', '--', ...
                     'LineWidth', 1, 'Color', 'green');
    set(gca, 'YTickMode', 'manual');
    set(gca, 'YTick', YTickPoints);
    set(gca, 'YTickLabel', num2str(100 .* get(gca, 'YTick')', '%1.0f%%'));              
    legend({'Log Returns', '2*GARCH Vol', '(-2)*O-GARCH Vol'}, ...
                         'location','best');            
    set(gcf, 'color', 'white');                  
end 

FigNo   = FigNo + 1;
% The O-garch is designed to test the common factor across securities and
% how it fits into this security, the GARCH is designed to fit the security
% well 

% Producing displays comparing correlations for pairs of securities from 
% from O-GARCH and moving window
figure(FigNo);
ChartLocation          = 1;

for i=1:(NStocks-1)
    for j=(i+1):NStocks    % i and j represent two different stocks
        % Extract O-GARCH correlations
        Variance_i     = squeeze(HT(i, i, :));  % extract the variance estimate of security i
        Variance_j     = squeeze(HT(j, j, :));
        Covariance     = squeeze(HT(i, j, :));
        CorrelOGARCH   = Covariance ./ sqrt(Variance_i .* Variance_j);
        
        % Computing moving window correlations
        MWCorr         = NaN(NDates, 1);
        CorrLegend     = [num2str(CorrWindow) '-Day'];
        CorrData       = StockPCAData(:, [i j]);
        
        for t=(CorrWindow+1):NDates
            CorrDates  = (t-CorrWindow):(t-1);
            Data       = CorrData(CorrDates,:);
            Temp       = corrcoef(Data);
            MWCorr(t, 1) = Temp(1, 2);   
        end 
        
        % Graph the correlation estimates
        subplot(3, 1, ChartLocation);
        ChartData      = [MWCorr CorrelOGARCH];
        ChartData      = ChartData((CorrWindow+1):end, :);
        AuxVar1        = DatesForReturns(RowsToUse(CorrWindow+1, 1), 1);
        StartDate      = (AuxVar1 - mod(AuxVar1, 10000))/10000;    
        TimeLabels     = linspace(StartDate, EndDate, NDates - CorrWindow);
        YTickPoints    = linspace(-1, 1, 9);

        % Charting the series
        figure(FigNo);
        titleStr = strcat({'Correlation - O-GARCH and Moving Window - '}, Tickers(1, Indices(i)), {' and '}, Tickers(1, Indices(j)));
        title(titleStr);
        xlabel('Dates');
        axis([StartDate EndDate -1.25 1.25]);
        grid on;
        hold on;
        plot(TimeLabels, ChartData(:, 1), 'LineStyle', '-' , ...
                          'LineWidth', 1, 'Color', 'blue');
        plot(TimeLabels, ChartData(:, 2), 'LineStyle', '--', ...
                         'LineWidth', 1, 'Color', 'red'); 
        set(gca, 'YTickMode', 'manual');
        set(gca, 'YTick', YTickPoints);
        set(gca, 'YTickLabel', num2str(get(gca, 'YTick')', '%.2f'));              
        legend({CorrLegend, 'O-GARCH'}, 'location','best');            
        set(gcf, 'color', 'white');

        ChartLocation  = ChartLocation + 1;
    end
end


% The correlations btw pairs of stocks, the model that captures the medium time trend
% in the blue line (100 day moving window correlation)is DCC, the OGARCH
% might not be a good model in the sense that:
% over/under estimate the correlations but not consistent
% sometimes there is huge gaps btw blue and red lines -- these have very
% high frequency components 
% to solve this, we can add time varying component for Var(residual)

