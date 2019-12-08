function [fitresult, gof] = fitFrontalParietal(fxi_deg, fyi_pdk, fxo_deg, fyo_pdk, pxi_deg, pyi_pdk, pxo_deg, pyo_pdk)
%CREATEFITS(FXI_DEG,FYI_PDK,FXO_DEG,FYO_PDK,PXI_DEG,PYI_PDK,PXO_DEG,PYO_PDK)
%  Create fits.
%
%  Data for 'frontal-id-poly' fit:
%      X Input : fxi_deg
%      Y Output: fyi_pdk
%  Data for 'frontal-id-gauss' fit:
%      X Input : fxi_deg
%      Y Output: fyi_pdk
%  Data for 'frontal-od-poly' fit:
%      X Input : fxo_deg
%      Y Output: fyo_pdk
%  Data for 'frontal-od-gauss' fit:
%      X Input : fxo_deg
%      Y Output: fyo_pdk
%  Data for 'parietal-id-poly' fit:
%      X Input : pxi_deg
%      Y Output: pyi_pdk
%  Data for 'parietal-id-gauss' fit:
%      X Input : pxi_deg
%      Y Output: pyi_pdk
%  Data for 'parietal-od-poly' fit:
%      X Input : pxo_deg
%      Y Output: pyo_pdk
%  Data for 'parietal-od-gauss' fit:
%      X Input : pxo_deg
%      Y Output: pyo_pdk
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Sep-2019 21:49:35

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 8, 1 );
gof = struct( 'sse', cell( 8, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'frontal-id-poly'.
[xData, yData] = prepareCurveData( fxi_deg, fyi_pdk );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'frontal-id-poly' );
h = plot( fitresult{1}, xData, yData );
legend( h, 'fyi_pdk vs. fxi_deg', 'frontal-id-poly', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'fxi_deg', 'Interpreter', 'none' );
ylabel( 'fyi_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'frontal-id-gauss'.
[xData, yData] = prepareCurveData( fxi_deg, fyi_pdk );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.941176470588235 6 3.30745686740786];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'frontal-id-gauss' );
h = plot( fitresult{2}, xData, yData );
legend( h, 'fyi_pdk vs. fxi_deg', 'frontal-id-gauss', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'fxi_deg', 'Interpreter', 'none' );
ylabel( 'fyi_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'frontal-od-poly'.
[xData, yData] = prepareCurveData( fxo_deg, fyo_pdk );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'frontal-od-poly' );
h = plot( fitresult{3}, xData, yData );
legend( h, 'fyo_pdk vs. fxo_deg', 'frontal-od-poly', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'fxo_deg', 'Interpreter', 'none' );
ylabel( 'fyo_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'frontal-od-gauss'.
[xData, yData] = prepareCurveData( fxo_deg, fyo_pdk );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.941176470588235 6 3.09862011363996];

% Fit model to data.
[fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'frontal-od-gauss' );
h = plot( fitresult{4}, xData, yData );
legend( h, 'fyo_pdk vs. fxo_deg', 'frontal-od-gauss', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'fxo_deg', 'Interpreter', 'none' );
ylabel( 'fyo_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'parietal-id-poly'.
[xData, yData] = prepareCurveData( pxi_deg, pyi_pdk );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult{5}, gof(5)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'parietal-id-poly' );
h = plot( fitresult{5}, xData, yData );
legend( h, 'pyi_pdk vs. pxi_deg', 'parietal-id-poly', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'pxi_deg', 'Interpreter', 'none' );
ylabel( 'pyi_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'parietal-id-gauss'.
[xData, yData] = prepareCurveData( pxi_deg, pyi_pdk );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.923076923076923 2 2.4950456034583];

% Fit model to data.
[fitresult{6}, gof(6)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'parietal-id-gauss' );
h = plot( fitresult{6}, xData, yData );
legend( h, 'pyi_pdk vs. pxi_deg', 'parietal-id-gauss', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'pxi_deg', 'Interpreter', 'none' );
ylabel( 'pyi_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'parietal-od-poly'.
[xData, yData] = prepareCurveData( pxo_deg, pyo_pdk );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult{7}, gof(7)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'parietal-od-poly' );
h = plot( fitresult{7}, xData, yData );
legend( h, 'pyo_pdk vs. pxo_deg', 'parietal-od-poly', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'pxo_deg', 'Interpreter', 'none' );
ylabel( 'pyo_pdk', 'Interpreter', 'none' );
grid on

%% Fit: 'parietal-od-gauss'.
[xData, yData] = prepareCurveData( pxo_deg, pyo_pdk );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.846153846153846 3 2.1072525928219];

% Fit model to data.
[fitresult{8}, gof(8)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'parietal-od-gauss' );
h = plot( fitresult{8}, xData, yData );
legend( h, 'pyo_pdk vs. pxo_deg', 'parietal-od-gauss', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'pxo_deg', 'Interpreter', 'none' );
ylabel( 'pyo_pdk', 'Interpreter', 'none' );
grid on

