%% Example using the standard prediction
clear; close all; clc

%load the TEX data
PhanSST_v001 = readtable("nutrient-effect-on-TEX/spreadsheets/published_data/PhanSST_v001.csv");

PhanTEX = PhanSST_v001(PhanSST_v001.ProxyType == "tex",:);

sel_site = PhanTEX.SiteName == "NIOP-C2_905_PC";
NIOP905_TEX_da = PhanTEX(sel_site,:);

% Extract relevant columns for analysis
TEX_values = NIOP905_TEX_da.ProxyValue;
site_lat = unique(NIOP905_TEX_da.ModLat);
site_lon = unique(NIOP905_TEX_da.ModLon);

%% Set the inputs for the prediction code
dats=TEX_values;
lon=site_lon;
lat=site_lat;
prior_std=6;

%select which model to use:
%SST
runname='subT';

%predict:
Output_Struct = bayspar_tex(dats, lon, lat, prior_std, runname);

ages = NIOP905_TEX_da.Age;
P50_pred_subT = Output_Struct.Preds(:,2);
% plot AGE vs pred

plot(ages,P50_pred_subT);