%% MainArtoExample

% This script demonstrates the Automatic Rejection of Temporally-invariant
% Outliers (ARTO) algorithm for background phase offset correction (BPO) in cardiac
% phase contrast MRI as described in the (currently under review by Magnetic 
% Resonance in Medicine) manuscript entitled "A method to correct background 
% phase offset for phase-contrast MRI in the presence of steady flow and 
% spatial wrap-around artifact." 
%
% Available BPO correction methods are weighted least squares (WLS),
% weighted least squares w/ ARTO (WLS+ARTO), weighted regularized least
% squares (WRLS), and weighted regularized least squares w/ ARTO
% (WRLS+ARTO). Parameters in the code reflect parameters described in the
% Theory and Methods sections of the main text. 
%
% Two example datasets are included: 1) Ascending Aorta and 2) Main
% pulmonary artery, both with spatial wrap-around artifact. These corresond
% to the examples shown in Figures 4 and 5 respectively.
%
% Written by Aaron Pruitt. Last updated: 07/08/2018
%
clear all; close all;
%% Options
opt.fit         = 'WRLS+ARTO';  % type of fit, ['WLS','WLS+ARTO','WRLS', or 'WRLS+ARTO']
opt.pOrd        = 2;            % polynomial order for fitting, [1, 2, 3, or 4]
opt.lam         = 50;           % regularization strength for L1-norm regularization
opt.mTh         = 0.04;         % thresholding for magnitude mask (fraction of maximum)
opt.midFOVFrac  = 0.5;          % Center fraction of pixels fit in midPE-FOV for the 0th iteration of ARTO.
opt.Kmax        = 2;            % Maximum number of ARTO iterations
opt.tau         = 3;            % Exclusion threshold for ARTO
opt.delta       = 2;            % Minimum separation b/w Gaussian components in GMM
opt.gmmComp     = 3;            % Number of components fit by GMM

%% Load dicoms
% Choose path to "Datasets" folder
dirName = 'C:\Users\pruit\Box Sync\Research\Code\Background Phase\Code\ARTO_Example\Datasets\';
% Dataset 1 (Figure 4)
fName{1} = [dirName 'Dataset 1\Magnitude\'];                   % Magnitude
fName{2} = [dirName 'Dataset 1\Phase\'];                       % Phase
cName    = [dirName 'Dataset 1\Contour\Contour_Dataset1.mat']; % Contour
% Dataset 2 (Figure 5)
% fName{1} = [dirName 'Dataset 2\Magnitude\'];                 % Magnitude
% fName{2} = [dirName 'Dataset 2\Phase\'];                     % Phase
% cName    = [dirName 'Dataset 2\Contour\Contour_Dataset2.mat'];% Contour 

% Read dicom series from selected directory
[mag_cine, ~] = readDicomFolder(fName{1});                     % Magnitude
[phase_cine, acqParam_mat] = readDicomFolder(fName{2});        % Phase

Contour = readContourSegment(cName);
opt.acqParam = getAcqParam(acqParam_mat);

% Express phase as a percentage of the range 
phase_cine = (double(phase_cine)-2048)/2048;

%% Pre-processing
% Average phase cine over time
phi = mean(phase_cine,3);
% Temporal standard deviation
sigma = std(phase_cine,0,3);
% Magnitude thresholding
mag_avg = mean(mag_cine,3);
M = ones(size(mag_avg));
M(mag_avg <= opt.mTh*(max(max(mag_avg)))) = 0;
M = logical(M);

%% BPO Correction
switch opt.fit
    case 'WLS'
        tic; [phi_Cor,cMap] = wls(phi,sigma,M,opt); compTime = toc;  
    case 'WLS+ARTO'      
        tic; [phi_Cor,cMap] = wls_arto(phi,sigma,M,opt); compTime = toc;       
    case 'WRLS'       
        tic; [phi_Cor,cMap] = wrls(phi,sigma,M,opt); compTime = toc;       
    case 'WRLS+ARTO'       
        tic; [phi_Cor,cMap] = wrls_arto(phi,sigma,M,opt); compTime = toc;      
end
%% FLOW QUANTIFICATION

% Uncorrected
NetFlow_UC = flow_quant(phase_cine, Contour, opt);

% After BPO Correction
phase_cine_corrected = bsxfun(@minus,phase_cine,cMap);
NetFlow_Cor = flow_quant(phase_cine_corrected, Contour, opt);


% Output final result
disp(['==================================']);
disp(['Net Volumetric Flow:']);
disp(['Uncorrected: ' num2str(round(NetFlow_UC,2)) ' mL']);
disp([opt.fit,':  ' num2str(round(NetFlow_Cor,2)) ' mL']);
disp(['Computation time for ' opt.fit ': ' num2str(round(compTime,2)) ' s']);
disp(['==================================']);


