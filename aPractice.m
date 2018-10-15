function aPractice(subjID,subjVert,oriMean,acq,relativeSurr)

% TI_psi_2int(subjID,subjVert,oriMean,acq,relativeSurr)
%
% e.g. aPractice('K',0,0,1,0);
%
% This experiment shows 2 intervals of the tilt illusion, using a psi adaptive staircase procedure. The 
% participant's task is to judge the interval in which the central orientation appears more vertical. The
% surround annulus is always oriented at ±15 degrees, and the centre orientation differs depending on the
% parameters set in oriMean and oriDiff.
%
% The psi adaptive staircase procedure used here is described in detail in the following paper:
%    Bayesian adaptive estimation of psychometric slope and threshold
%    - Kontsevich LL & Tyler CW (1999) Vision Res 39, 2729-2737.
%    And was originally written in MATLAB by Colin Clifford on 28.09.01 (based on Java code by Michael Pianta)
%    and subsequently modified CC 02.03.06 for subjective vertical estimation
%
% The remaining (non-psi) components of this code were written by Matt Patten on 26 May, 2014.
%
% Inputs:
%      subjID  - A string consisting of the initials of the participant. This is used to save a unique results file.
%     subjVert - The angle where the subject perceives to be (in deg). Positive value indicates to the right of vertical.
%      oriMean - In degrees, the midpoint between the two centre orientations presented in each interval. For a mean
%                of zero (vertical), the two orientations are located equally on either side of this. The magnitude of
%                such depends on the "dif" staircase intensity value.
%          acq - The acquisition number for this run.
% relativeSurr - 0 or 1. 0 indicates surround stim will be given an absolute value (e.g., ±15 deg), while 1 indicates 
%                that the surround orientation will vary on each trial corresponding to changes in the centre orientation.
%----------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check inputted parameters and directory structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surroundLabel = {'abs','rel'};
% output filename
root_dir = [fileparts(mfilename('fullpath')) '\']; %get base experiment directory (wherever this mfile is located)
save_dir = [root_dir 'results\' subjID '\'];
saveName = [save_dir subjID '_results_psi_2int_' surroundLabel{relativeSurr+1} '_mean' num2str(oriMean) '_' num2str(acq) '.mat']; %the results filename

%if results directory doesn't exist, create it
if ~exist(save_dir) 
    mkdir(save_dir);
end

%if results file already exists, ask if we want to overwrite
if exist(saveName,'file')
    tryAgain = true; %keep looping until we give permission to continue (by pressing 'y')
    while tryAgain
        overwriteResp = input('Results already exist for this subject and acquisition. Continue? (y/n)','s'); %get user input from keyboard
        if strcmp(overwriteResp,'y') %if they type y, then escape loop and overwrite file
            tryAgain = false; %break loop
        elseif strcmp(overwriteResp,'n') %if they type n, then abort program
            error('Aborted by user to avoid overwriting results file');
        end
    end         
end 

KbName('UnifyKeyNames');
ResponseLeft = KbName('F'); % arrows
ResponseRight = KbName('J'); %
ResponseExit = KbName('P');
RestrictKeysForKbCheck([ResponseLeft ResponseRight ResponseExit]);
firstorsecond = 1;
fristorsecodn = 1; 
keypress = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global usingBitsSharp;
global eventCounter;
global eventLog;
global the_experiment_start;

global buttonAlreadyPressed;
global response;
global responseTime;
global pickedLeftStim;

%jheapcl; %clean up java memory

%setup parameters
usingBitsSharp = false; %true; %only true if we have the bits# device attached to the monitor
overlay = false; %if using bits#, and want to present an overlay window to present any non-greyscale colours (e.g., red fixation point etc)
viewingDistance = 0.57; %in metres, please
screenNum = max(Screen('Screens')); %uses secondary monitor, if available
screenRect  = Screen(screenNum, 'rect');

%colour values
if usingBitsSharp %values when using bitsSharp need to be between 0-1, otherwise they need to be between 0-255 for psychtoolbox
    black = 0;
    grey = 0.5;
    red = [1 0 0];
else
    black = 0;
    grey = 128;
    red = [255 0 0];
end

% stimulus parameters
innerDiam.deg = 3; %in deg, please 
outerDiam.deg = 15; %in deg, please (was originally 256 in pixels)
sf.deg = 1; %cycles/deg (Tomassini uses 6.9)
sigma.deg = 1.15; %Tomassini uses 0.072
max_contrast = 0.5;
cosineWindowSize.deg = 1; %0.35; %the size of the cosine window, i.e., where the stimulus varies from full contrast to no contrast at the inner and outer edges of the stimulus.
ori_surround = 15; %orientation (in deg, ± from vertical) of the surrounds from vertical

%timing parameters
stim_full_contrast_time = 0.3; %in sec
fade_time = 0.1; %in sec
total_stim_time = fade_time + stim_full_contrast_time + fade_time;
inter_stimulus_interval_time = 0.5; %the length of time for the gap between first and second interval
preStimulusOnsetTime = 0.5; %time (in sec) after a button press is made before the next stimulus is shown
initial_fixation_time = 1; %in sec

%fixation parameters
fix_size.pix = 5;
fixRect = CenterRect([0 0 fix_size.pix fix_size.pix],screenRect);
fix_size = 5;
fix_frames = 60;

%other stuff
buttonAlreadyPressed = true; %don't accept any responses yet
response_labels = {'Staircase Num','Orientation mean','Min entropy value','Were intervals flipped?','Ppt Response','Reaction Time'};
stim_labels = {'Stim Diam (pix)','Spat freq','Sigma','Duration'}; 
int_labels = {'Staircase num','Ctr orient 1','Sur orient 1','Ctr orient 2','Sur orient 2','Ctr phase 1','Sur phase 1','Ctr phase 2','Sur phase 2'};
numIntervals = 2; %it's a 2 alternative forced choice task, so always leave this at 2. This is just for easier reading of code later on.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert angular sizes into physical size or pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert outer annulus from visual degrees to pixels
[outerDiam.m] = CalcPhysicalFromAngularSize(outerDiam.deg, viewingDistance);
[outerDiam.pix] = round(CalcPixelsFromPhysical(outerDiam.m)); %This process converts from angular to physical to pixels

%Convert inner aperture from visual degrees to pixels
[innerDiam.m] = CalcPhysicalFromAngularSize(innerDiam.deg, viewingDistance);
[innerDiam.pix] = round(CalcPixelsFromPhysical(innerDiam.m)); %This process converts from angular to physical to pixels

%Convert the gabor sigma from visual degrees to pixels
[sigma.m] = CalcPhysicalFromAngularSize(sigma.deg, viewingDistance);
[sigma.pix] = round(CalcPixelsFromPhysical(sigma.m)); %This process converts from angular to physical to pixels

%Convert the cosine window (where the stimulus varies from full contrast to no contrast at the inner and outer edges of the stimulus) to pixels
[cosineWindowSize.m] = CalcPhysicalFromAngularSize(cosineWindowSize.deg, viewingDistance);
[cosineWindowSize.pix] = round(CalcPixelsFromPhysical(cosineWindowSize.m)); %This process converts from angular to physical to pixels

%Convert the gabor spatial frequency from cycles/visual degree to cycles/pixel
[sf.m] = CalcPhysicalFromAngularSize(1/sf.deg, viewingDistance); %1/sf as the original unit is cycles/deg (not degrees/cycle)
[sf.pix] = round(CalcPixelsFromPhysical(sf.m)); %This process converts from angular to physical to pixels
sf.m = 1/sf.m; %convert it back to cycles/pix (not pix/cycles).
sf.pix = 1/sf.pix; %convert it back to cycles/pix (not pix/cycles).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Set psi (staircase) parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trialsPerStaircase = 1; %how many trials we will run per adaptive psi procedure
numStaircases = 2; %the amount of staircases we are running in this session

miss_rate = 0.04; %this is the recommended value in the Kontsevich & Tyler 1999 paper.

%intensity values (in this case, possible orientation differences)
xlevels = 281;
xstep = 0.1;
xbias = 0;

%lambda = (a,b) where a = threshold and b = slope of the psychometric function
%different possibilities of threshold values to examine
alevels = 41; %15; %test this many threshold levels
astep = 0.5; %1.0;  %with this difference between each threshold values
abias = 0;    %predicted threshold value to be examined around

%different slope parameters to compute
blevels = 9; %compare values for this many slope values
bratio = sqrt(2); %discriminability (d' = z*sqrt(2)) - the signal to noise ratio of the decision variable when a choice is made
bbase = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Initialize psi parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%creates a range of values on either side of the predicted bias
xmin = xbias - xstep*(xlevels-1)/2;
xmax = xbias + xstep*(xlevels-1)/2;
x = xmin:xstep:xmax; %increment between xmin and xmax given a specified step size

%creates a range of values on either side of the predicted bias
amin = abias - astep*(alevels-1)/2;
amax = abias + astep*(alevels-1)/2;
a = amin:astep:amax; %increment between amin and amax given a specified step size

logbmin = log(bbase) - log(bratio)*(blevels-1)/2;
logbmax = log(bbase) + log(bratio)*(blevels-1)/2;

b = zeros(1,blevels);
for i=1:blevels
    b(i) = exp(logbmin)*(bratio)^(i-1); %exponential value times by root 2 a tonne (n-1) of times 
end

%pre-allocate matrices for improved speed
pts_x = zeros(size(x)); % prob of success after next trial as fn of stim, x
psuccess_lx = zeros(alevels,blevels,xlevels); %probability of success given lambda and x
ptl_xsuccess = zeros(alevels,blevels,xlevels,numStaircases); %create an array of the probabilities for each possible threshold, slope, intensity value and for each staircase we're running
ptl_xfailure = zeros(alevels,blevels,xlevels,numStaircases); %create an array of the probabilities for each possible threshold, slope, intensity value and for each staircase we're running
minEntropy_idx = zeros(1,numStaircases); %small array to give the index for the value with the minimum entropy for each one of our staircases

%initialize output statistics for each staircase
a_hat = zeros(1,numStaircases);
a_err = zeros(1,numStaircases);
b_hat = zeros(1,numStaircases);
b_err = zeros(1,numStaircases);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Some initial calculations for the psi world
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up values for the prior
ptl = ones(alevels,blevels,numStaircases); %probability of t as a function of lambda (ptl)
ptl = ptl./(alevels*blevels); %uniform prior distribution (every value for each staircase (e.g., different thresholds and slopes) is equally likely)

% set up look-up table of conditional probabilities
for i = 1:alevels
    for j = 1:blevels
        %a 3D array that computes the probability of each stimulus value being considered as correct, given the intensity, threshold and slope of the psychometric function
        psuccess_lx(i,j,:) = logistic3([a(i) b(j) miss_rate],x); 
    end
end

%create order for staircases to be shown
totalTrials = trialsPerStaircase*numStaircases; %calculate number of stimuli per condition (num staircases x num trials per staircase)
responseArray = zeros(totalTrials,length(response_labels)); %create a blank response array for each of 4 different tidbits of data we want to keep track of
int_params = zeros(totalTrials,9); %save orientation and phase details for each trial and interval
psi_order = Shuffle(ceil([1:totalTrials]./trialsPerStaircase)); %randomize the order of presentation between the two staircases


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Open screen and set gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if usingBitsSharp
% 
%     changeMode('mono'); %ensure that we are in the correct mode.
%     PsychImaging('PrepareConfiguration'); %Initialize for Mono++ setting
%     PsychImaging('AddTask','General','FloatPoint32Bit'); %Initialize for Mono++ setting
% 
%     if overlay
%         PsychImaging('AddTask','General','EnableBits++Mono++OutputWithOverlay'); %Enable Mono++ drawing with an overlay window
%     else
%         PsychImaging('AddTask','General','EnableBits++Mono++Output'); %Enable Mono++ drawing but using greyscale only (no overlay window)
%     end
%     
%     %Apply a gamma correction to the screen
%     applyGammaCorrection('batmanLUT_20140430'); %can also use 'batman_20131015' or 'linear' setting
%     
%     %open display window(s) - we cannot use Screen('OpenWindow',...) command in this case
%     window = PsychImaging('OpenWindow', screenNum, grey, [], [], [], [], 0);
%     if overlay %open overlay window if necessary
%        overlaywindow = PsychImaging('GetOverlayWindow',window);
% 
%        %only matters if we are using bits# with an overlay window
%        cLUT_red = linspace(0,1,256); %equally spaced values to vary between black and yellow
%        cLUT_green = linspace(0,1,256); %equally spaced values to vary between black and yellow
%        cLUT_blue = zeros(1,256); %equally spaced values to vary between black and yellow
%        cLUT = [cLUT_red; cLUT_green; cLUT_blue]';
%        
%        %Load T-lock cLUT at next flip
%        Screen('LoadNormalizedGammaTable',window,cLUT,2); %the '2' indicates that the cLUT should be loaded to bits# and not the graphics card. The T-Lock is done automatically.
%     end
% 
% else
    
    %open display window - grey background
    window = Screen('OpenWindow', screenNum, grey);
    
%     %Gamma correction for Batman testing computer (latest as of 30/4/2014 - using the new ColorCal MKII)
%     R_gamma = 2.6678;
%     G_gamma = 2.6592;
%     B_gamma = 2.6497;
%     
%     % Create LUT
%     %lut_x = [0:1/255:1]';
%     %lut = [lut_x.^(1/R_gamma) lut_x.^(1/G_gamma) lut_x.^(1/B_gamma)];
%     %Screen(window,'LoadNormalizedGammaTable',lut);
% end

[I, ~, alpha] = imread(['C:\Kieran\honoursweirdstuff\dudepractice\1a.bmp']);
pracst1 = Screen('MakeTexture',window,I);
[I, ~, alpha] = imread(['C:\Kieran\honoursweirdstuff\dudepractice\1b.bmp']);
pracst2 = Screen('MakeTexture',window,I);
[I, ~, alpha] = imread(['C:\Kieran\honoursweirdstuff\dudepractice\2a.bmp']);
pracst3 = Screen('MakeTexture',window,I);
[I, ~, alpha] = imread(['C:\Kieran\honoursweirdstuff\dudepractice\2b.bmp']);
pracst4 = Screen('MakeTexture',window,I);

[screenWidth, screenHeight]=Screen('WindowSize', window);
SCREEN_X=screenWidth;
SCREEN_Y=screenHeight;

% texture containing grey screen
blank(1) = Screen('MakeTexture',window,127.*ones(SCREEN_Y, SCREEN_X));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Create 'Loading' Screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%draw 'loading' text
Screen('TextFont',window,'Arial');
Screen('TextSize',window,36);
%DrawFormattedText(window, 'Practice Trials', 'center', 'center');

%Screen('Flip',window,[],[],[],1); %display on screen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Do some initial calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scr.Hz = Screen('NominalFrameRate',screenNum); %get frame rate of the monitor
scr.frameLength = 1/scr.Hz; %how long each individual frame takes

%initialize variables for the event log
eventCounter = 1; %start writing to the event log on the first line
eventLog = {}; %start with empty cell array

responseStartTime = 0;  %Have to initialize this for the first trial

%pre-allocate stimulus matrices for faster processing
centreStim = zeros(outerDiam.pix, outerDiam.pix, numIntervals);
centreGrating = zeros(outerDiam.pix, outerDiam.pix, numIntervals);
surroundStim = zeros(outerDiam.pix, outerDiam.pix, numIntervals);
surroundGrating = zeros(outerDiam.pix, outerDiam.pix, numIntervals);

%create matrices which define locations for centre and surround apertures
centreAperture = createAperture(round(innerDiam.pix/2), outerDiam.pix, cosineWindowSize.pix, 0);
surroundAperture = createAperture(round(outerDiam.pix/2), outerDiam.pix, cosineWindowSize.pix, 1);
[ridx, cidx] = find(~isnan(centreAperture)); %find all non-NaN entries

%Initialise temporal window for stim fade on/off
fadeNumFrames = fade_time / scr.frameLength;
samples = (0:scr.frameLength:(fade_time-scr.frameLength)); %get equally space points between start time and time of duration (one contrast point for each frame that is presented)
temporalWindow_fadeOut = ((cos((pi*samples)/fade_time)*max_contrast)+max_contrast)/2; %create appropriate contrasts for each frame
temporalWindow_fadeIn = fliplr(temporalWindow_fadeOut);

%contrast 2
stim_durn_frames = 60; 
pdc = 1;
cont = pdc.*ones(1,stim_durn_frames); % initialize to peak dot contrast (defined above)
win_length = stim_durn_frames/2;     % define window length
cont(1:win_length) = pdc.*0.5.*(1-cos(pi*(1:win_length)./win_length)); % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length); % ... and down


%define surround orientations for each interval (1 left, 1 right of vertical)
bothSurrounds = [-ori_surround ori_surround];
bothMeans = [oriMean -oriMean] + subjVert; %compute experiment for both positive and negative mean values

disp(' ');
disp('Key experiment parameters:');
disp(['Orientation Difference: ' num2str(oriMean)]);
disp(['Subjective Vertical: ' num2str(subjVert)]);
disp(['Mean Orientations: ' num2str(bothMeans)]);
disp(' ');
pause(3);

delay_frames = 60; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off;
seed = sum(100*clock); %save the seed so we can recreate stimuli exactly at a later point in time (if necessary)
%rand('state',randomSeed); %apply random seed
randomSeed = rng(seed,'twister'); %rng('shuffle');
KbCheck; % run to load in mem
GetSecs; %takes longer on first call then on subsequent calls
Screen('Blendfunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %not entirely sure how this works - but it allows the 'globalAlpha' parameter in drawTexture to have control of the alpha for the whole texture
HideCursor; %hides the cursor during stimulus presentation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Create 'Ready' Screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%draw ready screen with instructions text
%Screen('TextFont',window,'Arial');
Screen('TextSize',window,26);
% Screen('TextStyle',window,1);
% DrawFormattedText(window,'BLOCK 3','center',350,black,[],[],[],1.65);

Screen('TextSize',window,20);
Screen('TextStyle',window,0);
DrawFormattedText(window,'On each trial, indicate whether the first or second interval\ncontains a grating that is oriented closer to vertical\n\nPress F for the 1st interval; press J for the 2nd interval.',...
                         'center',420,black,[],[],[],1.65);
Screen('Flip',window,[],[],[],1); %display on screen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Code to initiate experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startExp = false;
while (startExp == false)
    keypress = 0;
    while ~keypress
        [keypress, ~,~] = KbCheck;
        startExp = true;
    end
    %check for `board responses to start experiment
%     [keyIsDown,~,keyCode] = KbCheck;
%     if keyIsDown
%         if (keyCode(ResponseLeft | ResponseRight)) %if participant presses 's' or space bar, break loop and start experiment
%             startExp = true;
%         elseif keyCode(ResponseExit) % quit events - Q key or ESC
%             exitGracefully('User has aborted script with keyboard press');
%         end
%     end
end

the_experiment_start = GetSecs;
newEvent('Experiment Start',[],[]); %to keep a log of timing of the experiment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Initial fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show fixation only for initial fixation
drawFixationAndFlip(window,red,fixRect);

targetTime = GetSecs + initial_fixation_time;
while GetSecs < targetTime
    %run up the clock
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Trial Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:totalTrials
    
    %randomly select which psi (staircase) will be presented this trial
    psi = psi_order(t); 
    
    % 1. calculate conditional probability of response, r, after presenting stimulus intensity, x.
    for k = 1:xlevels %check the probability of each intensity level
        pts_x(k) = sum(sum(ptl(:,:,psi).*psuccess_lx(:,:,k))); %prior (uniform) * lookup table of likely probabilities for each value - then summed across each psychometric function
    end
    
    % 2. use Bayes rule to estimate posterior prob of each psychometric function
    for k = 1:xlevels
        ptl_xsuccess(:,:,k,psi) = (ptl(:,:,psi).*psuccess_lx(:,:,k))./pts_x(k);
        ptl_xfailure(:,:,k,psi) = (ptl(:,:,psi).*(1-psuccess_lx(:,:,k)))./(1-pts_x(k));
    end    
    
    % 3. estimate entropy of pdf as a fn of stim level and response
    for k = 1:xlevels
        HS(k) = -sum(sum(ptl_xsuccess(:,:,k,psi).*log(ptl_xsuccess(:,:,k,psi)))); %entropy if response is successful (for each intensity value)
        HF(k) = -sum(sum(ptl_xfailure(:,:,k,psi).*log(ptl_xfailure(:,:,k,psi)))); %entropy if response is incorrect (for each intensity value)
    end
    
    % 4. estimate expected entropy for each stim level
    for k = 1:xlevels
        EH(k) = HS(k).*pts_x(k) + HF(k).*(1-pts_x(k));
    end
    
    % 5. find stim level with minimum entropy
    [minEntropy_val, minEntropy_idx(psi)] = min(EH);
    
    % 6. Generate a grating with the centre stimulus at the minimum entropy level (as identified in step 5)

    imgArray = zeros(outerDiam.pix, outerDiam.pix, numIntervals); %reset stim matrix to be zeros before each trial
    
    if bothMeans(psi) > subjVert
        if firstorsecond == 1;
            centreOrient = [0 10]; %[bothMeans(psi)-(x(minEntropy_idx(psi))) bothMeans(psi)+(x(minEntropy_idx(psi)))]; %calculate orientation for centre stimulus - smallest value first
        else 
            centreOrient = [20 -5];
        end
    else
        if firstorsecond == 1;
            centreOrient = [0 10]; %[bothMeans(psi)-(x(minEntropy_idx(psi))) bothMeans(psi)+(x(minEntropy_idx(psi)))]; %calculate orientation for centre stimulus - smallest value first
        else 
            centreOrient = [20 -5];
        end
    end
    if relativeSurr
        surroundOrient = bothSurrounds + centreOrient; %get the order of the surround orientation for this staircase - smallest value first (so illusion is always in the direction of the mean)
    else
        surroundOrient = bothSurrounds; %get the order of the surround orientation for this staircase - smallest value first (so illusion is always in the direction of the mean)
    end
    centrePhase = [rand rand]; %randomize phase for each interval for centre stimulus
    surroundPhase = [rand rand]; %%randomize phase for each interval for surround stimulus
    
    %on randomly half the trials, flip order of presentation so leftward surround isn't always presented first
    flipIntervals = round(rand);
    if flipIntervals
        centreOrient = centreOrient;
        surroundOrient = surroundOrient;
    end
     
    for int=1:numIntervals %2-alternative forced choice task, thus need two stimuli

        %create centre grating
        centreStim(:,:,int) = MakeGrating(sf.pix,centreOrient(int),max_contrast,outerDiam.pix,centrePhase(int)); %make oriented grating to be used for the centre of stimulus
        centreGrating(:,:,int) = (centreStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to the central aperture only
        
        %create surround grating
        surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),0,outerDiam.pix,surroundPhase(int)); %NO SURROUND - swap this in for the line below to remove surround
        %surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),max_contrast,outerDiam.pix,surroundPhase(int)); %make oriented grating for stimulus surround
        surroundGrating(:,:,int) = (surroundStim(:,:,int).*surroundAperture) + ((1-surroundAperture).*128); %limit grating to within the outer edges of the surround aperture
        
        %combine gratings
        imgArray(:,:,int) = surroundGrating(:,:,int);
        for i=1:length(ridx) %for all values where the central grating exists
            imgArray(ridx(i),cidx(i),int) = centreGrating(ridx(i),cidx(i),int); %replace values of the surround stimulus with ones from the centre
        end
    end
%     for int=1:numIntervals %2-alternative forced choice task, thus need two stimuli
% 
%         %create centre grating
%         centreStim(:,:,int) = MakeGrating(sf.pix,centreOrient(int),max_contrast,outerDiam.pix,centrePhase(int)); %make oriented grating to be used for the centre of stimulus
%         centreGrating(:,:,int) = (centreStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to the central aperture only
%         
%         %create surround grating
%         %surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),max_contrast,outerDiam.pix,surroundPhase(int)); %make oriented grating for stimulus surround
%         surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),0,outerDiam.pix,surroundPhase(int)); %NO SURROUND - swap this in for the line below to remove surround
%         surroundGrating(:,:,int) = (surroundStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to within the outer edges of the surround aperture
%         
%         %combine gratings
%         imgArray(:,:,int) = surroundGrating(:,:,int);
%         for i=1:length(ridx) %for all values where the central grating exists
%             imgArray(ridx(i),cidx(i),int) = centreGrating(ridx(i),cidx(i),int); %replace values of the surround stimulus with ones from the centre
%         end
%     end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Present stimulus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for fix=1:fix_frames
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('DrawDots', window, [(SCREEN_X-fix_size)/2 (SCREEN_Y-fix_size)/2], fix_size, [0, 0, 0],[0,0],1);
        Screen('Flip',window);
    end
    for dd=1:30
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('Flip',window);
    end
    
    for int=1:numIntervals

        if int==1
            targetTime = GetSecs; %start timer now
        elseif int==2 %only start accepting responses/counting reaction time when 2nd stimulus is displayed
            responseStartTime = GetSecs;  %Reaction time measured from onset of test image
            buttonAlreadyPressed = false; %start accepting responses at this point
        end
        
        if usingBitsSharp
            imgArray(:,:,int) = imgArray(:,:,int) / 255; %everything needs to be between 0 and 1 for bits#
            imageTexture = Screen('MakeTexture',window,imgArray(:,:,int),[],[],2); %The 2 sends it to the bits# instead of the graphics card
        else
            imageTexture = Screen('MakeTexture',window,imgArray(:,:,int)); %turn matrix into useful texture that can be displayed on screen
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fading in

        startFadeTime = targetTime;
        targetTime = targetTime + fade_time; %specify when the trial presentation should end
        
        newEvent(sprintf('Begin interval %d',int),centreOrient(int),surroundOrient(int)); %record time of stim presentation, along with centre and surround orientations
        
        while GetSecs < targetTime %during stimulus presentation
            currentFrame = ceil((GetSecs - startFadeTime)/scr.frameLength); %find appropriate contrast for how long we are into the fading time
            if (currentFrame <= fadeNumFrames) %just to ensure it doesn't go over frames and cause an error.... just in case!
                Screen('DrawTexture',window,imageTexture,[],[],[],[],temporalWindow_fadeIn(currentFrame)); %draw all gabor patches, fading in
                drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
                handleKeyPress(responseStartTime,flipIntervals);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %full contrast
        
        Screen('DrawTexture',window,imageTexture); %draw stimulus image
        drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
        
        %forced stimulus presentation time
        targetTime = targetTime + stim_full_contrast_time;
        while GetSecs < targetTime
            handleKeyPress(responseStartTime,flipIntervals);
        end
        
        %GetImage call. Alter the rect argument to change the location of the screen shot
        screenshotMatrix = Screen('GetImage', window, [screenRect(3)/2-size(imgArray(:,:,1),1)/2 screenRect(4)/2-size(imgArray(:,:,1),2)/2 screenRect(3)/2+size(imgArray(:,:,1),1)/2 screenRect(4)/2+size(imgArray(:,:,1),2)/2]);
        imwrite(screenshotMatrix, 'new_interval_1.jpg')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fading out
        
        startFadeTime = targetTime; 
        targetTime = targetTime + fade_time; %specify when the trial presentation should end
        
        while GetSecs < targetTime %during stimulus presentation
            currentFrame = ceil((GetSecs - startFadeTime)/scr.frameLength); %find appropriate contrast for how long we are into the fading time
            if (currentFrame <= fadeNumFrames) %just to ensure it doesn't go over frames and cause an error.... just in case!
                Screen('DrawTexture',window,imageTexture,[],[],[],[],temporalWindow_fadeOut(currentFrame)); %draw all gabor patches, fading in
                drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
                handleKeyPress(responseStartTime,flipIntervals); %handle any responses/keypresses
            end
        end
        
        newEvent('Stim off',[],GetSecs - the_experiment_start - eventLog{eventCounter-1,1}); %get length of stimulus duration
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Inter-stimulus interval (1st int only)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Show blank screen after test (inter-stimulus interval)
        drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
        
        if int==1
            targetTime = targetTime + inter_stimulus_interval_time;

            %run up the clock while checking for key presses
            while GetSecs < targetTime
                handleKeyPress(0,flipIntervals);
            end
        end
        
    end %each interval

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Inter-trial interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %wait until user responds
    while (~buttonAlreadyPressed)
        handleKeyPress(responseStartTime,flipIntervals);
    end
    
    % 7. keep the posterior pdf that corresponds to the completed trial
    
    if pickedLeftStim==1 %if the participant responded to the interval corresponding to the stim with leftward surround
        ptl(:,:,psi) = ptl_xsuccess(:,:,minEntropy_idx(psi),psi); %we are correct, then update our probability function accordingly (for use in the next trial)
        
    else %if the participant responded to the interval corresponding to the stim with rightward surround
        ptl(:,:,psi) = ptl_xfailure(:,:,minEntropy_idx(psi),psi); %so update probability that we were incorrect this time
    end
    
    %save stimulus and response details
    responseArray(t,:) = [psi, bothMeans(psi), x(minEntropy_idx(psi)), flipIntervals, response, responseTime];
    int_params(t,:) = [psi centreOrient(1) surroundOrient(1) centreOrient(2) surroundOrient(2) centrePhase(1) surroundPhase(1) centrePhase(2) surroundPhase(2)];
    
    %forced gap after participant response and before start of next trial
    targetTime = GetSecs + preStimulusOnsetTime;
    while GetSecs < targetTime
    end
    firstorsecond = firstorsecond + 1
    
    if firstorsecond == 2;
        keypress = 0;
        while ~keypress
            [keypress, keysecs,keyCode] = KbCheck;
            if response == 1;
                DrawFormattedText(window, 'Correct Response \n \n The first grating was oriented closer to vertical, so the F key should have been pressed, ', 'center', 'center');
            else
                DrawFormattedText(window, 'Incorrect Response \n \n The first grating was oriented closer to vertical, so the F key should have been pressed', 'center', 'center');
            end
            Screen('Flip',window,[],[],[],1);
        end
    elseif firstorsecond == 3; 
        keypress = 0;
        while ~keypress
            [keypress, keysecs,keyCode] = KbCheck;
            if response == 2;
                DrawFormattedText(window, 'Correct Response \n \n The second grating was oriented closer to vertical, so the J key should have been pressed', 'center', 'center');
            else
                DrawFormattedText(window, 'Incorrect Response \n \n The second grating was oriented closer to vertical, so the J key should have been pressed', 'center', 'center');
            end
            Screen('Flip',window,[],[],[],1);
        end
        
    end
    for dd=1:delay_frames
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('Flip',window);
    end
        

end % end of stuff done every trial

for jj=1:stim_durn_frames
    Screen('TextSize',window,20);
    Screen('TextStyle',window,0);
    DrawFormattedText(window,'On each trial, indicate whether the first or second interval\ncontains a face that is oriented more directly at you\n\nPress F for the 1st interval; press J for the 2nd interval.',...
        'center',420,black,[],[],[],1.65);
    Screen('Flip',window,[],[],[],1); %display on screen
end

keypress = 0;
while ~keypress
    [keypress, ~,~] = KbCheck;
end
for dd=1:delay_frames
    Screen('DrawTexture', window, blank(1), [], [], [], 0);
    Screen('Flip',window);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Trial Loop 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:totalTrials
    
    %randomly select which psi (staircase) will be presented this trial
    psi = psi_order(t); 
    
    % 1. calculate conditional probability of response, r, after presenting stimulus intensity, x.
    for k = 1:xlevels %check the probability of each intensity level
        pts_x(k) = sum(sum(ptl(:,:,psi).*psuccess_lx(:,:,k))); %prior (uniform) * lookup table of likely probabilities for each value - then summed across each psychometric function
    end
    
    % 2. use Bayes rule to estimate posterior prob of each psychometric function
    for k = 1:xlevels
        ptl_xsuccess(:,:,k,psi) = (ptl(:,:,psi).*psuccess_lx(:,:,k))./pts_x(k);
        ptl_xfailure(:,:,k,psi) = (ptl(:,:,psi).*(1-psuccess_lx(:,:,k)))./(1-pts_x(k));
    end    
    
    % 3. estimate entropy of pdf as a fn of stim level and response
    for k = 1:xlevels
        HS(k) = -sum(sum(ptl_xsuccess(:,:,k,psi).*log(ptl_xsuccess(:,:,k,psi)))); %entropy if response is successful (for each intensity value)
        HF(k) = -sum(sum(ptl_xfailure(:,:,k,psi).*log(ptl_xfailure(:,:,k,psi)))); %entropy if response is incorrect (for each intensity value)
    end
    
    % 4. estimate expected entropy for each stim level
    for k = 1:xlevels
        EH(k) = HS(k).*pts_x(k) + HF(k).*(1-pts_x(k));
    end
    
    % 5. find stim level with minimum entropy
    [minEntropy_val, minEntropy_idx(psi)] = min(EH);
    
    % 6. Generate a grating with the centre stimulus at the minimum entropy level (as identified in step 5)

    imgArray = zeros(outerDiam.pix, outerDiam.pix, numIntervals); %reset stim matrix to be zeros before each trial
    
    if bothMeans(psi) > subjVert
        if firstorsecond == 1;
            centreOrient = [0 10]; %[bothMeans(psi)-(x(minEntropy_idx(psi))) bothMeans(psi)+(x(minEntropy_idx(psi)))]; %calculate orientation for centre stimulus - smallest value first
        else 
            centreOrient = [20 -5];
        end
    else
        if firstorsecond == 1;
            centreOrient = [0 10]; %[bothMeans(psi)-(x(minEntropy_idx(psi))) bothMeans(psi)+(x(minEntropy_idx(psi)))]; %calculate orientation for centre stimulus - smallest value first
        else 
            centreOrient = [20 -5];
        end
    end
    if relativeSurr
        surroundOrient = bothSurrounds + centreOrient; %get the order of the surround orientation for this staircase - smallest value first (so illusion is always in the direction of the mean)
    else
        surroundOrient = bothSurrounds; %get the order of the surround orientation for this staircase - smallest value first (so illusion is always in the direction of the mean)
    end
    centrePhase = [rand rand]; %randomize phase for each interval for centre stimulus
    surroundPhase = [rand rand]; %%randomize phase for each interval for surround stimulus
    
    %on randomly half the trials, flip order of presentation so leftward surround isn't always presented first
    flipIntervals = round(rand);
    if flipIntervals
        centreOrient = centreOrient;
        surroundOrient = surroundOrient;
    end
     
    for int=1:numIntervals %2-alternative forced choice task, thus need two stimuli

        %create centre grating
        centreStim(:,:,int) = MakeGrating(sf.pix,centreOrient(int),max_contrast,outerDiam.pix,centrePhase(int)); %make oriented grating to be used for the centre of stimulus
        centreGrating(:,:,int) = (centreStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to the central aperture only
        
        %create surround grating
        surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),0,outerDiam.pix,surroundPhase(int)); %NO SURROUND - swap this in for the line below to remove surround
        %surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),max_contrast,outerDiam.pix,surroundPhase(int)); %make oriented grating for stimulus surround
        surroundGrating(:,:,int) = (surroundStim(:,:,int).*surroundAperture) + ((1-surroundAperture).*128); %limit grating to within the outer edges of the surround aperture
        
        %combine gratings
        imgArray(:,:,int) = surroundGrating(:,:,int);
        for i=1:length(ridx) %for all values where the central grating exists
            imgArray(ridx(i),cidx(i),int) = centreGrating(ridx(i),cidx(i),int); %replace values of the surround stimulus with ones from the centre
        end
    end
%     for int=1:numIntervals %2-alternative forced choice task, thus need two stimuli
% 
%         %create centre grating
%         centreStim(:,:,int) = MakeGrating(sf.pix,centreOrient(int),max_contrast,outerDiam.pix,centrePhase(int)); %make oriented grating to be used for the centre of stimulus
%         centreGrating(:,:,int) = (centreStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to the central aperture only
%         
%         %create surround grating
%         %surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),max_contrast,outerDiam.pix,surroundPhase(int)); %make oriented grating for stimulus surround
%         surroundStim(:,:,int) = MakeGrating(sf.pix,surroundOrient(int),0,outerDiam.pix,surroundPhase(int)); %NO SURROUND - swap this in for the line below to remove surround
%         surroundGrating(:,:,int) = (surroundStim(:,:,int).*centreAperture) + ((1-centreAperture).*128); %limit grating to within the outer edges of the surround aperture
%         
%         %combine gratings
%         imgArray(:,:,int) = surroundGrating(:,:,int);
%         for i=1:length(ridx) %for all values where the central grating exists
%             imgArray(ridx(i),cidx(i),int) = centreGrating(ridx(i),cidx(i),int); %replace values of the surround stimulus with ones from the centre
%         end
%     end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Present stimulus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for fix=1:fix_frames
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('DrawDots', window, [(SCREEN_X-fix_size)/2 (SCREEN_Y-fix_size)/2], fix_size, [0, 0, 0],[0,0],1);
        Screen('Flip',window);
    end
    for dd=1:30
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('Flip',window);
    end
    
    for int=1:numIntervals
        if int==1
            targetTime = GetSecs; %start timer now
        elseif int==2 %only start accepting responses/counting reaction time when 2nd stimulus is displayed
            responseStartTime = GetSecs;  %Reaction time measured from onset of test image
            buttonAlreadyPressed = false; %start accepting responses at this point
        end
        
        if usingBitsSharp
            imgArray(:,:,int) = imgArray(:,:,int) / 255; %everything needs to be between 0 and 1 for bits#
            imageTexture = Screen('MakeTexture',window,imgArray(:,:,int),[],[],2); %The 2 sends it to the bits# instead of the graphics card
        else
            imageTexture = Screen('MakeTexture',window,imgArray(:,:,int)); %turn matrix into useful texture that can be displayed on screen
        end
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %fading in
% 
%         startFadeTime = targetTime;
%         targetTime = targetTime + fade_time; %specify when the trial presentation should end
%         
%         newEvent(sprintf('Begin interval %d',int),centreOrient(int),surroundOrient(int)); %record time of stim presentation, along with centre and surround orientations
%         
%         while GetSecs < targetTime %during stimulus presentation
%             currentFrame = ceil((GetSecs - startFadeTime)/scr.frameLength); %find appropriate contrast for how long we are into the fading time
%             if (currentFrame <= fadeNumFrames) %just to ensure it doesn't go over frames and cause an error.... just in case!
%                 Screen('DrawTexture',window,imageTexture,[],[],[],[],temporalWindow_fadeIn(currentFrame)); %draw all gabor patches, fading in
%                 drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
%                 handleKeyPress(responseStartTime,flipIntervals);
%             end
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %full contrast
%         
%         Screen('DrawTexture',window,imageTexture); %draw stimulus image
%         drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
%         
%         %forced stimulus presentation time
%         targetTime = targetTime + stim_full_contrast_time;
%         while GetSecs < targetTime
%             handleKeyPress(responseStartTime,flipIntervals);
%         end
%         
%         %GetImage call. Alter the rect argument to change the location of the screen shot
%         screenshotMatrix = Screen('GetImage', window, [screenRect(3)/2-size(imgArray(:,:,1),1)/2 screenRect(4)/2-size(imgArray(:,:,1),2)/2 screenRect(3)/2+size(imgArray(:,:,1),1)/2 screenRect(4)/2+size(imgArray(:,:,1),2)/2]);
%         imwrite(screenshotMatrix, 'new_interval_1.jpg')
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %fading out
%         
%         startFadeTime = targetTime; 
%         targetTime = targetTime + fade_time; %specify when the trial presentation should end
%         
%         while GetSecs < targetTime %during stimulus presentation
%             currentFrame = ceil((GetSecs - startFadeTime)/scr.frameLength); %find appropriate contrast for how long we are into the fading time
%             if (currentFrame <= fadeNumFrames) %just to ensure it doesn't go over frames and cause an error.... just in case!
%                 Screen('DrawTexture',window,imageTexture,[],[],[],[],temporalWindow_fadeOut(currentFrame)); %draw all gabor patches, fading in
%                 drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
%                 handleKeyPress(responseStartTime,flipIntervals); %handle any responses/keypresses
%             end
%         end
%         
%         newEvent('Stim off',[],GetSecs - the_experiment_start - eventLog{eventCounter-1,1}); %get length of stimulus duration
        
        if fristorsecodn == 1
            for jj=1:stim_durn_frames
                Screen('DrawTexture', window, pracst1, [], [], [], 0, cont(jj));
                Screen('Flip',window,[],[],[],1);
            end
        elseif fristorsecodn == 2
            for kk=1:stim_durn_frames
                Screen('DrawTexture', window, pracst2, [], [], [], 0, cont(kk));
                Screen('Flip',window,[],[],[],1);
            end
        elseif fristorsecodn == 3
            for jj=1:stim_durn_frames
                Screen('DrawTexture', window, pracst3, [], [], [], 0, cont(jj));
                Screen('Flip',window,[],[],[],1);
            end
        elseif fristorsecodn == 4
            for kk=1:stim_durn_frames
                Screen('DrawTexture', window, pracst4, [], [], [], 0, cont(kk));
                Screen('Flip',window,[],[],[],1);
            end
        end
        
        
        fristorsecodn = fristorsecodn + 1
        
        for dd=1:30
            Screen('DrawTexture', window, blank(1), [], [], [], 0);
            Screen('Flip',window);
        end
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Inter-stimulus interval (1st int only)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Show blank screen after test (inter-stimulus interval)
        drawFixationAndFlip(window,red,fixRect); %draw fixation and show current texture on screen
        
        if int==1
            targetTime = targetTime + inter_stimulus_interval_time;

            %run up the clock while checking for key presses
            while GetSecs < targetTime
                handleKeyPress(0,flipIntervals);
            end
        end
        
    end %each interval

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Inter-trial interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %wait until user responds
    while (~buttonAlreadyPressed)
        handleKeyPress(responseStartTime,flipIntervals);
    end
    
    % 7. keep the posterior pdf that corresponds to the completed trial
    
    if pickedLeftStim==1 %if the participant responded to the interval corresponding to the stim with leftward surround
        ptl(:,:,psi) = ptl_xsuccess(:,:,minEntropy_idx(psi),psi); %we are correct, then update our probability function accordingly (for use in the next trial)
        
    else %if the participant responded to the interval corresponding to the stim with rightward surround
        ptl(:,:,psi) = ptl_xfailure(:,:,minEntropy_idx(psi),psi); %so update probability that we were incorrect this time
    end
    
    %save stimulus and response details
    responseArray(t,:) = [psi, bothMeans(psi), x(minEntropy_idx(psi)), flipIntervals, response, responseTime];
    int_params(t,:) = [psi centreOrient(1) surroundOrient(1) centreOrient(2) surroundOrient(2) centrePhase(1) surroundPhase(1) centrePhase(2) surroundPhase(2)];
    
    %forced gap after participant response and before start of next trial
    targetTime = GetSecs + preStimulusOnsetTime;
    while GetSecs < targetTime
    end
    firstorsecond = firstorsecond + 1
    
    keypress = 0;
    while ~keypress
        [keypress, keysecs,keyCode] = KbCheck;
        if response == 2;
            DrawFormattedText(window, 'Correct Response \n \n The second face was more direct, so the J key should have been pressed', 'center', 'center');
        else
            DrawFormattedText(window, 'Incorrect Response \n \n The second face was more direct, so the J key should have been pressed', 'center', 'center');
        end
        Screen('Flip',window,[],[],[],1);
    end
%     
%     if firstorsecond == 2;
%         keypress = 0;
%         while ~keypress
%             [keypress, keysecs,keyCode] = KbCheck;
%             if response == 1;
%                 DrawFormattedText(window, 'Correct Response \n \n The first grating was oriented closer to vertical, so the F key should have been pressed, ', 'center', 'center');
%             else
%                 DrawFormattedText(window, 'Incorrect Response \n \n The first grating was oriented closer to vertical, so the F key should have been pressed', 'center', 'center');
%             end
%             Screen('Flip',window,[],[],[],1);
%         end
%     elseif firstorsecond == 3; 
%         keypress = 0;
%         while ~keypress
%             [keypress, keysecs,keyCode] = KbCheck;
%             if response == 2;
%                 DrawFormattedText(window, 'Correct Response \n \n The second grating was oriented closer to vertical, so the F key should have been pressed', 'center', 'center');
%             else
%                 DrawFormattedText(window, 'Incorrect Response \n \n The second grating was oriented closer to vertical, so the F key should have been pressed', 'center', 'center');
%             end
%             Screen('Flip',window,[],[],[],1);
%         end
%         
%     end
    for dd=1:delay_frames
        Screen('DrawTexture', window, blank(1), [], [], [], 0);
        Screen('Flip',window);
    end
        

end % end of stuff done every trial

for jj=1:stim_durn_frames
    Screen('TextSize',window,20);
    Screen('TextStyle',window,0);
    DrawFormattedText(window,'The practice trials are now complete\nNote that the experimental trials will be more difficult, and feedback will not be given\n\nPlease call a researcher to proceed.',...
        'center',420,black,[],[],[],1.65);
    Screen('Flip',window,[],[],[],1); %display on screen
end

keypress = 0;
while ~keypress
    [keypress, ~,~] = KbCheck;
end


newEvent('Experiment End',[],[]); %to keep a log of timing of the experiment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Calculate results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 8. Find new estimate of psychometric function
for psi=1:numStaircases
    a_hat(psi) = sum(a*ptl(:,:,psi)); %threshold
    b_hat(psi) = sum(b*ptl(:,:,psi)'); %slope
    
    a_mean = ones(size(a)).*a_hat(psi);
    b_mean = ones(size(b)).*b_hat(psi);
    
    a_err(psi) = sqrt(sum(((a-a_mean).^2)*ptl(:,:,psi)));
    b_err(psi) = sqrt(sum(((b-b_mean).^2)*ptl(:,:,psi)'));     % standard error
end

% record stimulus details ...
stim_params = [outerDiam.pix, sf.pix, sigma.pix, total_stim_time];
psy_fun = [];
for z = 1:numStaircases
    psy_fun = [psy_fun; a_hat(z) a_err(z) b_hat(z) b_err(z) 0];
end

clear centreGrating surroundGrating centreStim surroundStim imgArray centreAperture surroundAperture; %clear some of the larger matrices to reduce results file size
%save(saveName);  %,[stim_params;responseArray;psy_fun],'\t');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Draw figure of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%analysis_psi_2int(subjID,acq,relativeSurr)

%plot_psi_results(responseArray,psy_fun);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Function down...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return to command window ...
newEvent('Experiment End',[],[]);
ShowCursor;
Screen('CloseAll');
%reset_gamma;
warning on;
clear all;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function handleKeyPress(responseStartTime, flipIntervals)
KbName('UnifyKeyNames');
ResponseLeft = KbName('F'); % arrows
ResponseRight = KbName('J'); %
ResponseExit = KbName('P');
RestrictKeysForKbCheck([ResponseLeft ResponseRight ResponseExit]);

global buttonAlreadyPressed;
global response;
global responseTime;
global pickedLeftStim;

%response keys are:
% 1 - either numpad 1 or 1! indicates leftward of vertical
% 2 - either numpad 2 or 2@ indicates rightward of vertical

[keyIsDown, keySecs, keyCode] = KbCheck; %check which buttons on the keyboard are being pressed

%check whether a key has been pressed - if not, exit function.
if ~keyIsDown
    return;

else %if a key HAS been pressed

    %if trying to abort the program - close it
    if keyCode(ResponseExit) % quit events - Q key or ESC
        exitGracefully('User has aborted script with keyboard press');
    end
    
    if (~buttonAlreadyPressed)
        if (keyCode(49)==1 && keyCode(50)==1) %if participant has somehow managed to press both buttons (1! & 2@) at once
            beep; %beep at them! (and don't accept it as a response)
        elseif  (keyCode(97)==1 && keyCode(98)==1) %repeat this for numpad keys (1end & 2downarrow)
            beep;
            
        else %otherwise record the participant response
            if keyCode(ResponseLeft) %if '1' (left) button is pressed
                response = 1;
                if ~flipIntervals %if there was no flipping, then the left stim was presented and we responded leftward stim was more vertical
                    pickedLeftStim = 1;
                else
                    pickedLeftStim = 0;
                end
                responseTime = GetSecs - responseStartTime; % time in msecs
                buttonAlreadyPressed = true; %stop accepting any more responses
                newEvent('User response',response,responseTime);
                
            elseif keyCode(ResponseRight) %if '2' (right) button is pressed
                response = 2;
                if ~flipIntervals %if there was no flipping, then the left stim was presented and we responded leftward stim was more vertical
                    pickedLeftStim = 0;
                else
                    pickedLeftStim = 1;
                end
                responseTime = GetSecs - responseStartTime; % time in msecs
                buttonAlreadyPressed = true; %stop accepting any more responses
                newEvent('User response',response,responseTime);
            end
        end
    end

end % keyIsDown

end %handleKeyPress


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function newEvent(eventName,param1,param2)

    %keeps track of events during the experiment so we can check timing etc
    global eventCounter;
    global eventLog;
    global the_experiment_start;
    
    eventLog{eventCounter,1} = GetSecs - the_experiment_start;
    eventLog{eventCounter,2} = eventName;
    eventLog{eventCounter,3} = param1;
    eventLog{eventCounter,4} = param2;
    
    eventCounter = eventCounter + 1; %increase counter so next event is placed on the next line
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [raisedCosineMatrix] = createAperture(outerRadius, matrixSize, cosineWindowSize, surroundFlag)

if surroundFlag
    raisedCosineMatrix = zeros(matrixSize);
else
    raisedCosineMatrix = NaN(matrixSize); %repmat(999,matrixSize); %give an arbitrary unused value so we can separate stimulus entries to non-stimulus entries
end

[X,Y]  = meshgrid([-matrixSize/2+0.5:matrixSize/2-0.5],[-matrixSize/2+0.5:matrixSize/2-0.5]);
[th,r] = cart2pol(X,Y);

for col = 1:matrixSize %go through every entry in the matrix
    for row = 1:matrixSize
        
        if (r(row,col) <= outerRadius) %if this point is within the desired aperture
            
            if ((outerRadius - r(row,col)) <= cosineWindowSize) %and specifically the outer edge of the aperture
                if surroundFlag %if we are want to create a cosine window for this aperture
                    raisedCosineMatrix(row,col) = 1/2 * (1 - cos(2*pi*(outerRadius-r(row,col))/(2*cosineWindowSize))); %create cosine window at the outer edge
                else
                    raisedCosineMatrix(row,col) = 1; %otherwise put this part of the circle simply as 1 (and don't create a cosine window)
                end
            else %for all other parts inside the aperture
                raisedCosineMatrix(row,col) = 1; %define it as maximum contrast as per normal
            end
            
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function exitGracefully(errMsg)

global usingBitsSharp;

Screen('CloseAll'); %close open psychtoolbox textures
if usingBitsSharp
    applyGammaCorrection('LinearLUT');
else
    reset_gamma; %returns gamma to as it should be.
end
ShowCursor; %present the mouse pointer again
error(errMsg);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function drawFixationAndFlip(window, red, fixRect)

	%Draws fixation point and displays it on screen
   
	% Draw red round fixation with its centre at the middle of the screen
    %Screen('FillOval', window, red, fixRect);
    Screen('Flip',window);
    
end


%Small functions to convert to/from visual degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [size_m] = CalcPhysicalFromAngularSize(size_deg, viewingDistance) %answer in metres
	size_m = viewingDistance * tan(DegToRad(size_deg));
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [rad] = DegToRad(ang) %answer in radians
	rad = (ang / 360) * 2 * pi;
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [size_deg] = CalcAngularFromPhysicalSize(size_m, viewingDistance) %answer in degrees
	size_deg = RadToDeg(atan(size_m / viewingDistance));
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [deg] = RadToDeg(ang) %answer in degrees
	deg = ang / (2 * pi) * 360;
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [pix] = CalcPixelsFromPhysical(size_m)  
    %pixels per m.
    %pix = size_m * 3703.70; %BOLDscreen at scanner (screen size: 51.84cm x 32.4cm, resolution: 1920 x 1200)
    %pix = size_m * 4117.6; %My new lab computer (previously Tony's). Second smaller Dell LCD (screen size: 40.8cm x 25.5cm, resolution: 1680 x 1050);
    pix = size_m * 2844.44444444444; %BATMAN computer; left cubicle in room 519 (screen size: 36.0cm x 27.0cm, resolution: 1024 x 768);
    %pix = size_m * 3733.236548; %in metres. Based on measurements made using a ruler for my OLD work computer in room 519 (screen size 51.5cm x 32.1cm, resolution 1920 x 1200)
	return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [size_m] = CalcPhysicalFromPixels(pix)  
    %m per pixel
    %size_m = pix * 0.00027; %BOLDscreen at scanner (screen size: 51.84cm x 32.4cm, resolution: 1920 x 1200)
    %size_m = pix * 0.00024286; %My new lab computer (previously Tony's). Second smaller Dell LCD (screen size: 40.8cm x 25.5cm, resolution: 1680 x 1050);
    size_m = pix * 0.0003515625000549316; %BATMAN computer; left cubicle in room 519 (screen size: 36.0cm x 27.0cm, resolution: 1024 x 768);
    %size_m = pix * 0.0002678640871; %Based on measurements made using a ruler for my OLD work computer in room 519 (screen size 51.5cm x 32.1cm, resolution 1920 x 1200)
	return;
end