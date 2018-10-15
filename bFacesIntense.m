function bFacesIntense()
clear all

ID = 'tost2';
Condition = 'head' %'sizeDistance' %'presentationDelay' %'head';
% output filename
filename = ['C:\kieran\honoursweirdstuff\results\', Condition, ID, '.txt']; % CC

delay_frames = 120; %60 %120 %%for lab comp/timing condition
stim_durn_frames = 30; %60 %120
fix_frames = 60; %60 %120
num_images = 6;

fix_size = 5;
image_angle = [0 60 120 180 240 300]; %0:60:300; % degrees
image_dist = [90 140 140 90 140 140]; %[90 140 140 90 140 140]; %[105 155 155 105 155 155]; %[150 200 200 150 200 200] %[115 165 165 115 165 165]; %[150 200 200 150 200 200]; %[200 250 250 200 250 250] %fir second condition %300.*ones(1,6); % pixels

resize_factor = .155; %.175; %.225; %.3 %for second condition
resize_factorF = .14;
j_shift = [-10, 10];

h_pix = 0;
v_pix = 0;

black = [0 0 0]; % for fixation

% define response keys ...
KbName('UnifyKeyNames');
ResponseLeft = KbName('F'); % arrows
ResponseRight = KbName('J'); %
ResponseExit = KbName('P');
RestrictKeysForKbCheck([ResponseLeft ResponseRight ResponseExit]);
%----------------------------------

% Psi parameters ...

num_trials = 30; %0;
num_psi = 2;

miss_rate = 0.04;

xlevels = 201;
xstep = 0.1;
xbias = 0;

alevels = 201; %15;
astep = 0.1; %1.0;
abias = 0;

blevels = 9;
bratio = sqrt(2);
bbase = 1;

%----------------------------------

% intialise Psi's ...

%psi_counter = zeros(1,num_psi);

xmin = xbias - xstep*(xlevels-1)/2;
%xmax = xbias + xstep*(xlevels-1)/2;

x = zeros(1,xlevels);
x(1) = xmin;
for i = 2:xlevels
    x(i) = x(i-1) + xstep;
end

amin = abias - astep*(alevels-1)/2;
%amax = abias + astep*(alevels-1)/2;

a = zeros(1,alevels);
a(1) = amin;
for i = 2:alevels
    a(i) = a(i-1) + astep;
end

logbmin = log(bbase) - log(bratio)*(blevels-1)/2;
%logbmax = log(bbase) + log(bratio)*(blevels-1)/2;

b = zeros(1,blevels);
b(1) = exp(logbmin);
for i = 2:blevels
    b(i) = b(i-1)*bratio;
end

% set up values for the prior
ptl = ones(alevels,blevels,num_psi);
ptl = ptl./(alevels*blevels);   % uniform prior distribution

% set up look-up table of conditional probabilities
psuccess_lx = zeros(alevels,blevels,xlevels);

for i = 1:alevels
    for j = 1:blevels
        psuccess_lx(i,j,:) = logistic3([a(i) b(j) miss_rate],x);
    end
end

pts_x = zeros(size(x));    % prob of success after next trial as fn of stim, x

ptl_xsuccess = zeros(alevels,blevels,xlevels,num_psi);
ptl_xfailure = zeros(alevels,blevels,xlevels,num_psi);

%resp_arr = zeros(1,num_psi);    % set up response array

min_level = zeros(1,num_psi);

a_hat = zeros(1,num_psi);
a_err = zeros(1,num_psi);
b_hat = zeros(1,num_psi);
b_err = zeros(1,num_psi);

%----------------------------------------

% Set up graphics stuff ...

warning off;
rand('state',sum(100*clock));

% % initialization stuff from ContrastModulationDemo ...
AssertOpenGL;
% Open onscreen window on screen with maximum id:
%Screen('Preference','SkipSyncTests', 1); %remove later
screenid=max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
%i n t e n s e  antialiasing
%slows weaker computers
[win, winRect] = PsychImaging('OpenWindow', screenid, [], [], [], [], [], 16);

[screenWidth, screenHeight]=Screen('WindowSize', win);
SCREEN_X=screenWidth;
SCREEN_Y=screenHeight;

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% texture containing grey screen
blank(1) = Screen('MakeTexture',win,127.*ones(SCREEN_Y, SCREEN_X));
%blank(1,1) = Screen('MakeTexture',win,128.*ones(SCREEN_X, SCREEN_Y));

%Generate alpha layer blank texture:
%images = 128.*zeros(SCREEN_X, SCREEN_Y); %768);
%images(:,:,2) = 128;
%blank4ramp = Screen('MakeTexture',win,images);

% define raised cosine temporal window
pdc = 1;
cont = pdc.*ones(1,stim_durn_frames); % initialize to peak dot contrast (defined above)
win_length = stim_durn_frames/2;     % define window length
cont(1:win_length) = pdc.*0.5.*(1-cos(pi*(1:win_length)./win_length)); % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length); % ... and down
% keep an eye on timing ...
%vbl = zeros(1,stim_durn_frames);

%----------------------------------

HideCursor;

% Prompt for key press to start ...
%Screen('Blendfunction', win, GL_ONE, GL_ZERO);
Screen('FillRect', win, [128 128 128 0], [], 1);
%Screen('DrawTexture', win, blankTextures(1), [], [], [], 0);
Screen('DrawTexture',win, blank(1,1), [], [], [], 0);
DrawFormattedText(win, 'On each trial, indicate whether the middle face in the first or second interval\n is oriented more directly at you\n\nPress F for the 1st interval; press J for the 2nd interval.', 'center', 'center');
Screen('Flip',win);

% Wait for any button press to start ...
keypress = 0;
while ~keypress
    [keypress, keysecs,keyCode] = KbCheck;
end

% Show Initial Fixation
%Screen('DrawTexture',win, blank(1,1));
for fix=1:fix_frames
    Screen('DrawTexture', win, blank(1), [], [], [], 0);
    Screen('DrawDots', win, [(SCREEN_X-fix_size)/2 (SCREEN_Y-fix_size)/2], fix_size, [0, 0, 0],[0,0],1);
    Screen('Flip',win);
end

% Wait to flush out starting key press
Screen(win, 'WaitBlanking', 30);

FlushEvents('keyDown')
beep;

% %=
% something = imread(['C:\kieran\honoursweirdstuff\zawhy.bmp']);
% thing = Screen('MakeTexture',win,something);
% for i = 1:60
%     Screen('DrawTexture', win, thing, [], [], [], 0, []);
%     Screen('Flip', win);
% end

%----------------------------------------

num_stim = num_trials*num_psi;
ResponseArray = zeros(num_stim,6); %4);
psi_order = Shuffle(ceil([1:num_stim]./num_trials));
%usedtochoosefaces = randperm(6);
vector_pos = [cos(image_angle.*pi/180).*image_dist; sin(image_angle.*pi/180).*image_dist]';

for t = 1:num_stim
    
    % randomly select which psi to present ...
    psi = psi_order(t);
    
    t
    
    % 1. calculate conditional probability of response, r, given stimulus, x
    
    for k = 1:xlevels
        pts_x(k) = sum(sum(ptl(:,:,psi).*psuccess_lx(:,:,k)));
    end
    
    % 2. use Bayes rule to estimate posterior prob of each psycho fn ...
    
    for k = 1:xlevels
        ptl_xsuccess(:,:,k,psi) = (ptl(:,:,psi).*psuccess_lx(:,:,k))./pts_x(k);
        ptl_xfailure(:,:,k,psi) = (ptl(:,:,psi).*(1-psuccess_lx(:,:,k)))./(1-pts_x(k));
    end
    
    % 3. estimate entropy of pdf as a fn of stim level and response
    
    for k = 1:xlevels
        HS(k) = -sum(sum(ptl_xsuccess(:,:,k,psi).*log(ptl_xsuccess(:,:,k,psi))));
        HF(k) = -sum(sum(ptl_xfailure(:,:,k,psi).*log(ptl_xfailure(:,:,k,psi))));
    end
    
    % 4. estimate expected entropy for each stim level
    
    for k = 1:xlevels
        EH(k) = HS(k).*pts_x(k) + HF(k).*(1-pts_x(k));
    end
    
    % 5. find stim level with minimum entropy
    
    [min_val, min_level(psi)] = min(EH);
    
    % 6. RUN A TRIAL WITH STIMULUS min_level & ...
    
    % define stimulus image
    %psi
    test_head = sprintf('%2.1f',x(min_level(psi)));
    test_head = str2num(test_head);
    
    if psi == 1
        test_head = test_head + 5;
    elseif psi == 2
        test_head = test_head - 5;
    end
    
    if psi == 1
        test_head2 = 10 - test_head;
    elseif psi == 2
        test_head2 = -10 - test_head;
    end
    
    test_head3 = 0 - test_head;
    test_head4 = 0 - test_head2;
    
    test_head4 = sprintf('%2.1f',test_head4)    
    test_head3 = sprintf('%2.1f',test_head3)  
    test_head2 = sprintf('%2.1f',test_head2)
    test_head = sprintf('%2.1f',test_head)  %convoluted but works
    
    
    % initialise reaction timer
    FlushEvents('keyDown');
    keyIsDown = 0;
    responseStartTime = GetSecs;  % RT measured from onset of test image
    
    Screen('WaitBlanking', win);
    
    %merge
    usedtochoosering = randperm(6);
    usedtochoosefaces = randperm(6);
    whichway = 1; %randperm(2);
    flipIntervals = randperm(2);
    h_shift = randperm(2);
    v_shift = randperm(2);
    while isequal(h_shift,v_shift) == 1;
        h_shift = randperm(2);
        v_shift = randperm(2);
    end
    h_pix = j_shift(h_shift(1));
    v_pix = j_shift(v_shift(1));
    
    %pick surround faces
    for n = 1:num_images
        if usedtochoosering(1) == 1
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face1.png']);
        elseif usedtochoosering(1) == 2
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face2.png']);
        elseif usedtochoosering(1) == 3
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face3.png']);
        elseif usedtochoosering(1) == 4
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face4.png']);
        elseif usedtochoosering(1) == 5
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face5.png']);
        elseif usedtochoosering(1) == 6
            [I, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Face6.png']);
        end
        
        %pick mirrored surround faces
        if usedtochoosering(1) == 1
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace1.png']);
        elseif usedtochoosering(1) == 2
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace2.png']);
        elseif usedtochoosering(1) == 3
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace3.png']);
        elseif usedtochoosering(1) == 4
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace4.png']);
        elseif usedtochoosering(1) == 5
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace5.png']);
        elseif usedtochoosering(1) == 6
            [J, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\Ring\Gace6.png']);
        end
        
        I(:, :, 4) = alpha;
        J(:, :, 4) = alpha;
        stim(n) = Screen('MakeTexture',win,I);
        flop(n) = Screen('MakeTexture',win,J);
        
        if usedtochoosering(1) > 3
            dest_rect(n,:) = [SCREEN_X/2 SCREEN_Y/2 SCREEN_X/2 SCREEN_Y/2] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factor + [vector_pos(n,1) vector_pos(n,2) vector_pos(n,1) vector_pos(n,2)];
            dest_rect_shift(n,:) = [(SCREEN_X/2) (SCREEN_Y/2)+v_pix (SCREEN_X/2) (SCREEN_Y/2)+v_pix] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factor + [vector_pos(n,1) vector_pos(n,2) vector_pos(n,1) vector_pos(n,2)];
        else
            dest_rect(n,:) = [SCREEN_X/2 SCREEN_Y/2 SCREEN_X/2 SCREEN_Y/2] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factorF + [vector_pos(n,1) vector_pos(n,2) vector_pos(n,1) vector_pos(n,2)];
            dest_rect_shift(n,:) = [(SCREEN_X/2) (SCREEN_Y/2)+v_pix (SCREEN_X/2) (SCREEN_Y/2)+v_pix] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factorF + [vector_pos(n,1) vector_pos(n,2) vector_pos(n,1) vector_pos(n,2)];
        end
    end
    
    
    %ensure surround does not match center
    while usedtochoosering(1) == usedtochoosefaces(1);
        usedtochoosefaces = randperm(6);
    end
    
    
    if usedtochoosefaces(1) == 1
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude1\CEF1_h', test_head, '.bmp']);
    elseif usedtochoosefaces(1) == 2
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude2\CEF2_h', test_head, '.bmp']);
    elseif usedtochoosefaces(1) == 3
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude3\CEF3_h', test_head, '.bmp']);
    elseif usedtochoosefaces(1) == 4
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude4\CEM1_h', test_head, '.bmp']);
    elseif usedtochoosefaces(1) == 5
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude5\CEM2_h', test_head, '.bmp']);
    elseif usedtochoosefaces(1) == 6
        [appl, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude6\CEM3_h', test_head, '.bmp']);
    end
    
    if usedtochoosefaces(1) == 1
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude1\CEF1_h', test_head2, '.bmp']);
    elseif usedtochoosefaces(1) == 2
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude2\CEF2_h', test_head2, '.bmp']);
    elseif usedtochoosefaces(1) == 3
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude3\CEF3_h', test_head2, '.bmp']);
    elseif usedtochoosefaces(1) == 4
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude4\CEM1_h', test_head2, '.bmp']);
    elseif usedtochoosefaces(1) == 5
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude5\CEM2_h', test_head2, '.bmp']);
    elseif usedtochoosefaces(1) == 6
        [banan, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude6\CEM3_h', test_head2, '.bmp']);
    end
    
    if usedtochoosefaces(1) == 1
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude1\CEF1_h', test_head3, '.bmp']);
    elseif usedtochoosefaces(1) == 2
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude2\CEF2_h', test_head3, '.bmp']);
    elseif usedtochoosefaces(1) == 3
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude3\CEF3_h', test_head3, '.bmp']);
    elseif usedtochoosefaces(1) == 4
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude4\CEM1_h', test_head3, '.bmp']);
    elseif usedtochoosefaces(1) == 5
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude5\CEM2_h', test_head3, '.bmp']);
    elseif usedtochoosefaces(1) == 6
        [orang, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude6\CEM3_h', test_head3, '.bmp']);
    end
    
    if usedtochoosefaces(1) == 1
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude1\CEF1_h', test_head4, '.bmp']);
    elseif usedtochoosefaces(1) == 2
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude2\CEF2_h', test_head4, '.bmp']);
    elseif usedtochoosefaces(1) == 3
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude3\CEF3_h', test_head4, '.bmp']);
    elseif usedtochoosefaces(1) == 4
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude4\CEM1_h', test_head4, '.bmp']);
    elseif usedtochoosefaces(1) == 5
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude5\CEM2_h', test_head4, '.bmp']);
    elseif usedtochoosefaces(1) == 6
        [grap, ~, alpha] = imread(['C:\kieran\honoursweirdstuff\dude6\CEM3_h', test_head4, '.bmp']);
    end
    
    %orang = fliplr(appl);
    %grap = fliplr(banan);
    
    theface = Screen('MakeTexture',win,appl);
    theflip = Screen('MakeTexture',win,orang);
    theface2 = Screen('MakeTexture',win,banan);
    theflip2 = Screen('MakeTexture',win,grap);
    
    %female faces resized to match
    if usedtochoosefaces(1) > 3
        dest2rect = [SCREEN_X/2 SCREEN_Y/2 SCREEN_X/2 SCREEN_Y/2] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factor;
        dest2rect_shift = [(SCREEN_X/2) (SCREEN_Y/2)+v_pix (SCREEN_X/2) (SCREEN_Y/2)+v_pix] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factor;
    else
        dest2rect = [SCREEN_X/2 SCREEN_Y/2 SCREEN_X/2 SCREEN_Y/2] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factorF;
        dest2rect_shift = [(SCREEN_X/2) (SCREEN_Y/2)+v_pix (SCREEN_X/2) (SCREEN_Y/2)+v_pix] + [-size(I,2)/2+1 -size(I,1)/2+1 size(I,2)/2 size(I,1)/2].*resize_factorF;       
    end

    
    % Hold blank for a while to give processor a chance to catch up ...!
    %Screen(win, 'WaitBlanking', 60);
    
    vbl1 = Screen('Flip', win);
    %vbl_old = vbl1;
    % show image for set duration
    if whichway(1) == 1
        for jj=1:10
            Screen('DrawTexture', win, theface, [], dest2rect, [], 0, cont(jj));
            Screen('DrawTextures', win, stim, [], dest_rect', [], 0, cont(jj));
            Screen('DrawTexture', win, blank(1), [], [], [], 0);
            vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
            %vbl_old = vbl_new;
        end
    elseif whichway(1) == 2
        for jj=1:10
            Screen('DrawTexture', win, theflip, [], dest2rect, [], 0, cont(jj));
            Screen('DrawTextures', win, flop, [], dest_rect', [], 0, cont(jj));
            Screen('DrawTexture', win, blank(1), [], [], [], 0);
            vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
            %vbl_old = vbl_new;
        end
    end
    Screen('DrawTexture', win, blank(1), [], [], [], 0);
    vbl2 = Screen('Flip', win);
    
    %durn1 = vbl2-vbl1
    
    % ------------
    
    vbl1 = Screen('Flip', win);
    %vbl_old = vbl1;
    % show image for set duration
    if flipIntervals(1) == 1
        if whichway(1) == 1
            for jj=1:stim_durn_frames
                Screen('DrawTexture', win, theface, [], dest2rect, [], 0, cont(jj));
                Screen('DrawTextures', win, stim, [], dest_rect', [], 0, cont(jj));
                vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
                %vbl_old = vbl_new;
            end
            test_head
        elseif whichway(1) == 2
            for jj=1:stim_durn_frames
                Screen('DrawTexture', win, theflip, [], dest2rect, [], 0, cont(jj));
                Screen('DrawTextures', win, flop, [], dest_rect', [], 0, cont(jj));
                vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
                %vbl_old = vbl_new;
            end
            test_head3;
        end
    elseif flipIntervals(1) == 2
        if whichway(1) == 1
            for kk=1:stim_durn_frames
                Screen('DrawTexture', win, theface2, [], dest2rect_shift, [], 0, cont(kk));
                Screen('DrawTextures', win, flop, [], dest_rect_shift', [], 0, cont(kk));
                Screen('Flip',win);
            end
            test_head2
        elseif whichway(1) == 2
            for kk=1:stim_durn_frames
                Screen('DrawTexture', win, theface2, [], dest2rect_shift, [], 0, cont(kk));
                Screen('DrawTextures', win, stim, [], dest_rect_shift', [], 0, cont(kk));
                Screen('Flip',win);
            end
            %test_head2  
        end
    end
    
    Screen('DrawTexture', win, blank(1), [], [], [], 0);
    vbl2 = Screen('Flip', win);
    
    durn1 = vbl2-vbl1 % should be stim_durn_frames + 1 ...
    
    % ---------------------------
    
    usedtochoosefaces(1)
    usedtochoosering(1)
    
    for dd=1:delay_frames
        Screen('DrawTexture', win, blank(1), [], [], [], 0);
        Screen('Flip',win);
    end
    
    %Screen('WaitBlanking', win);
    vbl1 = Screen('Flip', win);
    
    % show image for set duration
    if flipIntervals(1) == 1    
        if whichway(1) == 1
            for kk=1:stim_durn_frames
                Screen('DrawTexture', win, theface2, [], dest2rect_shift, [], 0, cont(kk));
                Screen('DrawTextures', win, flop, [], dest_rect_shift', [], 0, cont(kk));
                Screen('Flip',win);
            end
            test_head2;
        elseif whichway(1) == 2
            for kk=1:stim_durn_frames
                Screen('DrawTexture', win, theface2, [], dest2rect_shift, [], 0, cont(kk));
                Screen('DrawTextures', win, stim, [], dest_rect_shift', [], 0, cont(kk));
                Screen('Flip',win);
            end
            %test_head2
        end
    elseif flipIntervals(1) == 2
        if whichway(1) == 1
            for jj=1:stim_durn_frames
                Screen('DrawTexture', win, theface, [], dest2rect, [], 0, cont(jj));
                Screen('DrawTextures', win, stim, [], dest_rect', [], 0, cont(jj));
                vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
                %vbl_old = vbl_new;
            end
            test_head;
        elseif whichway(1) == 2
            for jj=1:stim_durn_frames
                Screen('DrpawTexture', win, theflip, [], dest2rect, [], 0, cont(jj));
                Screen('DrawTextures', win, flop, [], dest_rect', [], 0, cont(jj));
                vbl_new = Screen('Flip',win); %,vbl_old+(1/120));
                %vbl_old = vbl_new;
            end
            test_head3;
        end
    end

    
    vbl2 = Screen('Flip', win);
    
    durn2 = vbl2-vbl1
    
    Screen('DrawTexture', win, blank(1), [], [], [], 0);
    
    Screen('Flip',win);
    
    ID %wtf i changed nothing
    
    % Collect response
    while ~keyIsDown
        [keyIsDown, keysecs,keyCode] = KbCheck;
    end
    beep;
    
    if keyIsDown
        if keyCode(ResponseLeft)
            Response = 0; %1;
            ResponseTime = keysecs - responseStartTime; % time in msecs
        elseif keyCode(ResponseRight)
            Response = 1; %0;
            ResponseTime = keysecs - responseStartTime; % time in msecs
        elseif keyCode(ResponseExit)
            % break out of program on any other response
            Screen('CloseAll');
            warning on;
            clear hidecursor;
            crash_out_here;
        end
    end
    
    % Hold blank for a while ...
    Screen('DrawTexture', win, blank(1), [], [], [], 0);
    Screen('Flip',win);
            
    Screen(win, 'WaitBlanking', 120);
    
    % ... then display fixation
    Screen('DrawDots', win, [(SCREEN_X-fix_size)/2 (SCREEN_Y-fix_size)/2], fix_size, black,[0,0],1);
    Screen('Flip',win);
    
    Screen(win, 'WaitBlanking', 60);
    
    % 7. keep the posterior pdf that corresponds to the completed trial
    
    %     if (Response)
    %         ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
    %     else
    %         ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
    %     end
    
%     if (Response) % CC - recoded responses for psi==2
%         if (psi==1)
%             ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
%         elseif (psi==2)
%             ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
%         end
%     else
%         if (psi==1)
%             ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
%         elseif (psi==2)
%             ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
%         end
%     end
    
    if (Response) % CC - recoded responses
        if (psi==1)
            if (flipIntervals(1)==1)
                ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
            elseif (flipIntervals(1)==2)
                ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
            end
        elseif (psi==2)
            if (flipIntervals(1)==1)
                ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
            elseif (flipIntervals(1)==2)
                ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
            end
        end    
    else
        if (psi==1)
            if (flipIntervals(1)==1)
                ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
            elseif (flipIntervals(1)==2)
                ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
            end
        elseif (psi==2)
            if (flipIntervals(1)==1)
                ptl(:,:,psi) = ptl_xsuccess(:,:,min_level(psi),psi);
            elseif (flipIntervals(1)==2)
                ptl(:,:,psi) = ptl_xfailure(:,:,min_level(psi),psi);
            end
        end
    end
    
    % write stimulus & response details to ResponseArray ...
    running_a_hat = sum(a*ptl(:,:,psi)); % this is just for saving ...
    ResponseArray(t,:) = [psi, x(min_level(psi)), Response, ResponseTime, running_a_hat, (flipIntervals(1)-1)];
    
    % tidy up after each trial - CC
    %Screen('Close',[stim, flop, theface, theflip, theface2, theflip2]);
    
end % end of stuff done every trial

% 8. Find new estimate of psychometric function
for psi=1:num_psi
    a_hat(psi) = sum(a*ptl(:,:,psi));
    b_hat(psi) = sum(b*ptl(:,:,psi)');
    
    a_mean = ones(size(a)).*a_hat(psi);
    b_mean = ones(size(b)).*b_hat(psi);
    
    a_err(psi) = sqrt(sum(((a-a_mean).^2)*ptl(:,:,psi)));
    b_err(psi) = sqrt(sum(((b-b_mean).^2)*ptl(:,:,psi)'));     % standard error
end

a_hat
a_err

TI = 0.5*(a_hat(2)-a_hat(1))
TI_err = sqrt(((a_err(1)^2)+(a_err(2)^2))/2)

% write ResponseArray to file ...
psy_fun = [];
for z = 1:num_psi
    psy_fun = [psy_fun; a_hat(z) a_err(z) b_hat(z) b_err(z) 0 0];
end
%dlmwrite(filename,[ResponseArray;psy_fun],'\t');

plot_psi_results(ResponseArray,psy_fun);

dlmwrite(filename,[ResponseArray;psy_fun],'\t');

for jj=1:stim_durn_frames
    Screen('TextSize',win,20);
    Screen('TextStyle',win,0);
    DrawFormattedText(win,'This portion of the study is now complete\n\nPlease call a researcher to proceed.',...
        'center',420,black,[],[],[],1.65);
    Screen('Flip',win,[],[],[],1); %display on screen
end

keypress = 0;
while ~keypress
    [keypress, ~,~] = KbCheck;
end

% return to command window ...
ShowCursor;
Screen('CloseAll');
warning on;
clear HideCursor;