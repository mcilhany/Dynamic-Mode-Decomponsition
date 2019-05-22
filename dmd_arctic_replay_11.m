%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Read Arctic Sea Ice lifetimes for last 25 years from NOAA
%%%%    "Report Card" from Youtube
%%%%    Kevin Mcilhany  09-JAN-2019
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clf
colormap jet(100)

FileNamee = 'arctic_25_years.mp4'
% Color is converted using either BLUE only or all colors to grayscale
% Land is removed, then put back later in order to reduce matrix size
for blahinitialize = 1:1
    framespersecond     = 30;           % Video frame rate  
    deltat              = 1/framespersecond;  % Time in seconds of each frame
    dt                  = deltat;
    starttimee0         = 12.0667;      % Time in seconds when the video begins showing valid data
    endtimee            = 56.66667;     %  Time in seconds when the video ends showing valid data
    numyears            = 25.75;        %  Number of years showing valid data (JAN1990 - SEP2015)
    numyears            = 1;
    xmin                = 340;
    xmax                = 1140;
    ymin                = 60;
    ymax                = 600;
    pixsize             = 740;
    
    dmdperiodicityyears = 1;
    framesperyear       = 52;           % Set the DMD time window to ONE year (52 frames)
    yearlongtime        = framesperyear*deltat;
    dmdframesize        = ceil(dmdperiodicityyears*framesperyear);
    lastframe           = 21000;
    lastframe           = ceil(min(lastframe,dmdframesize*numyears));
    movietimesize       = (lastframe-1)*deltat;         % Set the DMD time window to ONE year (52 frames)
    timee               = linspace(0,movietimesize,lastframe);   % Time used for the prediction video (2 years) 
    yearstart           = 1;            %  Which year to start the analysis on (1990 = #1)
    yearend             = 1;            %  Which year to end the analysis on  
%    yearend             = ceil(numyears);  
    if yearend<yearstart; yearend = yearstart;   end 
    starttimee    = starttimee0 + yearlongtime*(yearstart-1);   % Starting time in video of the first analysis year

    maxmodenumber       = 51             % Maximum number of DMD modes to apply 
%    maxmodenumber      = 10000;        % Apply as many as possible
    if maxmodenumber>(dmdframesize-1); maxmodenumber = dmdframesize-1; end
    colornormalizee     = 1;
    badmodescut         = -0.60;        % Cut out DMD modes with an abs(eigvalss)<threshold = badmodescut
    framemaskid         = 1;            % Mask with NOTHING removed

    watchframes = 0;
    makeimages = 0;
    saveasonee = 0;

    for blahframemask=1:1
        if framemaskid==1
            framee1crit     = 0;              % The framemaskid is an attempt to remove "dead" pixels from the video 
            framee2crit     = 0;              %      so that they are not treated as part of the data used in DMD
            framee3crit     = 0;              %      Examples of this are the Land and the Text in the video
            framee1thresh   = 0;              % frameecrit = is the threshold to cross for when to count a color value change 
            framee2thresh   = 0;              % frameethresh = overall years, per pixel, how many pixels failed to change color
            framee3thresh   = 0;              %                by more than this amount determines if a pixel is "dead" 
            bluegray        = 0;              % bluegray = 1 uses only BLUE,  = 0 uses all colors to convert to grayscale
        elseif framemaskid==2
            framee1crit     = 1;              % The framemaskid is an attempt to remove "dead" pixels from the video 
            framee2crit     = 1;              %      so that they are not treated as part of the data used in DMD
            framee3crit     = 1;              %      Examples of this are the Land and the Text in the video
            framee1thresh   = 2;              % frameecrit = is the threshold to cross for when to count a color value change 
            framee2thresh   = 2;              % frameethresh = overall years, per pixel, how many pixels failed to change color
            framee3thresh   = 2;
            bluegray        = 0;              % bluegray = 1 uses only BLUE,  = 0 uses all colors to convert to grayscale
        elseif framemaskid==3
            framee1crit     = 2;              % The framemaskid is an attempt to remove "dead" pixels from the video 
            framee2crit     = 2;              %      so that they are not treated as part of the data used in DMD
            framee3crit     = 2;              %      Examples of this are the Land and the Text in the video
            framee1thresh   = 3;              % frameecrit = is the threshold to cross for when to count a color value change 
            framee2thresh   = 3;              % frameethresh = overall years, per pixel, how many pixels failed to change color
            framee3thresh   = 3;
            bluegray        = 0;              % bluegray = 1 uses only BLUE,  = 0 uses all colors to convert to grayscale
        elseif framemaskid==4
            framee1crit     = 1;              % The framemaskid is an attempt to remove "dead" pixels from the video 
            framee2crit     = 1;              %      so that they are not treated as part of the data used in DMD
            framee3crit     = 1;              %      Examples of this are the Land and the Text in the video
            framee1thresh   = 42;              % frameecrit = is the threshold to cross for when to count a color value change 
            framee2thresh   = 42;              % frameethresh = overall years, per pixel, how many pixels failed to change color
            framee3thresh   = 22;
            bluegray        = 0;              % bluegray = 1 uses only BLUE,  = 0 uses all colors to convert to grayscale
        elseif framemaskid==5
            framee1crit     = 1;              % The framemaskid is an attempt to remove "dead" pixels from the video 
            framee2crit     = 1;              %      so that they are not treated as part of the data used in DMD
            framee3crit     = 1;              %      Examples of this are the Land and the Text in the video
            framee1thresh   = 42;              % frameecrit = is the threshold to cross for when to count a color value change 
            framee2thresh   = 42;              % frameethresh = overall years, per pixel, how many pixels failed to change color
            framee3thresh   = 22;
            bluegray        = 1;              % bluegray = 1 uses only BLUE,  = 0 uses all colors to convert to grayscale
        end
        framemask = sprintf('_DATA/framemask_%02d_%02d_%02d_%02d_%03d_%03d_%03d_%01d.mat',framemaskid,framee1crit,framee2crit,framee3crit,framee1thresh,framee2thresh,framee3thresh,bluegray)
        load(framemask)
    end
    
    pixaspect = 10.7/8;
    handf = figure('Units','pixels','PaperPosition',[1 1 pixsize pixsize*pixaspect],'PaperSize',[pixsize,pixsize*pixaspect],'Position',[1 1 pixsize pixsize*pixaspect])
    set(handf,'Visible','on','PaperPositionMode', 'manual','Units','inches','PaperPosition',[0,0,6,12],'PaperSize',[6,12])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)]; 
end


yearfirst = 1;                          % Flag for the first pass through the analysis
for yeari = yearstart:yearend           % Analysis loop over all years

    yearstop = min(yeari,numyears);     % Determine when to stop looking at the video, for the last year,                          
    endtimee      = starttimee + yearlongtime*yearstop;   %  it may stop short, if so, adjust the number of DMD modes available as well    

    for blahvideoprocess = 1:1    
        arctic = VideoReader(FileNamee);     % Read in the video
        arctic.CurrentTime = starttimee;
        currAxes = axes;
        framei = 1;
        while arctic.CurrentTime<=endtimee    % Play the video for the current year
            vidFrame = readFrame(arctic);
            if bluegray==1
                framee = vidFrame(ymin:ymax,xmin:xmax,3);
            else
                framee = rgb2gray(vidFrame(ymin:ymax,xmin:xmax,:));
            end
            framee(frameenochangeind) = 0;    %     %   Apply MASK    

            if bluegray==1              % shift the gray scale image which is too saturated in white
                image((framee-55), 'Parent', currAxes)     % Base shift for BLUE only
            else
                image((framee-45), 'Parent', currAxes)     % Base shift for grayscale from all colors
            end

            currAxes.Visible = 'off';
            drawnow
            framei = framei+1;
        end
        totalframes = framei-1;           %  Count of the video frames for the current year  (caution it may have stopped short)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%   Process the Video
        %%%%   Prepartation for  Dynamic Mode Decomposition
        %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [xgrid,ygrid] = meshgrid(1:size(framee,2),1:size(framee,1));
        ygrid = flipud(ygrid);
        xgrid(frameenochangeind) = [];    %   Apply MASK to x,y grid
        ygrid(frameenochangeind) = [];    %   Apply MASK to x,y grid
        clf
        arctic = VideoReader(FileNamee);
        arctic.CurrentTime = starttimee;
        vidFrame = readFrame(arctic);
        framee = rgb2gray(vidFrame(ymin:ymax,xmin:xmax,:));
        framee(frameenochangeind) = [];   %   Apply MASK to data
        scatter(xgrid,ygrid,1,framee)     %   Show the data with the mask applied

        currAxes = axes;
        framei = 1;
        bigx     = zeros(size(framee,2),totalframes);     % Setup a matrix to contain ALL of the data for ALL frames
        bigxnext = zeros(size(framee,2),totalframes);     % Setup a matrix to contain ALL of the data for ALL frames
        while arctic.CurrentTime<=endtimee
            vidFrame = readFrame(arctic);
            if bluegray==1
                framee = vidFrame(ymin:ymax,xmin:xmax,3);
            else
                framee = rgb2gray(vidFrame(ymin:ymax,xmin:xmax,:));
            end 
            framee(frameenochangeind) = [];          %   Apply MASK   
            bigx(:,framei) = framee(:);              %   Store data in matrix bigx
            framei = framei+1;
        end
        totalframes = framei-1;
        if maxmodenumber>totalframes;   maxmodenumber=totalframes;  end   % Make sure you only solve for the number of DMD that 

        framei = 1;                     %   Also, store the data into another matrix for the comparison video at the end, bigxnext
        if yeari<numyears               %   Only do this for years before the last year
            endtimenext      = starttimee0 + yearlongtime*(yeari+1);
            while arctic.CurrentTime<endtimenext
                vidFrame = readFrame(arctic);
                if bluegray==1
                    framee = vidFrame(ymin:ymax,xmin:xmax,3);
                else
                    framee = rgb2gray(vidFrame(ymin:ymax,xmin:xmax,:));
                end
                framee(frameenochangeind) = [];    %   Apply MASK  
                bigxnext(:,framei) = framee(:);
                framei = framei+1;
            end
            nobigxnext = 0;
        else 
            nobigxnext = 1;
        end
        bigxx = ((bigx-56)>0).*ceil((bigx-56)/20);
        bigxxnext = ((bigxnext-56)>0).*ceil((bigxnext-56)/20);
        bigx = bigxx;
        bigxnext = bigxxnext;
    end

    for blahDMD = 1:1    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%   Dynamic Mode Decomposition
    %%%%    Here it is (finally) - just three  lines 
    %%%%    once you have the BigX matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colormap jet(100)
    bigxx = ((bigx-56)>0).*ceil((bigx-56)/20);
    X1 = bigx(:,1:end-1);           %   The following is the exact code used from Brunton's book for DMD execution
    X2 = bigx(:,2:end);
    r = maxmodenumber; % rank truncation
    [U, S, V] = svd(X1, 'econ');
    Ur = U(:, 1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:, 1:r);
    %% Build Atilde and DMD Modes
    Atilde = Ur'*X2*Vr/Sr;
    [W, D] = eig(Atilde);
    Phi = X2*Vr/Sr*W;  % DMD Modes

    %% DMD Spectra
    lambda = diag(D);
    omega = log(lambda)/dt;
    eigmat = W;
    eigvals = D;
    eigvalss = diag(eigvals);
    uu = U;
    end

    for blahplots = 1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DONE!
    %%%  Save LOTS of pictures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        clf
        plot(real(eigvalss),imag(eigvalss),'*','MarkerSize',8)
        rectangle('Position',[-1 -1 2 2],'Curvature',[1,1],'linewidth',2);
        daspect([1,1,1])
        axis equal 
        for eigi = 1:maxmodenumber
            text(real(eigvalss(eigi)),imag(eigvalss(eigi)),num2str(eigi))
        end
        title(sprintf('DMD Arctic 25 Eigenvalues Year:%02d  %03d Modes',yeari,maxmodenumber))
        xlabel('Real')
        ylabel('Imaginery')
        axis tight equal
        %    saveas(gcf,sprintf('_FIGS/dmd_arctic_year_%02d_eigenvalues.jpg',yeari));
        figureimage=getframe(gcf);
        imwrite(figureimage.cdata,sprintf('_FIGS/dmd_arctic_year_%02d_eigenvalues_%03dmodes_%02dmask.tiff',yeari,maxmodenumber,framemaskid),'compression','jpeg','rowsperstrip',8)

        clf
        plot(omega,'*','MarkerSize',8)
        for eigi = 1:maxmodenumber
            text(real(omega(eigi)),imag(omega(eigi)),num2str(eigi))
        end
        title(sprintf('DMD Arctic 25    Ritz values Year:%02d    %03d Modes',yeari,maxmodenumber))
        xlabel('Real')
        ylabel('Imaginery')
        axis auto 
        %    saveas(gcf,sprintf('_FIGS/dmd_arctic_year_%02d_ritzvalues.jpg',yeari));
        figureimage=getframe(gcf);
        imwrite(figureimage.cdata,sprintf('_FIGS/dmd_arctic_year_%02d_ritzvalues_%03dmodes_%02dmask.tiff',yeari,maxmodenumber,framemaskid),'compression','jpeg','rowsperstrip',8)

        for blahimages = 1:1
            if makeimages==1
                for dmdi = 1:maxmodenumber
                    dmdmoder(:,dmdi) = real(Phi(:,dmdi));
                    dmdmodei(:,dmdi) = imag(Phi(:,dmdi));
                end
                cminr = min(min(dmdmoder));
                cmaxr = max(max(dmdmoder));
                cmini = min(min(dmdmodei));
                cmaxi = max(max(dmdmodei));
                clf
                subplot(2,1,1)
                dmdi = 1;
                scatter(xgrid,ygrid,1,dmdmoder(:,dmdi))
                title(sprintf('DMD Year:%02d   Real Component Mode #%04d',yeari,dmdi));
                axis tight equal
                xlabel('Color Normalize = n')
                subplot(2,1,2)
                scatter(xgrid,ygrid,1,dmdmoder(:,dmdi))
                caxis([cminr cmaxr]);
                title(sprintf('DMD Year:%02d   Real Component Mode #%04d',yeari,dmdi));
                axis tight equal
                xlabel('Color Normalize = y')
                figureimage=getframe(gcf);
                imwrite(figureimage.cdata,sprintf('_FIGS/dmd_arctic_year_%02d_real_colornormalize.tiff',yeari),'compression','jpeg','rowsperstrip',8)

                for dmdi = 1:maxmodenumber
                    clf
                    scatter(xgrid,ygrid,1,dmdmoder(:,dmdi))
            %        pcolor(dmdmoder(:,dmdi))
                    if colornormalizee==1;        caxis([cminr cmaxr]);   end
                    title(sprintf('DMD Arctic 25 Year Sea Ice Lifetime Year:%02d   Real Component Mode #%04d  %03d Modes',yeari,dmdi,maxmodenumber));
                    axis tight equal 
            %        saveas(gcf,sprintf('_FIGS/dmd_arctic_year_%02d_real_%04d.jpg',yeari,dmdi));
                    set(gcf,'Visible','on','PaperPositionMode', 'manual','Units','inches','PaperPosition',[0,0,8.2,8.0],'PaperSize',[8.2,8.0])
                    set(gca,'Position',[0.03 0.03 0.94 0.93])
                    figureimage=getframe(gcf);
                    imwrite(figureimage.cdata,sprintf('_FIGS/dmd_arctic_year_%02d_real_%03d_modes.tiff',yeari,maxmodenumber),'compression','jpeg','rowsperstrip',8)
                    clf
                    scatter(xgrid,ygrid,1,dmdmodei(:,dmdi))
            %        pcolor(dmdmodei(:,dmdi))
                    if colornormalizee==1;        caxis([cmini cmaxi]);   end
                    title(sprintf('DMD Arctic 25 Year Sea Ice Lifetime Year:%02d   Imaginery Component Mode #%04d %03d Modes',yeari,dmdi,maxmodenumber));
                    axis tight equal
            %        saveas(gcf,sprintf('_FIGS/dmd_arctic_year_%02d_imaginery_%04d.jpg',yeari,dmdi));
                    set(gcf,'Visible','on','PaperPositionMode', 'manual','Units','inches','PaperPosition',[0,0,8.2,8.0],'PaperSize',[8.2,8.0])
                    set(gca,'Position',[0.03 0.03 0.94 0.93])
                    figureimage=getframe(gcf);
                    imwrite(figureimage.cdata,sprintf('_FIGS/dmd_arctic_year_%02d_imaginery_%03d_modes.tiff',yeari,maxmodenumber),'compression','jpeg','rowsperstrip',8)
                end
            end
        end
    end

    dmdmodes = Phi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Save the data to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(sprintf('_DATA/dmd_arctic_year_%02d_modes_%02dmask.mat',yeari,framemaskid),'bigx','eigvalss','dmdmodes','omega','timee','bigxnext','totalframes','nobigxnext','dmdframesize','xgrid','ygrid','yeari','starttimee','endtimee','deltat','framemaskid');
    starttimee    = endtimee;
    Phi = dmdmodes;

    for blahpredict = 1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Rebuild the data from the DMD modes.  Begin by calculating
    %%%   the amplitues from the first frame of the video projected onto
    %%%   each of the DMD modes, yielding,  b0
    %%%   Then, from the Arnoldi technique, evolve the system for each frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    badmodesi = find(abs(eigvalss)<badmodescut);
    Phi(:,badmodesi) = [];
    omega(badmodesi) = [];
    maxmodenumber = size(Phi,2)
    b0 = Phi\bigx(:, 1);
    X_pred = zeros(numel(omega), length(timee));
    for tt = 1:length(timee)     %  Prediction window extends across the full year (hindcasting) and adds another full year (forecasting)
        X_pred(:, tt) = b0 .* exp(omega .* timee(tt));
    end;
    X_pred = Phi * X_pred;
    dmdpredictionreal = real(X_pred);
    end

    for blahvideo = 1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%   Play the data back from the BigX matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = VideoWriter(sprintf('_FIGS/dmd_arctic_year_%02d_playback_%03dmodes_%6.4fdeltat_%02dmask.avi',yeari,maxmodenumber,deltat,framemaskid));
    open(v);
    clf
    for framei = 1:lastframe
        yearj = yeari+floor(framei/dmdframesize);
        subplot(2,1,1)
        scatter(xgrid,ygrid,1,dmdpredictionreal(:,framei))
        title(sprintf('DMD Reconstructed Simulation using %03d Modes:  Year: %02d  frame (52/year):%03d',maxmodenumber,yearj,framei))
        axis tight equal
        set(gca,'nextplot','replacechildren'); 
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        leftt = outerpos(1) + ti(1);
        bottomm = (outerpos(2) + ti(2))*0.92;
        ax_width = (outerpos(3) - ti(1) - ti(3))*0.95;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [leftt bottomm ax_width ax_height];
        set(gca,'Position',[0.09 0.52 0.87 0.450])

        subplot(2,1,2)
        scatter(xgrid,ygrid,1,bigx(:,framei))
        title(sprintf('Actual Video using %03d Modes:  Year: %02d  frame (52/year):%03d',maxmodenumber,yearj,framei))
        axis tight equal
        set(gca,'nextplot','replacechildren'); 
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        leftt = outerpos(1) + ti(1);
        bottomm = (outerpos(2) + ti(2))*0.92;
        ax_width = (outerpos(3) - ti(1) - ti(3))*0.95;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [leftt bottomm ax_width ax_height];
        set(gca,'Position',[0.09 0.02 0.87 0.450])
       frame = getframe(gcf);
       writeVideo(v,frame);
       %pause
    end
    close(v);
    end

end

















