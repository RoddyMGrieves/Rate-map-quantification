
%% mapping settings
    mapset                                              = struct;
    mapset.padding                                      = 0; % in mm
    mapset.mindwell                                     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist                                      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.maxdist                                      = 640; % (mm, default 50) used by kadaptive, adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.srate                                        = 50; % (default 50) sampling rate of data in Hz, used to calculate time
    mapset.steps                                        = 32; % the number of convolution size steps to use for kadaptive    
    mapset.kern                                         = 'biweight'; % kernel
    mapset.smethod                                      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing
    mapset.methods                                      = {'histogram','ash','ksde','kyadaptive','kadaptive','fyhn'}; % histogram ash ksde kyadaptive kadaptive fyhn
    mapset.methods                                      = {'histogram'};    
    mapset.twindow                                      = 0.25; % time window (s) over which to estimate instantaneous firing for temporal methods

    close all;

    % get position data (random walks)
    if 0
        nwalks = config.nwalks;
        all_pos = cell(nwalks,1);
        for ww = 1:nwalks
            wlength = 64;
            wdir = [scratch_space '\' config.environs{1} '_' num2str(wlength) '_walk' num2str(ww) '.mat'];
            load(wdir,'pox','poy','pot','epoly','emat','wlength'); % load walk data       
            
            pox = fillmissing(pox(:),'linear');
            poy = fillmissing(poy(:),'linear');
            
            duration = numel(pox)*(1/50);
            all_pos(ww,1) = {[pox(:) poy(:) pot(:)]};
        end
    end
    
    if 0
        sdir = [scratch_space '\' ename '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata');     
    end


% figure
% for ii = 1:25
%     subplot(5,5,ii)
%         imagesc(sdata.pmap_1mm{ii})
%         daspect([1 1 1])
%         colormap('turbo')
%         title(sprintf('%d',ii))
% end










%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % create a figure
    figure('Position',[100 50 900 850]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','none'); % makes the background colour white

    qimg = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\PLoS Comp Bio submission 2\striking\Picture1.png';
    fig_dir = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\outputs'; % desk computer + bonsai computer
    

xmin = 50;
xmax = 850;
psiz = 62;
xbuff = 5;
xvec = xmin : (psiz+xbuff) : xmax;
ymin = 50;
ymax = 750;
ybuff = 5;
yvec = ymax : -(psiz+ybuff) : ymin;
[xx,yy] = meshgrid(xvec,yvec);
xx = xx';
yy = yy';

    ncol = length(xvec);
    nrow = length(yvec);
    nplots = ncol*nrow;
    mainp = [55];
    rmap_colormap = turbo;
    
    maps_for_plotting = cell(length(nplots),1);
    fname = [fig_dir '\striking_maps.mat'];
    override_maps = 0;
    if exist(fname,'file') && ~override_maps
        mem = 1;
        disp(sprintf('\t\t...loading old maps'))
        load(fname,'maps_for_plotting')
    else
        disp(sprintf('\t\t...generating maps'))            
        mem = 0;
    end

    mnames = {
            'histogram',[1],[70],[5];
            'histogram',[2],[60],[8];
            'histogram',[5],[30],[9];
            'histogram',[10],[40],[17];
            'histogram',[20],[40],[18];
            'histogram',[30],[30],[19];
            'histogram',[40],[40],[23];
            'histogram',[50],[40],[29];
            'histogram',[60],[0],[23];
            'histogram',[100],[160],[23];
            'ash',[60],[20],[38];
            'ash',[65],[18],[35];
            'ash',[70],[16],[34];
            'ash',[75],[14],[33];
            'ash',[80],[12],[25];
            'ash',[85],[10],[22];
            'ash',[90],[8],[21];
            'ash',[80],[8],[25];   
            'ash',[160],[4],[21];
            'ash',[120],[8],[25];             
            'kadaptive',[2],[16000],[12];
            'kadaptive',[5],[4000],[10];
            'kadaptive',[10],[3000],[25];
            'kadaptive',[20],[1500],[22];
            'kadaptive',[15],[2000],[21];
            'kadaptive',[3],[8000],[11];
            'kadaptive',[40],[4000],[12];
            'kadaptive',[2],[4000],[12];   
            'kadaptive',[60],[2000],[12];
            'kadaptive',[2],[256],[12];             
            'kyadaptive',[1],[1],[21];
            'kyadaptive',[2],[2],[21];
            'kyadaptive',[3],[3],[21];
            'kyadaptive',[5],[4],[21];
            'kyadaptive',[10],[1],[21];
            'kyadaptive',[15],[2],[21];
            'kyadaptive',[20],[3],[21];
            'kyadaptive',[5],[4],[21];
            'kyadaptive',[10],[1.5],[21];
            'kyadaptive',[75],[2],[21];
            'kyadaptive',[65],[2],[21];
            'ksde',[7],[30],[21];
            'ksde',[9],[25],[21];
            'ksde',[11],[30],[21];
            'ksde',[13],[35],[21];
            'ksde',[15],[15],[21];
            'ksde',[17],[20],[21];
            'ksde',[50],[25],[21];
            'fyhn',[7],[35],[21];
            'fyhn',[9],[30],[21];
            'fyhn',[11],[35],[21];
            'fyhn',[13],[35],[21];
            'fyhn',[15],[20],[21];
            'fyhn',[17],[25],[21];
            'fyhn',[20],[30],[21];
            'fyhn',[40],[90],[21];
            'fyhn',[65],[30],[21];            
            };
    
mnames2 = {'histogram',[100],[300]};

q_pixels = [30 31 41 44 56 67 78 102];
q_pixels = [];

    rr = randperm(size(sdata,1),nplots);
    for ii = 1:nplots
        disp(sprintf('\t\t...%d',ii))            
        
        % ax = subaxis(nrow,ncol,ii,'SpacingHoriz',0.005,'SpacingVert',0,'Holdaxis'); % 'SpacingHoriz',0.1,'SpacingVert',0.1,
        ax = axes('Units','pixels','Position',[xx(ii) yy(ii) psiz psiz],'color','none');
            axis square
            % excell = mnames{ii,4};
            % excell = 19;
            excell = rr(ii);

            % cut the position and spike data to the trial length
            % we want
            pindx = sdata.pos_index(excell);
            pos_now = all_pos{pindx,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3); 

            spk = sdata.spk{excell};
            spx = spk(:,1);
            spy = spk(:,2);
            spt = spk(:,3);

            tcut = 16*60;
            ppox = pox(pot<tcut);
            ppoy = poy(pot<tcut);
            ppot = pot(pot<tcut);            
            pspx = spx(spt<tcut);
            pspy = spy(spt<tcut);
            pspt = spt(spt<tcut);            
            pos = [ppox ppoy]; % positions in mm
            spk = [pspx pspy]; % spikes in mm
% continue
            mindx = randi(size(mnames,1),[1 1]);
            if ismember(ii,q_pixels)
                mindx = randi(size(mnames2,1),[1 1]);
            end
            % mindx = 7;
            if mem & ~ismember(ii,q_pixels)
                ratemap = maps_for_plotting{ii,1};
            else                    
                rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                if ismember(ii,q_pixels)
                    rmset.method = mnames2{mindx,1};
                    rmset.binsize = mnames2{mindx,2};
                    rmset.ssigma = mnames2{mindx,3};
                else
                    rmset.method = mnames{mindx,1};
                    rmset.binsize = mnames{mindx,2};
                    rmset.ssigma = mnames{mindx,3};
                end

% rmset.method = mnames{1,1};
% rmset.binsize = mnames{1,2};
% rmset.ssigma = mnames{1,3};

                rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
                rmset.pot = ppot;
                rmset.spt = pspt;
                if strcmp(rmset.method,'fyhn')
                    rmset.maxdist = 320;
                    rmset.mindist = 50;                                
                end
                [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                maps_for_plotting{ii,1} = ratemap;
            end

            im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
            daspect([1 1 1])
            axis xy off tight square
            colormap(gca,rmap_colormap)     
            % text(0.5,0.5,sprintf('%d',ii),'HorizontalAl','left','Units','normalized','FontSize',8,'color','w')
            ax.CLim = [0 max(ratemap(:),[],'omitnan')];

    end

    ax = axes('Units','pixels','Position',[xvec(6) yvec(6) psiz psiz]);
        qi = imread(qimg);
        J = ~logical(im2gray(qi));
        d = bwdist(~J,'euclidean');
        d(d>0) = d(d>0)+normrnd(0,180,size(d(d>0)));
        d(d==0) = d(d==0)+normrnd(0,1,size(d(d==0)));
        
        % ds = normpdf(d,0,10); 
        ds = imgaussfilt(d,20);
        ds = imresize(ds,0.07,'nearest');
        % ds = padarray(ds,[80 0],0,'post'); % move up
        % ds = padarray(ds,[0 100],0,'pre'); % right
        % ds = padarray(ds,[30 0],0,'pre'); % move down
        % ds = padarray(ds,[0 10],0,'post'); % left

        ds = rot90(ds,-1);
        im = imagesc(ds,'alphadata',~isnan(ds)); hold on;
        daspect([1 1 1])
        axis ij off tight square
        colormap(gca,rmap_colormap)  
        ax.CLim = [0 max(ds(:),[],'omitnan')];

        exportgraphics(gcf,['plos_strike.eps'],'BackgroundColor','none','Colorspace','rgb','Resolution',350);  

        save(fname,'maps_for_plotting')









