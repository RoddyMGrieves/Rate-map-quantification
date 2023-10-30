



c_selected = {
    'r926_09042018_t12_c1','926_09042018.mat';
    'r926_09042018_t12_c2','926_09042018.mat';
    'r926_09042018_t12_c5','926_09042018.mat';
    'r852_11082017_t7_c3','852_11082017.mat';
    'r852_16032017_t8_c8','852_16032017.mat';
    'r852_16032017_t8_c9','852_16032017.mat';    
    'r896_11012018_t13_c4','896_11012018.mat';
    'r896_16012018_t15_c2','896_16012018.mat';
    'r896_16012018_t14_c7','896_16012018.mat';
    'r926_13042018_t12_c4','926_13042018.mat';
    'r926_13042018_t14_c5','926_13042018.mat';
};

dat = table;
for ii = 1:size(c_selected,1)
    disp(sprintf('%s',c_selected{ii,1}))
    if ii>1 && strcmp(c_selected{ii,2},c_selected{ii-1,2})
        % do nothing
    else
        load(c_selected{ii,2});
    end
    dat.rat(ii,1) = { sdata.rat_num };
    dat.date(ii,1) = { sdata.date };
    dat.uci(ii,1) = c_selected(ii,1);

    mname = 'square1_3d';        
    spx = sdata.(c_selected{ii,1}).(mname).spx;
    spy = sdata.(c_selected{ii,1}).(mname).spy;
    spt = sdata.(c_selected{ii,1}).(mname).spt;
    pox = sdata.(mname).pox;
    poy = sdata.(mname).poy;
    pot = sdata.(mname).pot;

    dat.pos(ii,1) = { [pox(:) poy(:) pot(:)] };
    dat.spk(ii,1) = { [spx(:) spy(:) spt(:)] };
end

save('example_pcells.mat','dat');



















if 0
    x = fieldnames(sdata);
    idx = contains(x,['r' sdata.rat_num]);
    idx2 = contains(x,'_c0');
    cnames = x(idx & ~idx2);
    
    
    fig1 = figure('visible','on','Position',[50,60,1400,800]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    
    xplots = 7;
    yplots = 7;
    max_cells = 100;
    ncells = length(cnames);
    tot_cells = min([ncells max_cells]);
    for cc = 1:tot_cells
        mname = 'square1_3d';    
        pox = sdata.(mname).pox;
        poy = sdata.(mname).poy;
        pot = sdata.(mname).pot;
        spx = sdata.(cnames{cc}).(mname).spx;
        spy = sdata.(cnames{cc}).(mname).spy;
        spt = sdata.(cnames{cc}).(mname).spt;
    
        rmset = struct; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
        rmset.method = 'histogram';
        rmset.binsize = 20;
        rmset.ssigma = 40;
        [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);   
    
        ax = subplot(xplots,yplots,cc);
        imagesc(ratemap);
        axis xy
        daspect([1 1 1])
        colormap(ax,turbo)
        title(sprintf('%s',cnames{cc}));
end
end
































