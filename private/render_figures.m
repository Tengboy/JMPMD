function audit=render_figures(spec,out,fontName)
% Publication exports only: no estimator or frozen data is modified.
colors=[0.000 0.337 0.620;0.800 0.270 0.090;0.220 0.520 0.320;0.120 0.120 0.120;0.490 0.300 0.630;0.680 0.200 0.290;0.320 0.450 0.700];
markers={'o','s','^','none','d','v','x'}; styles={'-','--','-.',':','--','-.',':'};
audit=struct('figures',numel(spec),'curves',0,'points',0,'titles_empty',true,'entries',[]);
entries=cell(numel(spec),1);
for fi=1:numel(spec)
    f=spec(fi); rows=ceil(numel(f.panels)/f.cols);
    fig=figure('Visible','off','Color','w','Units','centimeters',...
        'Position',[2 2 f.width f.height],'PaperUnits','centimeters',...
        'PaperPosition',[0 0 f.width f.height],'PaperSize',[f.width f.height],...
        'Renderer','painters','InvertHardcopy','off');
    marginLeft=1.55; marginRight=.32; gapX=1.45;
    if f.cols==1, marginLeft=2; marginRight=.65; end
    if f.cols==3, marginLeft=1.40; gapX=1.30; end
    marginBottom=1.85; marginTop=.65; gapY=2.05;
    aw=(f.width-marginLeft-marginRight-gapX*(f.cols-1))/f.cols;
    ah=(f.height-marginBottom-marginTop-gapY*(rows-1))/rows;
    nf=0; np=0; legendAudit=cell(numel(f.panels),1);
    for pi=1:numel(f.panels)
        p=f.panels(pi); ri=floor((pi-1)/f.cols); ci=mod(pi-1,f.cols);
        pos=[marginLeft+ci*(aw+gapX),marginBottom+(rows-1-ri)*(ah+gapY),aw,ah];
        ax=axes(fig,'Units','centimeters','Position',pos,'FontName','Times New Roman',...
            'FontSize',9,'LineWidth',.65,'TickDir','in','Box','on',...
            'XColor',[.1 .1 .1],'YColor',[.1 .1 .1],'Layer','top');
        hold(ax,'on'); grid(ax,'on');
        ax.GridColor=[.68 .68 .68];ax.GridAlpha=.40;ax.GridLineStyle='-';
        ax.MinorGridColor=[.76 .76 .76];ax.MinorGridAlpha=.25;ax.MinorGridLineStyle=':';
        if p.log,ax.YScale='log';ax.YMinorGrid='on';end
        hh=gobjects(numel(p.curves),1); labels=cell(numel(p.curves),1);
        for j=1:numel(p.curves)
            c=p.curves(j); idx=c.style;
            assert(all(isfinite(c.x))&&all(isfinite(c.y)),'Nonfinite figure input.');
            if p.log,assert(all(c.y>0),'Nonpositive logarithmic data.');end
            args={'Color',colors(idx,:),'LineStyle',styles{idx},'LineWidth',1.2,...
                'Marker',markers{idx},'MarkerSize',4,'MarkerFaceColor','w'};
            if isempty(c.error)
                hh(j)=plot(ax,c.x,c.y,args{:});
                if numel(c.x)>30,hh(j).MarkerIndices=unique([1 10:10:numel(c.x) numel(c.x)]);end
            else
                hh(j)=errorbar(ax,c.x,c.y,c.error,args{:},'CapSize',4);
            end
            assert(isequal(hh(j).XData(:),c.x(:))&&isequal(hh(j).YData(:),c.y(:)),...
                'Exported line differs from frozen summary.');
            labels{j}=c.label;nf=nf+1;np=np+numel(c.x);
        end
        if ~isempty(p.xticks),xticks(ax,p.xticks);end
        if ~isempty(p.xlim),xlim(ax,p.xlim);end
        if ~isempty(p.ylim)
            ylim(ax,p.ylim);
        else
            ly=ylim(ax);
            if p.log
                span=log10(ly(2)/ly(1));ylim(ax,ly.*[10^(-.035*span) 10^(.035*span)]);
            else
                span=diff(ly);ylim(ax,ly+[-.025 .04]*span);
            end
        end
        if contains(p.ylabel,'差值')
            yline(ax,0,':','Color',[.35 .35 .35],'HandleVisibility','off','LineWidth',.7);
        end
        if ~isempty(p.vline)
            xline(ax,p.vline(1),':','Color',[.5 .5 .5],'HandleVisibility','off','LineWidth',.8);
        end
        xlabel(ax,p.xlabel,'FontName',fontName,'FontSize',10,'Interpreter','tex');
        ylabel(ax,p.ylabel,'FontName',fontName,'FontSize',10,'Interpreter','tex');
        title(ax,'');
        leg=legend(ax,hh,labels,'Location',p.legend,'FontName',fontName,'FontSize',7.5,...
            'Interpreter','none','Box','on','Color','w','EdgeColor',[.25 .25 .25]);
        leg.ItemTokenSize=[17 10];leg.LineWidth=.55;
        legendAudit{pi}=pad_legend(ax,leg,labels,hh);
        if numel(f.panels)>1
            text(ax,.5,-1.35/ah,['(' char('a'+pi-1) ')'],'Units','normalized',...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontName','Times New Roman','FontSize',9,'Clipping','off');
        end
        assert(isempty(ax.Title.String),'Plot title is prohibited.');
    end
    drawnow;
    stem=fullfile(out,f.name);
    print(fig,[stem '.pdf'],'-dpdf','-painters');
    print(fig,[stem '.png'],'-dpng','-r300');
    savefig(fig,[stem '.fig']);
    close(fig);
    entries{fi}=struct('name',f.name,'panels',numel(f.panels),'curves',nf,'points',np,...
        'data_exact',true,'titles_empty',true,'panel_labels_below',true,'formats',{{'pdf','png','fig'}},...
        'width_cm',f.width,'height_cm',f.height,'legends',[legendAudit{:}]);
    audit.curves=audit.curves+nf;audit.points=audit.points+np;
    fprintf('%02d/%02d %s: %d curves, %d points\n',fi,numel(spec),f.name,nf,np);
end

audit.entries=[entries{:}];
fid=fopen(fullfile(out,'figure_render_audit.json'),'w','n','UTF-8');
fprintf(fid,'%s\n',jsonencode(audit));fclose(fid);
fprintf('PASS: all %d figures exported without titles; %d curves/%d points exact.\n',audit.figures,audit.curves,audit.points);
end

function audit=pad_legend(ax,leg,labels,hh)
% Measure the actual CJK labels and reserve physical padding for both exports.
% Automatic TeX legend sizing left too little room at the right border.
drawnow;leg.Units='centimeters';old=leg.Position;ap=ax.Position;
widths=zeros(numel(labels),1);heights=widths;
for k=1:numel(labels)
    t=text(ax,0,0,labels{k},'Units','centimeters','FontName',leg.FontName,...
        'FontSize',leg.FontSize,'Interpreter','none','Visible','off');
    ex=t.Extent;widths(k)=ex(3);heights(k)=ex(4);delete(t);
end
token=leg.ItemTokenSize(1)*2.54/72;gap=.12;pad=.22;
requiredWidth=max(widths)+token+gap+2*pad;
requiredHeight=numel(labels)*max(heights)+.24;
sz=max(old(3:4)+[.28 .16],[requiredWidth requiredHeight]);
assert(sz(1)<ap(3)&&sz(2)<ap(4),'Legend does not fit its panel.');
% Re-evaluate placement using the padded box, including interpolated segments.
% MATLAB's earlier 'best' placement only knew the smaller automatic box.
xy=[];
for k=1:numel(hh)
    x=hh(k).XData(:);y=hh(k).YData(:);
    if strcmp(ax.YScale,'log'),y=log10(y);yl=log10(ax.YLim);else,yl=ax.YLim;end
    x=(x-ax.XLim(1))/diff(ax.XLim)*ap(3)+ap(1);
    y=(y-yl(1))/diff(yl)*ap(4)+ap(2);
    for j=1:numel(x)-1
        xy=[xy;linspace(x(j),x(j+1),101)' linspace(y(j),y(j+1),101)']; %#ok<AGROW>
    end
    if isa(hh(k),'matlab.graphics.chart.primitive.ErrorBar')
        lo=hh(k).YData(:)-hh(k).YNegativeDelta(:);
        hi=hh(k).YData(:)+hh(k).YPositiveDelta(:);
        for j=1:numel(x)
            ey=linspace(lo(j),hi(j),101)';ey=(ey-yl(1))/diff(yl)*ap(4)+ap(2);
            xy=[xy;repmat(x(j),101,1) ey]; %#ok<AGROW>
        end
    end
end
right=ap(1)+ap(3)-sz(1)-.10;left=ap(1)+.10;
top=ap(2)+ap(4)-sz(2)-.10;bottom=ap(2)+.10;
candidates=[right top;left bottom;left top;right bottom;(left+right)/2 top;(left+right)/2 bottom];
hits=zeros(size(candidates,1),1);
for k=1:size(candidates,1)
    b=candidates(k,:);
    hits(k)=nnz(xy(:,1)>=b(1)-.045&xy(:,1)<=b(1)+sz(1)+.045&...
        xy(:,2)>=b(2)-.045&xy(:,2)<=b(2)+sz(2)+.045);
end
[bestHits,best]=min(hits);new=[candidates(best,:) sz];
leg.Position=new;drawnow;
actual=leg.Position;
assert(actual(3)>=requiredWidth-1e-6&&actual(4)>=requiredHeight-1e-6);
audit=struct('labels',{labels},'interpreter','none','automatic_box_cm',old,...
    'padded_box_cm',actual,'max_text_width_cm',max(widths),...
    'required_width_cm',requiredWidth,'required_height_cm',requiredHeight,...
    'padding_each_side_cm',pad,'measured_fit',true,'sampled_curve_overlap',bestHits);
end
