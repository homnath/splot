% Marker size for plot3 is a point
% Marker aize for scatter3 is a square point, e.g., 8 in plot3 is equivalent to 64 in scatter3.
% However, with savefigure I get the larger size for scatter3 figure, therefore I use 16 for scatter 3 and 8 for plot3.
% On the toher hand for normal saving I get the same size.
function plot_deresult(preinfo,deplot)
% The function plot_deresult plot the various figures for DE result stored
% in a mat file.
% infname: input file name
% HISTORY: HNG, Dec 02, 2009

% Code below is only for plotting the figure
fig_prop1='visible';
fig_val1='on';
if preinfo.plot==0
    fig_val1='off';
end

% Initialize
gen_pop=[];
gen_val=[];
allmem=[];
srcx=[];
load(deplot.mfile);

% Set plot range to DE range
xmin=deplot.range(1,1); xmax=deplot.range(1,2);
ymin=deplot.range(2,1); ymax=deplot.range(2,2);
zmin=-deplot.range(3,2); zmax=-deplot.range(3,1); % Map to z positive up
tmin=deplot.range(4,1); tmax=deplot.range(4,2);

axlim(1,:)=([xmin xmax]);
axlim(2,:)=([ymin ymax]);
axlim(3,:)=([zmin zmax]);
axlim(4,:)=([tmin tmax]);

xtick1=round10(xmin,deplot.rtick(1)); xtick2=round10(xmax,deplot.rtick(1));
ytick1=round10(ymin,deplot.rtick(2)); ytick2=round10(ymax,deplot.rtick(2));
ztick1=round10(zmin,deplot.rtick(3)); ztick2=round10(zmax,deplot.rtick(3));
ttick1=round10(tmin,deplot.rtick(4)); ttick2=round10(tmax,deplot.rtick(4));

axtick_loc{1}=linspace(xtick1,xtick2, deplot.ntick(1));
axtick_loc{2}=linspace(ytick1,ytick2, deplot.ntick(2));
axtick_loc{3}=linspace(ztick1,ztick2, deplot.ntick(3));
axtick_loc{4}=linspace(ttick1,ttick2, deplot.ntick(4));

axtick_lab{1}=axtick_loc{1};
axtick_lab{2}=axtick_loc{2};
axtick_lab{3}=-axtick_loc{3};  % Map to z positive down
axtick_lab{4}=axtick_loc{4};

axlab{1}=sprintf('x (%s)',strtrim(deplot.lunit));
axlab{2}=sprintf('y (%s)',strtrim(deplot.lunit));
axlab{3}=sprintf('z (%s)',strtrim(deplot.lunit));
axlab{4}=sprintf('t_0 (s)');

deplot.npop=size(gen_val,1);
npar=size(deplot.range, 1);

% All coordinates mapped to positive Z-axis up only for plotting purpose
% However, we label axes as if Z axis is positive down!
srcx(3)=-srcx(3);
allmem(:,3)=-allmem(:,3);
for i_gen=1:ngen
    gen_pop{i_gen}(:,3)=-gen_pop{i_gen}(:,3);
end

% Plot last generation population
if deplot.lastpop>0
    fprintf(1,'plotting last generation population...');
    last_pop=gen_pop{end};
    mean_srcx=mean(last_pop);
    
    lred = [1.0 0.7 0.7]; % light red
    dred = [0.5 0.0 0.0];
    
    figure(fig_prop1,fig_val1)    
    set(gca,'Fontsize',14)
    set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])
    %plot3(xr,yr,-zr,'^g','markersize',6,'markerfacecolor','g');
    %hold on
    %for i_rec=1:nrec
    %    text(xr(i_rec),yr(i_rec),-zr(i_rec),num2str(i_rec));
    %end
    %hold on
    %plot3(xs0(1),xs0(2),-xs0(3),'o','markersize',8,'MarkerFaceColor',blue,'MarkerEdgeColor',blue);
    plot3(last_pop(:,1),last_pop(:,2),last_pop(:,3),'o','markersize',8,'MarkerFaceColor',lred,'MarkerEdgeColor',lred);
    hold on
    plot3(mean_srcx(1),mean_srcx(2),mean_srcx(3),'o','markersize',8,'markerFaceColor','r','MarkerEdgeColor','r');
    hold on  
    plot3(srcx(1),srcx(2),srcx(3),'ob','markersize',8,'MarkerFaceColor',dred,'MarkerEdgeColor',dred);
    
    set(gca,'box','on','xlim',axlim(1,:),'ylim',axlim(2,:),'zlim',axlim(3,:));
    set(gca,'xtick',axtick_loc{1},'xticklabel',axtick_lab{1});
    set(gca,'ytick',axtick_loc{2},'yticklabel',axtick_lab{2});
    set(gca,'ztick',axtick_loc{3},'zticklabel',axtick_lab{3});    

    xlabel(axlab{1});
    ylabel(axlab{2});
    zlabel(axlab{3});
    if preinfo.save>0
        outf=strcat(preinfo.fighead,'_last_population');
        savefigure(outf);
    end
    fprintf(1,'complete\n');
end

% Plot all populations
if deplot.allpop>0
    fprintf(1,'plotting all population...');
    % Color by allval
    minval=min(allval);
    maxval=max(allval);
    figure(fig_prop1,fig_val1)
    set(gca,'Fontsize',14)
    set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])   
    
    plot3_color(allmem(:,1),allmem(:,2),allmem(:,3),allval,[minval maxval],{'o',8})
    
    set(gca,'box','on','xlim',axlim(1,:),'ylim',axlim(2,:),'zlim',axlim(3,:));
    set(gca,'xtick',axtick_loc{1},'xticklabel',axtick_lab{1});
    set(gca,'ytick',axtick_loc{2},'yticklabel',axtick_lab{2});
    set(gca,'ztick',axtick_loc{3},'zticklabel',axtick_lab{3});    

    xlabel(axlab{1});
    ylabel(axlab{2});
    zlabel(axlab{3});
    if preinfo.save>0    
        outf=strcat(preinfo.fighead,'_all_population');
        savefigure(outf);
    end
    fprintf(1,'complete\n')
end % 

% Plot generation vs objective function
% We shift the stair plot add additional value to put the points in the
% middle
if deplot.gen_of>0
    fprintf(1,'plotting generation vs objective function...');    
    fig_genof=figure(fig_prop1,fig_val1); %fig_genof=figure(fig_prop1,fig_val1);
    hold on
    hvect=1:ngen;
    stairs([0.5 hvect+0.5],[gen_bestval(1,1:ngen) gen_bestval(1,ngen)],'-k');
    child_handles=findobj(fig_genof,'type','line');
    leg(1)=child_handles(1);
    plot(hvect,gen_bestval(1,1:ngen),'ok','markersize',5,'MarkerFaceColor','k');
    %child_handles=findobj(fig_genof,'type','line');
    %leg(2)=child_handles(2);
    
    stairs([0.5 hvect+0.5],[gen_meanval(1,1:ngen) gen_meanval(1,ngen)],'-r');
    child_handles=findobj(fig_genof,'type','line');
    leg(2)=child_handles(2);
    plot(hvect,gen_meanval(1,1:ngen),'or','markersize',5,'MarkerFaceColor','r');
    %child_handles=findobj(fig_genof,'type','line');
    %leg(4)=child_handles(4);
    
    legend(leg,{'Best','Mean'},'Location','SouthEast');
    legend(gca,'boxoff')
    xlabel('Generation number');
    ylabel('S({\bfx})');
    if preinfo.save>0
        outf=sprintf('%s_objective_evolution',preinfo.fighead);
        savefigure(outf);
    end
    fprintf(1,'complete\n')
end

% Plot generation vs standard deviation
if deplot.gen_sd>0
    fprintf(1,'plotting generation vs standrad deviation...');
    fig_gensd=figure(fig_prop1,fig_val1); %fig_gensd=figure(fig_prop1,fig_val1);
    hold on
    hvect=1:ngen;
    stairs([0.5 hvect+0.5],[gen_sdpar(1,1:ngen) gen_sdpar(1,ngen)],'-k');
    child_handles=findobj(fig_gensd,'type','line');
    leg(1)=child_handles(1);
    plot(hvect,gen_sdpar(1,1:ngen),'ok','markersize',5,'MarkerFaceColor','k');
    %child_handles=findobj(fig_gensd,'type','line');
    %leg(2)=child_handles(2);
    
    stairs([0.5 hvect+0.5],[gen_sdpar(2,1:ngen) gen_sdpar(2,ngen)],'-r');
    child_handles=findobj(fig_gensd,'type','line');
    leg(2)=child_handles(2);
    plot(hvect,gen_sdpar(2,1:ngen),'or','markersize',5,'MarkerFaceColor','r');
    %child_handles=findobj(fig_gensd,'type','line');
    %leg(4)=child_handles(4);
    
    
    stairs([0.5 hvect+0.5],[gen_sdpar(3,1:ngen) gen_sdpar(3,ngen)],'-g');
    child_handles=findobj(fig_gensd,'type','line');
    leg(3)=child_handles(3);
    plot(hvect,gen_sdpar(3,1:ngen),'og','markersize',5,'MarkerFaceColor','g');
    %child_handles=findobj(fig_gensd,'type','line');
    %leg(6)=child_handles(6);
    
    stairs([0.5 hvect+0.5],[gen_sdpar(4,1:ngen) gen_sdpar(4,ngen)],'-b');
    child_handles=findobj(fig_gensd,'type','line');
    leg(4)=child_handles(4);
    plot(hvect,gen_sdpar(4,1:ngen),'ob','markersize',5,'MarkerFaceColor','b');
    %child_handles=findobj(fig_gensd,'type','line');
    %leg(8)=child_handles(8);
    %plot(hvect,stdval,'-om','MarkerFaceColor','m');
    xlabel('Generation number');
    ylabel('Normalized standard deviation');    
    legend(leg,{'x','y','z','t_o'});
    legend(gca,'boxoff')
    if preinfo.save>0
        outf=sprintf('%s_std_evolution',preinfo.fighead);
        savefigure(outf);
    end
    fprintf(1,'complete\n')
end

% Plot parameter vs objective function
if deplot.par_of>0
    fprintf(1,'plotting parameter vs objective function...');
    % plot one parameter into one figure    
    for i_par = 1:npar        
        [allmemSorted, sortIndex] = sort(allmem(:,i_par));
        figure(fig_prop1,fig_val1)
        plot(allmemSorted, allval(sortIndex),'o','markersize',8,'MarkerFaceColor','b');

        set(gca,'xlim',axlim(i_par,:));
        set(gca,'xtick',axtick_loc{i_par},'xticklabel',axtick_lab{i_par});
        xlabel(axlab{i_par});        
        ylabel('S({\bfX})');        

        if preinfo.save>0
            outf=sprintf('%s_all_parameter%d',preinfo.fighead,i_par);
            savefigure(outf);
        end
    end
    fprintf(1,'complete\n')
end

% Plot parameter vs generation
if deplot.par_gen>0
    fprintf(1,'plotting parameter vs generation...');    
    % Color by genval
    minval=min(min(gen_val));
    maxval=max(max(gen_val));
    for i_par = 1:npar
        gvect=zeros(deplot.npop,1);
        figure(fig_prop1,fig_val1)
        hold on
        for i_gen=1:ngen
            gvect(:)=i_gen;                                      
            plot_color(gen_pop{i_gen}(:,i_par),gvect,gen_val(:,i_gen),[minval maxval],{'o',8});    
        end 
        
        yt=get(gca,'ytick');
        yt(1)=1; yt(end)=ngen;
        set(gca,'ytick',yt,'ylim',[1 ngen],'ydir','reverse');
        set(gca,'xlim',axlim(i_par,:));
        set(gca,'xtick',axtick_loc{i_par},'xticklabel',axtick_lab{i_par});
        xlabel(axlab{i_par});
        ylabel('Generation number');        
        if preinfo.save>0
            outf=sprintf('%s_generation_parameter%d',preinfo.fighead,i_par);
            savefigure(outf);
        end
    end
    fprintf(1,'complete\n')
end

% Plot parameter vs parameter
if deplot.par_par>0
    fprintf(1,'plotting parameter vs parameter...');
    % Color by allval    
    minval=min(allval);
    maxval=max(allval);    

    % plot 2 parameters into one figure    
    for i_par1=1:npar-1
        for i_par2=i_par1+1:npar           

            figure(fig_prop1,fig_val1)
            axes('FontSize', 14);
            
            plot_color(allmem(:,i_par1), allmem(:,i_par2),allval,[minval maxval],{'o',8});         

            xlabel(axlab{i_par1});
            ylabel(axlab{i_par2});      

            set(gca,'xlim',axlim(i_par1,:),'ylim',axlim(i_par2,:));
            set(gca,'xtick',axtick_loc{i_par1},'xticklabel',axtick_lab{i_par1});
            set(gca,'ytick',axtick_loc{i_par2},'yticklabel',axtick_lab{i_par2});            

            if preinfo.save>0
                outf=sprintf('%s_parameter%d_parameter%d',preinfo.fighead,i_par1,i_par2);
                savefigure(outf);
            end

        end
    end
    fprintf(1,'complete\n')
end

% Save DE animation 
if deplot.anim>0
    % Color by gen_val
    minval=min(min(gen_val));
    maxval=max(max(gen_val));
    
    % This segment is necessary only for formatting of progress display
    ndig=ceil(log10(ngen+1));
    form=sprintf('saving DE animation...%%%dd/%%%dd',ndig,ndig);
    dums=sprintf(form,ngen,ngen);
    nbs=length(dums); % Number of backspaces required for overwriting
    bs(1:2:2*nbs-1)='\';
    bs(2:2:2*nbs)='b';
    nform=sprintf('%s%s',bs,form); % Format including backspaces
    %---------------------------    
    
    for i_gen=1:ngen        
        figure('visible','off')
        set(gca,'Fontsize',12)
        set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])  
        
        plot3_color(gen_pop{i_gen}(:,1),gen_pop{i_gen}(:,2),gen_pop{i_gen}(:,3),gen_val(:,i_gen),[minval maxval],{'o',8})
                
        set(gca,'box','on','xlim',axlim(1,:),'ylim',axlim(2,:),'zlim',axlim(3,:));
        set(gca,'xtick',axtick_loc{1},'xticklabel',axtick_lab{1});
        set(gca,'ytick',axtick_loc{2},'yticklabel',axtick_lab{2});
        set(gca,'ztick',axtick_loc{3},'zticklabel',axtick_lab{3});    

        xlabel(axlab{1});
        ylabel(axlab{2});
        zlabel(axlab{3});
        
        outf_anim=strcat(preinfo.fighead,'_de_animation_',num2str(i_gen));
        print('-dtiff','-r600',outf_anim);
        
        % Display progress
        if (i_gen==1)
            fprintf(1,form,i_gen,ngen); % Print without backspaces
        else
            fprintf(1,nform,i_gen,ngen); % Print with backspaces
        end        
        %-----------------
    end
    fprintf(1,'\n');
end % if deplot.anim

fprintf(1,'-----------------------------------\n');
%---------------------------------------------------

    