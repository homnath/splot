function delocation(fheader,model,rec,hdata,compute,delocate,iloop)
% Revision
%   HNG, Nov 11,2009
clear all
t_start=clock; % Record starting clock

if delocate.savefig>0
    addpath('./savefigure');
end

if isempty(hdata)
    error('hdata is empty!');
end

% Scale everything to model origin
rec.x(:,1)=rec.x(:,1)-model.ox(1);
rec.x(:,2)=rec.x(:,2)-model.ox(2);
rec.x(:,3)=rec.x(:,3)-model.ox(3);

delocate.range(1:3,1)=delocate.range(1:3,1)-model.ox;
delocate.range(1:3,2)=delocate.range(1:3,2)-model.ox;

delocate.init(1:3,1)=delocate.init(1:3,1)-model.ox;

model.ox(1:3)=0;

% Add path to folder de
addpath ./de

% Differential evolution global minimization
% set title
optimInfo.title = 'DE location';
optimInfo.fheader = fheader;

% specify objective function
objFctHandle = @deobject;

% define parameter names, ranges, quantizations and initial values :
paramDefCell = {'',delocate.range, delocate.quant, delocate.init};

% no additional function parameters needed
objFctSettings = {};

% no parameter vector needed
objFctParams = [];

% get default DE parameters
DEParams = getdefaultparams;

% set number of population members (often 10*D is suggested; here we use
% more as we know that the Foxholes functions has many local minima).
DEParams.NP = delocate.np; %50

DEParams.F = delocate.f;
DEParams.CR = delocate.c;
DEParams.strategy = delocate.type;

% do not use slave process here
DEParams.feedSlaveProc = 0;

% set times
DEParams.maxiter       = delocate.ngen;
DEParams.maxtime       = delocate.tmax;
DEParams.maxclock      = [];

% set display options
DEParams.refreshiter   = 1;
DEParams.refreshtime   = 10;  % in seconds
DEParams.refreshtime2  = 20;  % in seconds
DEParams.refreshtime3  = 40;  % in seconds
DEParams.minimizeValue = 0; % We want to maximize
DEParams.maxsds        = delocate.sds;
DEParams.maxsdt        = delocate.sdt;
DEParams.plotfig       = delocate.plotfig;
DEParams.savefig       = delocate.savefig;
if ~exist('../temp','dir')
    stat = mkdir('../temp');
    if ~stat
        error('directory ''../temp'' cannot be created!');
    end
end

outf_temp=strcat('../temp/detemp');
save(outf_temp,'fheader','model','rec','hdata','compute');

% set random state in order to always use the same population members here
rand('state', 1);
%deobject([490.1767; 187.5114; 305.202; 2.2])
%deobject([490.1767; 518.; 518.; 2.1])
% start differential evolution
[srcx,sd_srcx,sd_space,sd_time,stackval,ngen,nfun,gen_pop,gen_val,gen_sdpar,gen_meanval,gen_bestval,allmem,allval,minMaxSign] = ...
    differentialevolution(DEParams, paramDefCell, objFctHandle, ...
    objFctSettings, objFctParams,optimInfo);

t_total=etime(clock,t_start);
fprintf(1,'    Total elapsed time (s): %.4f\n',t_total);
fprintf(1,'    Elapsed time per generation (s): %.4f\n',t_total/ngen);

rmpath('./de');

if compute.nmonte>0
    outdat=sprintf('%s_optimization_result_monte_carlo%d',fheader,iloop); % .mat file
    outfname=strcat(fheader,'_optimization_result_summary_monte_carlo');
elseif compute.njack>0
    outdat=sprintf('%s_optimization_result_jack_knife%d',fheader,iloop); % .mat file
    outfname=strcat(fheader,'_optimization_result_summary_jack_knife');
else
    outdat=sprintf('%s_optimization_result',fheader); % .mat file
    outfname=strcat(fheader,'_optimization_result_summary');
end

save(outdat,'srcx','sd_srcx','sd_space','sd_time','stackval','ngen','nfun','gen_pop','gen_val','gen_sdpar','gen_meanval','gen_bestval','allmem','allval','minMaxSign');

last_pop=gen_pop{end};
mean_srcx=mean(last_pop);

% Uncertainty for 95% confidence interval
p=95;
ipop=allval >= 0.01*p*max(allval);
sd_p=std(allmem(ipop,:));
sd_ps=sqrt(sum(sd_p(1:3).^2));

% Write optimization result summary
if iloop==1
    outf=fopen(outfname,'w');
    fprintf(outf,'best_src(1:4), mean_srcx(1:4), sd_src(1:4), sd_space, sd_src%d%%, sd_space%d%%, ngen, nfun\n',round(p),round(p));
else
    outf=fopen(outfname,'a');
end
fprintf(outf,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d\n',srcx,mean_srcx,sd_srcx,sd_space,sd_p,sd_ps,ngen,nfun);
fclose(outf);
%---------------------------

fprintf(1,'95%% uncertainty: %f %f %f %f\n',sd_p);
fprintf(1,'95%% uncertainty sd_space: %f\n',sd_ps);

fprintf(1,'Standard deviation space-time: %f %f %f %f\n',sd_srcx);
fprintf(1,'Standard deviation space: %f\n',sd_space)

fprintf(1,'Best estimated source location: x=%f y=%f z=%f t0=%f\n',srcx(1),srcx(2),-srcx(3),srcx(4));
fprintf(1,'Mean estimated source location: x=%f y=%f z=%f t0=%f\n',mean_srcx(1),mean_srcx(2),-mean_srcx(3),mean_srcx(4));

if delocate.plotfig>0 || delocate.savefig>0
    % Code below is only for plotting the figure
    fig_prop1='visible';
    fig_val1='on';
    if delocate.plotfig==0
        fig_val1='off';
    end

    %blue = [0.0 0.4 1.0];
    %lblue = [0.65 0.85 1.0]; % light blue
    lred = [1.0 0.7 0.7]; % light red
    dred = [0.5 0.0 0.0];


    if delocate.savefig>0    
        outfres=strcat(fheader,'_result');
    end

    figure(fig_prop1,fig_val1)
    set(gca,'Fontsize',14)
    set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])
    %plot3(xr,yr,-zr,'^g','markersize',8,'markerfacecolor','g');
    %hold on
    %for i_rec=1:nrec
    %    text(xr(i_rec),yr(i_rec),-zr(i_rec),num2str(i_rec));
    %end
    %hold on
    %plot3(xs0(1),xs0(2),-xs0(3),'o','markersize',8,'MarkerFaceColor',blue,'MarkerEdgeColor',blue);
    plot3(last_pop(:,1),last_pop(:,2),-last_pop(:,3),'o','markersize',8,'MarkerFaceColor',lred,'MarkerEdgeColor',lred);
    plot3(mean_srcx(1),mean_srcx(2),-mean_srcx(3),'o','markersize',8,'markerFaceColor','r','MarkerEdgeColor','r');
    %plot3(estmem(:,1),estmem(:,2),-estmem(:,3),'oc','markersize',8,'markerfacecolor','c');    
    plot3(srcx(1),srcx(2),-srcx(3),'ob','markersize',8,'MarkerFaceColor',dred,'MarkerEdgeColor',dred);         

    % Set plot range to entire model
    %xmin=model.ox(1); xmax=model.ox(1)+(model.nx(1)-1)*model.dh;
    %ymin=model.ox(2); ymax=model.ox(2)+(model.nx(2)-1)*model.dh;
    %zmin=model.ox(3); zmax=model.ox(3)+(model.nx(3)-1)*model.dh;
    
    % Set plot range to DE range
    xmin=delocate.range(1,1); xmax=delocate.range(1,2);
    ymin=delocate.range(2,1); ymax=delocate.range(2,2);
    zmin=-delocate.range(3,2); zmax=-delocate.range(3,1); % Map to z positive up

    xlim=([xmin xmax]);
    ylim=([ymin ymax]);
    zlim=([zmin zmax]);

    xtick1=round10(xmin,2); xtick2=round10(xmax,2);
    ytick1=round10(ymin,2); ytick2=round10(ymax,2);
    ztick1=round10(zmin,2); ztick2=round10(zmax,2);

    set(gca,'box','on','xlim',xlim,'ylim',ylim,'zlim',zlim);
    set(gca,'xtick',xtick1:100:xtick2,'xticklabel',xtick1:100:xtick2);
    set(gca,'ytick',ytick1:100:ytick2,'yticklabel',ytick1:100:ytick2);
    set(gca,'ztick',ztick1:100:ztick2,'zticklabel',-ztick1:-100:-ztick2); % But we label as if z positive down

    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');   

    if delocate.savefig>0
        savefigure(outfres);
    end

    % Plot all members in scaled color
    rawval=minMaxSign*allval;
    
    % Get color map
    % Extrema (max or min) is in red.
    % Change only 2 lines below to flip the color
    % or to change the color for max and min e.g. max always red and min
    % always blue
    mp0=jet(256);
    mp=mp0(end:-1:1,:);
    % Change 3 lines below to display best members on top  
    %[rawval]=sort(rawval,'ascend');
    minval=min(rawval); % rawval(1);
    maxval=max(rawval); % rawval(end);

    figure(fig_prop1,fig_val1)
    set(gca,'Fontsize',14)
    set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])
    %plot3(xr,yr,-zr,'^g','markersize',8,'markerfacecolor','g');
    %hold on
    %for i_rec=1:nrec
    %    text(xr(i_rec),yr(i_rec),-zr(i_rec),num2str(i_rec));
    %end
    %hold on
    %plot3(srcx(1),srcx(2),-srcx(3),'o','markersize',8, 'Color','b','markerfacecolor','b');
    %hold on
    %plot3(xs0(1),xs0(2),-xs0(3),'or','markersize',8,'MarkerFaceColor',blue,'MarkerEdgeColor',blue);

    %scatter3(allmem(isort,1), allmem(isort,2), -allmem(isort,3),12,rawval,'o','filled');
    scatter3(allmem(:,1), allmem(:,2), -allmem(:,3),12,rawval,'o','filled');
    colormap(mp);
    caxis([minval maxval]);
    grid off
    set(gca,'box','on','xlim',xlim,'ylim',ylim,'zlim',zlim);
    set(gca,'xtick',xtick1:100:xtick2,'xticklabel',xtick1:100:xtick2);
    set(gca,'ytick',ytick1:100:ytick2,'yticklabel',ytick1:100:ytick2);
    set(gca,'ztick',ztick1:100:ztick2,'zticklabel',-ztick1:-100:-ztick2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    if delocate.savefig>0    
        outfres_all=strcat(fheader,'_result_all');
        savefigure(outfres_all);
    end
    if delocate.savefig>0
        rmpath('./savefigure');
    end
end % if delocate.plotfig>0 || delocate.savefig>0

if delocate.saveanim    
    % Save DE animation
    % Set plot range to DE range
    xmin=delocate.range(1,1); xmax=delocate.range(1,2);
    ymin=delocate.range(2,1); ymax=delocate.range(2,2);
    zmin=-delocate.range(3,2); zmax=-delocate.range(3,1); % Map to z positive up

    xlim=([xmin xmax]);
    ylim=([ymin ymax]);
    zlim=([zmin zmax]);

    xtick1=round10(xmin,2); xtick2=round10(xmax,2);
    ytick1=round10(ymin,2); ytick2=round10(ymax,2);
    ztick1=round10(zmin,2); ztick2=round10(zmax,2);
    
    rawval=minMaxSign*reshape(gen_val,[numel(gen_val) 1]); % Changed to one dimensional array 
    
    % This segment is necessary only for formatting of progress display
    ndig=ceil(log10(ngen+1));
    form=sprintf('saving DE animation...%%%dd/%%%dd',ndig,ndig);
    dums=sprintf(form,ngen,ngen);
    nbs=length(dums); % Number of backspaces required for overwriting
    bs(1:2:2*nbs-1)='\';
    bs(2:2:2*nbs)='b';
    nform=sprintf('%s%s',bs,form); % Format including backspaces
    %---------------------------
    
    % Get color map
    % Extrema (max or min) is in red.
    % Change only 2 lines below to flip the color
    % or to change the color for max and min e.g. max always red and min
    % always blue
    mp0=jet(256);
    mp=mp0(end:-1:1,:);
    % Change 3 lines below to display best members on top  
    %[rawval]=sort(rawval,'ascend');
    minval=min(rawval); % rawval(1);
    maxval=max(rawval); % rawval(end);
    
    for i_gen=1:ngen        
        figure('visible','on')
        set(gca,'Fontsize',12)
        set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])
        
        ind1=(i_gen-1)*delocate.np+1;
        pop_index=ind1:1:ind1+delocate.np-1;
        scatter3(gen_pop{i_gen}(:,1), gen_pop{i_gen}(:,2), -gen_pop{i_gen}(:,3),12,rawval(pop_index),'o','filled');
        colormap(mp);
        caxis([minval maxval]);
        grid off
        set(gca,'box','on','xlim',xlim,'ylim',ylim,'zlim',zlim);
        set(gca,'xtick',xtick1:100:xtick2,'xticklabel',xtick1:100:xtick2);
        set(gca,'ytick',ytick1:100:ytick2,'yticklabel',ytick1:100:ytick2);
        set(gca,'ztick',ztick1:100:ztick2,'zticklabel',-ztick1:-100:-ztick2);
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('z (m)');
        
        outf_anim=strcat(fheader,'_de_animation_',num2str(i_gen));
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
end % if delocate.saveanim
%------------------------------------------------------
