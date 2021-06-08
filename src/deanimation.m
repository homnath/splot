% Save DE animation
function deanimation(infname)
% infname: DE optimization result file
    
    
% Need to be changed for general (changed)
[path,fheader]=fileparts(infname);
fheader=fullfile('../output/',fheader);

load(infname,'gen_pop','gen_val');
ngen=length(gen_pop);
NP=size(gen_val,1);
% get color map
mp      = jet(256);
%mp = mp0(end:-1:1,:); % Reverse the colormap (for red: max, and darkblue: min)
%mp(1,:) = [0 0 0.4]; % set dark blue for lowest value
Nmp     = size(mp, 1);
% scale cost values
a = 10;
p = 0.9;
rawval=reshape(gen_val,[numel(gen_val) 1]); % Changed to one dimensional array
maxVal = quantile2(rawval, p);
if maxVal > 0
    rawvalScaled = maxVal*log(a*rawval/maxVal+1)/log(a+1);
    index = rawval > maxVal;
    rawvalScaled(index) = maxVal;
else
    rawvalScaled = rawval;
end

% compute bins for cost value color indices
minValScaled = min(rawvalScaled);
maxValScaled = max(rawvalScaled);
if minValScaled ~= maxValScaled	     
    edges = linspace(minValScaled, maxValScaled, Nmp+1);
    edges(1)   = -inf;
    edges(end) = inf;
end

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

    ind1=(i_gen-1)*NP+1;
    pop_index=ind1:1:ind1+NP-1;

    % divide members into bins and plot
    if minValScaled ~= maxValScaled
        for n=length(edges)-1:-1:1
            ipop = find(rawvalScaled(pop_index) >= edges(n) & rawvalScaled(pop_index) < edges(n+1));                
            plot3(gen_pop{i_gen}(ipop,1), gen_pop{i_gen}(ipop,2), -gen_pop{i_gen}(ipop,3),'o','markersize',2, 'Color', mp(n,:),'markerfacecolor',mp(n,:));
            hold on
            pop_index=setdiff(pop_index,ipop);
            if isempty(pop_index)
                break;
            end
        end
    else
        plot3(gen_pop{i_gen}(:,1), gen_pop{i_gen}(:,2), -gen_pop{i_gen}(:,3),'o','markersize',2, 'Color', mp(1,:),'markerfacecolor',mp(1,:));
    end

    set(gca,'box','on','xlim',xlim,'ylim',ylim,'zlim',zlim);
    set(gca,'xtick',xtick1:100:xtick2,'xticklabel',xtick1:100:xtick2);
    set(gca,'ytick',ytick1:100:ytick2,'yticklabel',ytick1:100:ytick2);
    set(gca,'ztick',ztick1:100:ztick2,'zticklabel',-ztick1:-100:-ztick2);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    if delocate.savefig>0    
        outf_anim=strcat(fheader,'_animation_',num2str(i_gen));
        print('-dtiff','-r600',outf_anim);
    end
    % Display progress
    if (i_gen==1)
        fprintf(1,form,i_gen,ngen); % Print without backspaces
    else
        fprintf(1,nform,i_gen,ngen); % Print with backspaces
    end        
    %-----------------
end
fprintf(1,'\n');