
function bien=denoise(se,qq,nume)

% Function for denoising the signal by wavelet packets

% se -> signal 
% qq -> sample of noise 
% nume -> maximum number of decomposition levels


% Decompose the signal on all the nodes of the wavelet packet tree
t=wpdec(se,nume,'db12','shannon');

% Obtain the best decomposition for this signal through the Shannon entropy
[t1,e,d1]=besttree(t);
nume=treedpth(t1);

% Copy the data structure
d3=d1;

% Decompose the sample of noise on all the nodes of the wavelet packet tree
t2=wpdec(qq,nume,'db12','shannon');

% Determine the terminal nodes associated with the signal decomposition
n=leaves(t1);

% Apply the same decomposition tree of the signal to the noise tree
t2=wpjoin(t2,n);

% Number of terminal nodes
nu=length(n);

% Put some variables to 0 
vario1=0;
vario2=0;

% Denoising loop
for i=1:nu

        % Obtain the coefficients associated with the original signal 
		coefa=read(t1,'data',n(i));

        % Coeficientes correspondientes al ruido
		coefr=read(t2,'data',n(i));

        % Obtain the variance of the signal coefficients of this node.
		vario1=std(coefa).^2;

        % Obtain the variance of the noise coefficients of this node.
		vario2=std(coefr).^2;
        
        % Determine the threshold
        thet=sqrt(vario2*2*log(length(coefa)));
              
        % Make a copy of the obtained coefficients for the signal
		coefb=coefa;

	   % Apply a Shoft-Thresholding  Denoising
		lok=length(coefb);
		for ty=1:lok
			if abs(coefb(ty))<thet
				coefb(ty)=0;
            else
            % SHOFT THRESHOLDING
            coefb(ty)=coefb(ty)-thet*sign(coefb(ty));
            
            % HARD THRESHOLDING
            %coefb(ty)=coefb(ty);
            
            % NON-NEGATIVE GARROTE
            %coefb(ty)=coefb(ty)-((thet/coefb(ty))^2);
            
            % NON-NEGATIVE GARROTE MODIFICADO
            %coefb(ty)=coefb(ty)-((thet/coefb(ty))^2)*sign(coefb(ty));

			end
		end	
			
        % Modify the coefficients associated with currently used node by the new ones,
        % the filtered coefficients.
        t1=write(t1,'data',n(i),coefb);

% Loop end
end

% Reconstruction of the denoised signal
%plot(t1);

bien=wprec(t1);

% Plot the new signal
%figure;
%gcf4=axes('units','normalized','pos',[0.11 0.15 0.84 0.75],'Fontsize',8);

%%%%NO%%%%plot(ti,bien,'k','LineWidth',0.5);grid on;

%plot(bien,'k','LineWidth',0.5);grid on;

%%%%NO%%%%%%%%axis([0,10,-9,9]);

%title('Filtered signal with wavelet packet (Node dependent threshold)','FontSize',10);

%%%%NO%%%%%%%55%xlabel('Time (s)','FontSize',8);

%xlabel('Samples','FontSize',8);
%ylabel('Amplitude (levels)','FontSize',8);
%set(gcf,'position',[1 29 800 504]);
