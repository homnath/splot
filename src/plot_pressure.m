inpf=fopen('/home/hgharti/tigress/groningen/specfem3d/EXAMPLES/moving_source/OUTPUT_FILES/HPHONE.Y001X001.FXP.semp','r');
data=fscanf(inpf,'%f %f',[2 inf])';

figure
plot(data(:,1),data(:,2));

dt=data(2,1)-data(1,1)
datap=data(:,2);
forder=4; ftype=2; ffreq=[10 200];
datap=buttwo(forder,ffreq,dt,ftype,1,datap);

figure
plot(data(:,1),datap);