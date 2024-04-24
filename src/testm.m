
t=linspace(0,3,3000000);
x1=25;
x2=75;
fr= 50;

%P1= pressure(x1);
%P2= pressure(x2);
%disp(P1)
%disp(P2)



figure
plot(t,pressure(x1))
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('Original Data','Low Pass','Band Pass','High Pass')


function [P]= pressure(x)
t= linspace(0,3,3000000);
fr= 50;
v1=4;
P = (1./(4.*pi.*abs(x-x_s(t,v1))).*p_src(fr,t));
%disp(P)
end 

function xs = x_s(t,v)
%moving source function 
xs = v.*t+50;
end

function ps= p_src(f,t)
%pressure source function 
ps = 0.1.*sin(2.*pi.*f.*t);
end 
