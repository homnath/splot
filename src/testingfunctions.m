cx1=25;
x2=75;


P= pressure(x1);
P2= pressure(x2);
disp(P)
disp(P2)

%ps= p_src(fr,t);
%disp(ps

%VErsion 1 all in 1
function P= P1(x,f,t)
v=4;
x_s = v.*t+50;
ps = 0.1.*sin(2.*pi.*f.*t);
P = (1./(4.*pi.*abs(x-x_s)).*ps);

end

%version 2 defining xs and ps in P
function ps= p_src(f,t)
ps = 0.1.*sin(2.*pi.*f.*t);
end 

function xs = x_s(t)
v=4;
xs = v.*t+50;
end 

function P= P1(x,f,t)
p_src1=p_src(f,t);
x_sv= x_s(t);
P = (1./(4.*pi.*abs(x-x_s(t))).*p_src1);
end 



%version 3 each function seperate  
function [P]= pressure(x)
t= linspace(0,3,0.000001);
fr= 50;
v1=4;
P = (1./(4.*pi.*abs(x-x_s(t,v1))).*p_src(fr,t));
disp(P)
end 

function xs = x_s(t,v)
%moving source function 
xs = v.*t+50;
end

function ps= p_src(f,t)
%pressure source function 
ps = 0.1.*sin(2.*pi.*f.*t);
end 


%version 4 nested functions 

function [P]= pressure(x)
t= linspace(0,3,0.000001) ;
fr= 50;
v1=4;
P = (1./(4.*pi.*abs(x-x_s(t,v1))).*p_src(fr,t));
    function xs = x_s(t,v)
    xs = v.*t+50;
    end
    function ps= p_src(f,t)
    %pressure source function 
    ps = 0.1.*sin(2.*pi.*f.*t);
    end 
end





