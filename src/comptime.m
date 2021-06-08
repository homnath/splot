function comptime(niter)
for i_iter=1:niter
    t1=clock;
    roundn(123456789,7);%fprintf(1,'roundn(123456789,4): %d\n',roundn(123456789,7));
    t2=clock;
    tf1(i_iter)=etime(t2,t1);
    %fprintf(1,'Total time: %f\n',etime(t2,t1));

    t1=clock;
    round10(123456789,7); %fprintf(1,'round10(123456789,4): %d\n',round10(123456789,7));
    t2=clock;
    tf2(i_iter)=etime(t2,t1);
    %fprintf(1,'Total time: %f\n',etime(t2,t1));
end

fprintf(1,'Mean computation time: roundn=%f, round10=%f\n',mean(tf1),mean(tf2));
fprintf(1,'SD of computation time: roundn=%f, round10=%f\n',std(tf1),std(tf2));
figure;
hold on
plot(tf1,'-k');
plot(tf2,'-r');
legend('roundn','round10');
xlabel('computation number')
ylabel('computation time (s)')

a=rand(3011,1);
t0=clock;
abs(hilbert(a));
t1=clock;
abs(hilbert(a(1:100)));
t2=clock;
format long g
et1=etime(t1,t0)
et2=etime(t2,t1)
