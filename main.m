%{Runge-Kutta Methods applied to Radioactivity}%

% Description: Using Rung Kutta methods to solve radiocarbon problem by
% plot it in different time steps. Besides, this program calculates the
% error (the difference between exact values and appoximate values)for
% each method. 

clc, clear all, close all; 

ti = 0;
tf =16;

dt = [1,0.5,0.25]; %array hold three different time steps

err_1= zeros(3,3);% "error" array for RK1
err_2= zeros(3,3);% "error" array for RK2
err_4= zeros(3,3);% "error" array for RK4

time = [5,10,15];

for i= 1:length(dt)
        
    nt = (tf - ti)/dt(i);%number of time steps
    t = linspace (ti, tf, nt);
    
    
    % initial "y" array for each method
    y1= zeros(1, nt);
    y2= zeros(1, nt);
    y4= zeros(1, nt);
    
    %initial condition
    y1(1) = 1;
    y2(1) = 1;
    y4(1) = 1;
    
    for k = 1:nt-1
        y1(k+1) =  RungeKutta(1,y1(k),dt(i));
        y2(k+1) =  RungeKutta(2,y2(k),dt(i));
        y4(k+1) =  RungeKutta(3,y4(k),dt(i));
            
    end
    
    %Calculate y exact by using given function
    y_exact = zeros (1,nt);
    for m = 1:nt
        y_exact(m) = exp(-2.45^(-1)*log(2)*t(m));
    end    
    
    %Using ode45
    tspan = [ti tf];
    y0 = 1;
    [t_o,y_o] = ode45(@(t_o,y_o) (-log(2)/2.45)*y_o, tspan, y0);
    
    
   
    % Plot for three time steps
    figure(i);
    p = plot (t,y1,'b',t,y2,'r',t,y4,'y', t_o, y_o,'-o',t,y_exact, 'g');
    xlabel ('time (s)');
    ylabel ('y');
    title (sprintf('RK1, RK2, RK4 with the time step %4.2f s',dt(i)));
    xlim ([0 16])
    set(p, 'LineWidth', 2);
    set(gcf, 'Position', [20 20 1200 600]);
    set(gca, 'LineWidth', 2, 'FontSize', 15);
    legend('RK1','RK2','RK4', 'ode45','exact');

    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on'; 
    
    
    
    %Calculate error
    %Loop through the time array holding three time values 5s, 10s, 15s
    for n = 1:length(time)
         ex = interp1(t,y_exact,time(n));
         y_rk1 = interp1(t,y1,time(n)); %interp1 function extracts value of y from the plot at given time t
         err_1(i,n)= abs(y_rk1 - ex);
         
         y_rk2 = interp1(t,y2,time(n));
         err_2(i,n)= abs(y_rk2 - ex);
         
         y_rk4 = interp1(t,y4,time(n));
         err_4(i,n)= abs(y_rk4 - ex);
         
         
    end  
  
end

%Display the output
outputResults(1, err_1, dt, 3);
outputResults(2, err_2, dt, 3);
outputResults(3, err_4, dt, 3);

function [ykp1] = RungeKutta(method,yk,dt)
% RUNGEKUTTA holds three methods RK1, RK2, RK4 that returns y at next time step
%   method is either 1,2 or 3 corresponding to RK1, RK2, RK4
%   input takes current y value and the time step size dt

    switch (method)
        case 1
        %% RK1
            c1 = dt*func(yk);
            ykp1 = yk + c1;
       
        case 2
         %% RK2
            c1 = dt*func(yk);
            c2 = dt*func(yk + 0.5*c1);
            ykp1 = yk + c2;
        
        case 3
        %% RK4
            c1 = dt*func(yk);
            c2 = dt*func(yk + 0.5*c1);
            c3 = dt*func(yk + 0.5*c2);
            c4 = dt*func(yk + c3);
            ykp1 = yk + (c1 + 2*c2 + 2*c3 + c4)/6;
  
    end   
      
end
function [out] = func(yk)
%FUNC is the function which calculate the amount of C-15 over time

    t =2.45; %half life
    out = (-log(2)/t)*yk;
end

function outputResults(method, err,dt, maxIt)
%OUTPUTRESULTS is to display the error output for each method RK
    
    switch (method)
        case 1
            %RK1
            
            fprintf("               ==RK1==\n");
            fprintf ("                  5s        10s       15s\n");
            for i = 1:maxIt
                fprintf(sprintf('dt = %4.2f s:',dt(i)));
               
                disp (err(i,:));
               
            end
         case 2
             %RK2

            fprintf("               ==RK2==\n");
            fprintf ("                  5s        10s       15s\n");
            for i = 1:maxIt
                fprintf(sprintf('dt = %4.2f s:',dt(i)));
               
                disp (err(i,:));
               
            end 
         case 3
             %RK4
            fprintf("               ==RK4==\n");
            fprintf ("                  5s        10s       15s\n");
            for i = 1:maxIt
                fprintf(sprintf('dt = %4.2f s:',dt(i)));
               
                disp (err(i,:));
               
            end      
    end        
     
end
