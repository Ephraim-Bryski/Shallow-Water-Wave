% Ephraim Bryski
% CIV 524
% Final Project
% This code includes all simulations discussed in the report. (Some figures
% however were created in the command window). Change the variable "demo"
% to change the type of simulation.

clearvars
close all


write_vid=false; % whether it records and saves a video, directory must be changed


% choices for demo: reflect, shoal, per, superpos, pebble, pebble per (not gonna be used), diffract, refract
demo='refract';

g=10;
T=.5;
diffract_gap=0.2; % fraction of length open for waves
refract_prof=@(x,y,M,N) 1-x/M; % used to generate profile and to get depth at specific points to check refraction

% directory MUST be changed to record video:
if write_vid
    vid=VideoWriter(['G:\My Drive\CIV 524\Final Project\',demo,'.mp4']);
    open(vid)
end

% set up simulation conditions:
switch demo
    case {'reflect','per','shoal','superpos'}
        M=100;       % region x distance
        N=.2;        % region y distance
        dx=10^-1;    % position step
        dt=3*10^-3;  % time step
        t_tot=100;   % total simulation time
    otherwise
        M=10;
        N=10;
        dx=.1;
        dt=0.001;
        t_tot=25;
end
time_gap=50; % no effect on simulation, number of steps between replotting, make larger to make faster

% selects boundary conditions:
switch demo
    case 'pebble per'
        BC_x='fixed';
        BC_y='per';
    otherwise
        BC_x='fixed';
        BC_y='fixed';
end

% set up discretization of region and time:
m=M/dx;             
n=N/dx;
n_t_steps=t_tot/dt;

u=zeros(n,m+1);
v=zeros(n+1,m);
[X,Y] = meshgrid(linspace(0,M,m),linspace(0,N,n));
Z=zeros(n,m);

% set up bathymetry:
switch demo
    case 'refract'
        h0=refract_prof(X,Y,M,N);
    case 'shoal'
        h0=(M-0.8*X)/M*2.5;
    otherwise
        h0=2.5;
end

% set up visualization:
hold on
switch demo
    case {'reflect','per','shoal','superpos'}
        surface=plot3(linspace(0,M,m),zeros(1,m),zeros(1,m));
        view([0,-1,0])
        bounds=[-1,3];
        xlabel('Position (m)')
        zlabel('\eta (m)')
    case 'diffract'
        cla
        plot3([0,0,M,M,0,0],[N/2*(1+diffract_gap),N,N,0,0,N/2*(1-diffract_gap)],2+zeros(1,6),'k','LineWidth',5)
        hold on
        surface=surf(X,Y,Z);
        surface.EdgeColor='none';
        xlabel('X Position (m)')
        ylabel('Y Position (m)')
        bounds=[-1.8,2.8];
        view(2)
        c=colorbar;
        c.Label.String="Water Surface (m)";
        hold off
    case 'refract'
        surface=surf(X,Y,Z);
        plot3([0,M,M,0],[N,N,0,0],1+zeros(1,4),'k','LineWidth',5)
        surface.EdgeColor='none';
         xlabel('X Position (m)')
        ylabel('Y Position (m)')
        bounds=[-.5,1];
        view([.01,-1,0])
        c=colorbar;
        c.Label.String="Water Surface (m)";
%         camroll(-45)
    case {'pebble','pebble per'} 
        surface=surf(X,Y,Z);
        surface.EdgeColor='none';
        xlabel('X Position (m)')
        ylabel('Y Position (m)')
        zlabel('\eta (m)')
        bounds=[-.5,1];
        view(3)
end
% camlight('right')
zlim(bounds)
caxis(bounds)
axis square



% set up initial water profile to produce waves:
switch demo
    case {'reflect','per','shoal'}
        etta=2*(exp(-(0.1*(X)).^2));
    case 'superpos'
        etta=2*(exp(-(0.2*(X)).^2))+exp(-(0.1*(X-M)).^2);
    case {'diffract','refract'}
        etta=0*X; % waves are instead produced by having a boundary condition with a varying velocity
    case {'pebble','pebble per'}
        etta=2*(exp(-(0.2*(X-2*M/3)).^2-(0.2*(Y-3*M/4)).^2));
end


for t_step=1:n_t_steps 

    t=t_step*dt;

    % set up the bounding velocity:
    switch demo
        case 'diffract'
            v_forcing=[zeros((1-diffract_gap)/2*n,1);
                    8*sin(2*pi/T*t)*ones(diffract_gap*n,1);
                    zeros((1-diffract_gap)/2*n,1)];
        case 'refract'
            if t<T
                v_forcing=[1*sin(2*pi/T*t)*ones(n,1)];
            else
                v_forcing=zeros(n,1);
            end
        otherwise
            v_forcing=zeros(n,1);
    end

    
    h=h0;

    % for the periodic demo, the BC is initially fixed to produce the wave
    % and is switched to periodic when it travels away from the wall:
    if t_step==2000 && strcmp(demo,'per')
        BC_x='per';
    end
        

    % discretization of du/dt=-g*dn/dx:
    if strcmp(BC_x,'fixed') 
        dettadx=(etta(:,2:end)-etta(:,1:end-1))/dx;
        dudt=[v_forcing,-g*dettadx,zeros(n,1)];
    elseif strcmp(BC_x,'per')
        etta_x=[etta,etta(:,1)];
        dettadx=(etta_x(:,2:end)-etta_x(:,1:end-1))/dx;
        dudt=-g*dettadx;
        dudt=[dudt(:,end),dudt];
    else
        error('Unknwon BC_x name')
    end

    % discretization of dv/dt=-g*dn/dy:
    if strcmp(BC_y,'fixed')   
        dettady=(etta(2:end,:)-etta(1:end-1,:))/dx;
        dvdt=[zeros(1,m);-g*dettady;zeros(1,m)];
    elseif strcmp(BC_y,'per')
        etta_y=[etta;etta(1,:)];
        dettady=(etta_y(2:end,:)-etta_y(1:end-1,:))/dx;
        dvdt=-g*dettady;
        dvdt=[dvdt(end,:);dvdt];
    else
        error('Unknwon BC_y name')
    end
    
    % discretization of dh/dt=-h0*(du/dx+dv/dy):
    u=u+dudt*dt;
    v=v+dvdt*dt;
    h_ext_x=[h(:,1),(h(:,1:end-1)+h(:,2:end))/2,h(:,end)];
    h_ext_y=[h(1,:);(h(1:end-1,:)+h(2:end,:))/2;h(end,:)];
    dhudx=(h_ext_x(:,2:end).*u(:,2:end)-h_ext_x(:,1:end-1).*u(:,1:end-1))/dx;
    dhvdy=(h_ext_y(2:end,:).*v(2:end,:)-h_ext_y(1:end-1,:).*v(1:end-1,:))/dx;
    dhdt=-(dhudx+dhvdy);
    etta=etta+dhdt*dt;

   
    % update the visual:
    if mod(t_step,time_gap)==0
        if write_vid
            writeVideo(vid,getframe(gcf))
        end
        pause(10^-6)
        switch demo
            case {'reflect','per','shoal','superpos'}
                surface.ZData=etta(1,:);
            otherwise
                surface.ZData=etta;
        end
    end


end

if write_vid
    close(vid)
end



% refraction check, the points were determined us ginput() and copied on to
% the script:

% x_refract=[
%     2.0889
%     2.4933
%    41.7116
%    43.7332
%    72.3046
%    77.1563];
%    
% y_refract =[
%    88.3423
%    76.8868
%    87.2642
%    75.6739
%    89.4205
%    77.0216];
% 
% for i=1:3
%     x1=x_refract(i*2-1);
%     x2=x_refract(i*2);
%     y1=y_refract(i*2-1);
%     y2=y_refract(i*2);
%     plot3([x1,x2],[y1,y2],[1;1],'k','LineWidth',5)
%     angle(i)=atan((y2-y1)/(x2-x1));
%     h(i)=refract_prof((x1+x2)/2,(y1+y2)/2,M,N);
%     c(i)=sqrt(g*h(i));
% end
% figure
% scatter(c,sin(-pi/4-angle),'filled')
% hold on
% ylabel('sin(\theta)')
% xlabel('Wave Speed (m/s)')
% title('sin(\theta) vs Wave Speed During Refraction')
% xlim([0,7.5])
% ylim([0,.7])
% coeffs=polyfit(c,sin(-pi/4-angle),1);
% m=coeffs(1);
% b=coeffs(2);
% plot([0,7.5],m*[0,7.5]+b)
