%% Title: Class Project(Problem1)
 % Assumptions: Prismatic channel
 %              fix lateral flow
 %              Rough bed (By only adding Tolh to initial condition)
 %              Trapezoidal cross-section
 
%% Setting
 clear
 clc
 tic
 
 S0=0.0005;         % Bed Slope
 Tl=200;            % Lateral Flow Time [min]
 Hn1(1)=0.33;       % Lateral Flow Rate [mm/min]
 Sideslope=1;       % Slope of the inclined walls of the channel
 n=0.02;            % Maning Coefficient
 L=400;             % Channel Length [m]
 Bottom=1;               % Bottom width
 
 Totaltime=300;     % Total Simulation Time 
 dt=1;              % Time Step
 ElapsedTime=0;     % Elapsed Time
 Tolh=1d-7;         % Water depth threshold below which the control volume is assumed to be dry
 
%% Mesh Generation
 nxc=400;           % Number of Cell Center
 nxv=nxc+1;         % Number off vertex
 dx=L/nxc;          % length step
 
 xv(1:nxv-1)=0;
 for i=1:nxv-1
     xv(i+1)=xv(i)+dx;
 end
 
 xc(1:nxc)=0;
 for i=1:nxc 
   xc(i)=0.5*(xv(i+1)+xv(i));
 end
 
 elec(1:nxc)=0;
 elev(:)=5-S0*xv(:);
 for i=1:nxc 
   elec(i)=0.5*(elev(i+1)+elev(i));
 end

%% Geometry Calculation (Multi-dimensional Codes)
% Empty

%% Initial Condition

Vn1(1:nxc)=0;  % Velocity vector at the current time step (located at the cell centers) [m/s]
Vn(1:nxc)=0;   % Velocity vector at previous time step (located at the cell centers) [m/s]

Hn1(1:nxc)=0;  % Pressure head vector at the current time step (located at the cell centers) [m]
Hn(1:nxc)=0;   % Pressure head vector at previous time step (located at the cell centers) [m]

%% Governing Equation (Processing)

Dstr(2:nxc-1)=0;
Vstr(2:nxc-1)=0;
Sf(2:nxc-1)=0;
while ElapsedTime < Totaltime
    ElapsedTime=ElapsedTime+dt; 

    % Left Boundary
    Vn1(1)= (1/n)*(Hn1(1)^(2/3))*sqrt(S0);
        
        if ElapsedTime>200
            Hn(1)=0;
        end
        
    for i=2:nxc-1
        Dstr(i)=0.5*(Hn(i+1)+Hn(i-1));
        Vstr(i)=0.5*(Vn(i+1)+Vn(i-1));
        Hn1(i)=0.5*(Hn(i+1)+Hn(i-1))-0.5*(dt/dx)*Dstr(i)*(Vn(i+1)-Vn(i-1))...
              -0.5*(dt/dx)*Vstr(i)*(Hn(i+1)-Hn(i-1))+0.00033;           %Continuity Eq.
        if Hn1(i)<Tolh
            Hn1(i)=0;
            Vn1(i)=0;
        else
            Sf(i) = ((n^2)*Vn1(i)*abs(Vn1(i))/(Hn1(i)/(1+2*Hn1(i))^(2/3)) );
            Vn1(i)=0.5*(Vn(i+1)+Vn(i-1))-0.5*(dt/dx)*9.81*(Hn(i+1)-Hn(i-1))...
                  -0.5*(dt/dx)*Vstr(i)*(Vn(i+1)-Vn(i-1))+9.81*dt*(S0-Sf(i)); %Momentum Eq.
        end
      
    end
    % Right Boundary
        Hn1(nxc)=Hn1(nxc-1);
        Vn1(nxc)=Vn1(nxc-1);
        
        if ElapsedTime<Totaltime
            if (ElapsedTime+dt>Totaltime)
                dt=Totaltime-ElapsedTime;  
            end
        end
    
    
        Vn(:)=Vn1(:);
        Hn(:)=Hn1(:);

    disp(['Elapsed Time: ' num2str(ElapsedTime)])

end

disp(['Elapsed Time: ' num2str(ElapsedTime)])

Q=Hn1(i+1)*Vn;

time=zeros(1,Totaltime);
for i=1:nxc-1
    time(i+1)=time(i)+dt;
end

%% Results (Post Processing) 
% Final results
plot(time,Q)
% plot(xv,elev,'k',xc,Hn1+elec,'b')
% legend('Bed Elevation(m)','Water Head(m)')
toc