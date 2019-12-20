%% Title: Class Project(Problem1)
 % Assumptions: Prismatic channel
 %              fix lateral flow
 %              Rough bed (By only adding Tolh to initial condition)
 %              Trapezoidal cross-section
 
%% Setting
 clear
 clc
 tic
 
 S0=0.0005;                          % Bed Slope
 Tl_in=[200,150,120,100];            % Lateral Flow Time [min]
 HnL=[0.33,0.44,0.55,0.66];       % Lateral Flow Rate [mm/min]
 Sideslope=1;                        % Slope of the inclined walls of the channel
 n=0.02;                             % Maning Coefficient
 L=400;                              % Channel Length [m]
 Bottom=1;                           % Bottom width
 
 Totaltime=300;                      % Total Simulation Time 
 Tolh=1d-7;                          % Water depth threshold below which the control volume is assumed to be dry
 
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



%% Governing Equation (Processing)

Q1=zeros(4,nxc);
Q2=zeros(4,nxc);
for j=1:4
    Vn1=zeros(1,nxc);  % Velocity vector at the current time step (located at the cell centers) [m/s]
    Vn=zeros(1,nxc);   % Velocity vector at previous time step (located at the cell centers) [m/s]

    Hn1=zeros(1,nxc);  % Pressure head vector at the current time step (located at the cell centers) [m]
    Hn=zeros(1,nxc);   % Pressure head vector at previous time step (located at the cell centers) [m]

    Vn3=zeros(1,nxc);  % Velocity vector at the current time step (located at the cell centers) [m/s]
    Hn3=zeros(1,nxc);  % Pressure head vector at the current time step (located at the cell centers) [m]
    Dstr=zeros(2,nxc-1);
    Vstr=zeros(2,nxc-1);
    Dstr2=zeros(2,nxc-1);
    Vstr2=zeros(2,nxc-1);
    Sf=zeros(2,nxc-1);
    Sf3=zeros(2,nxc-1);
    Tl=Tl_in(j);
    dt=1/60;                            % Time Step
    ElapsedTime=0.0;                      % Elapsed Time

    while ElapsedTime < Totaltime
        ElapsedTime=ElapsedTime+dt; 

        % Left Boundary
        Hn1(1)=0;
        Vn1(1)= (1/n)*(Hn1(1)^(2/3))*sqrt(S0);

        Vn2(:)=Vn1(:);
        Hn2(:)=Hn1(:);

       if ElapsedTime<Tl
            for i=2:nxc-1
                Dstr(i)=0.5*(Hn(i+1)+Hn(i-1));
                Vstr(i)=0.5*(Vn(i+1)+Vn(i-1));
                Hn1(i)=0.5*(Hn(i+1)+Hn(i-1))-0.5*(dt/dx)*Dstr(i)*(Vn(i+1)-Vn(i-1))...
                      -0.5*(dt/dx)*Vstr(i)*(Hn(i+1)-Hn(i-1))+(HnL(j)/60000);                %Continuity Eq.
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

            Vn(:)=Vn1(:);
            Hn(:)=Hn1(:);

        elseif ElapsedTime>Tl
            for i=2:nxc-1
                Dstr2(i)=0.5*(Hn2(i+1)+Hn2(i-1));
                Vstr2(i)=0.5*(Vn2(i+1)+Vn2(i-1));
                Hn3(i)=0.5*(Hn2(i+1)+Hn2(i-1))-0.5*(dt/dx)*Dstr2(i)*(Vn2(i+1)-Vn2(i-1))...
                      -0.5*(dt/dx)*Vstr2(i)*(Hn2(i+1)-Hn2(i-1));                %Continuity Eq.


                Sf3(i) = ((n^2)*Vn2(i)*abs(Vn2(i))/(Hn2(i)/(1+2*Hn2(i))^(2/3)));
                Vn3(i)=0.5*(Vn2(i+1)+Vn2(i-1))-0.5*(dt/dx)*9.81*(Hn2(i+1)-Hn2(i-1))...
                      -0.5*(dt/dx)*Vstr2(i)*(Vn2(i+1)-Vn2(i-1))+9.81*dt*(S0-Sf3(i)); %Momentum Eq.

            end

            % Right Boundary

            Hn3(nxc)=Hn3(nxc-1);
            Vn3(nxc)=Vn3(nxc-1);

            Vn2(:)=Vn3(:);
            Hn2(:)=Hn3(:);

        end              

            if ElapsedTime<Totaltime
                if (ElapsedTime+dt>Totaltime)
                    dt=Totaltime-ElapsedTime;  
                end
            end   

        disp(['Elapsed Time: ' num2str(ElapsedTime)])

    end
    Q1(j,1:nxc)=transpose(Hn1(:).*Vn1(:));
    Q2(j,1:nxc)=flip(Hn3(:).*Vn3(:));
end

time1=zeros(1,Totaltime);
dt1=0.50125;
for i=1:nxc-1
    time1(i+1)=time1(i)+dt1;
end
if time1(nxc-1)+dt1<200
   time1(nxc) = 200;
end

time2=zeros(1,Totaltime);
time2(1)=200;
for i=1:nxc-1
    time2(i+1)=time2(i)+0.2;
end
%% Results (Post Processing) 
% Final results
plot(time1,Q1(1,:),'b',time2,Q2(1,:),'--b',time1,Q1(2,:),'y',time2,Q2(2,:),'--y'...
    ,time1,Q1(3,:),'c',time2,Q2(3,:),'--c',time1,Q1(4,:),'g',time2,Q2(4,:),'--g')
legend('rain=0.33','dep=con#1','rain=0.44','dep=con#2','rain=0.55','dep=con#3','rain=0.66','dep=con#4')
toc