function [XY,kernel,transition]=AxSymMLN(mach_exit,Gamma,N_lin,N_comp,index_comp,TOL,AR)
%Ludovico Fossà 01/2021

%algorithm based on 
%--Anderson, J.D. "Modern Compressible Flow", McGrawHill Education, Third
%Edition, 2003
%--Argrow, B.M. and Emanuel, G., "Comparison of Minimum Length Nozzles", 
%Journal of Fluids Engineering, 1988
%--Foelsch, K. "The analytical design of an axially symmetric laval nozzle 
%for aparallel and uniform jet", Journal of the Aeronautical Sciences, 
%1949.
%--Ying-Nien Yu "A summary of design techniques for axisymmetric hypersonic
%wind tunnels", Technical report, NATO-Science&TechnologyOrganization,1958

global th gamma options NEXT
options = optimset('Display','off');

%% INPUT
% mach_exit=2.68;
% gamma=1.4;
% N_lin=20; %linear kernel
% N_comp=5; %compressed kernel - set 1 for uncompressed kernel
% index_comp=5; %compression exponent
% plot_text=false;
% th=1e-7; %tolerance
% AR=1; %aspect ratio for the transition region
gamma=Gamma;
th=TOL;

%% OUTPUT VALUES
theta_wmax=prandtl_meyer(mach_exit,gamma)/4; %AXSYM flows
delta_theta=theta_wmax/N_lin; %linear kernel
steps=N_lin+N_comp-1;

%number of points - KERNEL REGION
n_kernel=0;
for i=0:steps
    n_kernel=n_kernel+steps-i;
end
kernel=zeros(n_kernel,4);
sel=@(x) steps+1+x.*(steps-0.5*(x+1)); %new series
figure(1)
grid
daspect([1 1 0.2778])
hold on

%% KERNEL REGION
% Compute points on the first characteristic line (origin)
% Compression region
% Compute point 1 - 3RD UNIT PROCESS
theta1=(1/N_comp).^index_comp*delta_theta; %theta
mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
r1=1; %ORIGIN
x1=0; %ORIGIN
kernel(1,:)=unit_process3([theta1,mach1,x1,r1]);
line([x1 kernel(1,3)],[r1 kernel(1,4)])
% Compute point 2 - 2ND UNIT PROCESS
if(N_comp>=2)
    theta1=(2/N_comp).^index_comp*delta_theta;
    mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
    r1=1; %ORIGIN
    x1=0; %ORIGIN
    kernel(2,:)=unit_process2([theta1,mach1,x1,r1],kernel(1,:));
    line([x1 kernel(2,3)],[r1 kernel(2,4)])
end
% Compute internal kernel points - 1ST UNIT PROCESS
for i=3:N_comp
    theta1=(i/N_comp).^index_comp*delta_theta;
    mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
    r1=1; %ORIGIN
    x1=0; %ORIGIN
    kernel(i,:)=unit_process1([theta1,mach1,x1,r1],kernel(i-1,:));
    line([x1 kernel(i,3)],[r1 kernel(i,4)])
end
% Linear region
if(N_comp>=2)
    for i=2:N_lin
        theta1=i*delta_theta;
        mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
        r1=1; %ORIGIN
        x1=0; %ORIGIN
        kernel(N_comp+i-1,:)=unit_process1([theta1,mach1,x1,r1],kernel(N_comp+i-2,:));
        line([x1 kernel(N_comp+i-1,3)],[r1 kernel(N_comp+i-1,4)])
    end
else
    i=2; % run unit process 2
    theta1=i*delta_theta;
    mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
    r1=1; %ORIGIN
    x1=0; %ORIGIN
    kernel(N_comp+i-1,:)=unit_process2([theta1,mach1,x1,r1],kernel(N_comp+i-2,:));
    line([x1 kernel(N_comp+i-1,3)],[r1 kernel(N_comp+i-1,4)])
    for i=3:N_lin
        theta1=i*delta_theta;
        mach1=fsolve(@(x) prandtl_meyer(x,gamma)-theta1,1.1,options); %mach
        r1=1; %ORIGIN
        x1=0; %ORIGIN
        kernel(N_comp+i-1,:)=unit_process1([theta1,mach1,x1,r1],kernel(N_comp+i-2,:));
        line([x1 kernel(N_comp+i-1,3)],[r1 kernel(N_comp+i-1,4)])
    end
end
% Compute the remaining characteristics
% Compression region
for j=0:N_comp-1
    %Centerline - 3rd unit process
    kernel(sel(j),:)=unit_process3(kernel(sel(j-1)+1,:));
    %line([kernel(sel(j-1),3) kernel(sel(j-1)+1,3)],[kernel(sel(j-1),4) kernel(sel(j-1)+1,4)])
    line([kernel(sel(j-1),3) kernel(sel(j-1)+1,3)],[kernel(sel(j-1),4) kernel(sel(j-1)+1,4)])
    line([kernel(sel(j-1)+1,3) kernel(sel(j),3)],[kernel(sel(j-1)+1,4) kernel(sel(j),4)])
    %Adjacent - 2nd unit process
    for i=1:N_comp-j-1
        kernel(sel(j)+1,:)=unit_process2(kernel(sel(j-1)+2,:),kernel(sel(j),:));
        line([kernel(sel(j-1)+2,3) kernel(sel(j)+1,3)],[kernel(sel(j-1)+2,4) kernel(sel(j)+1,4)])
        line([kernel(sel(j-1)+1,3) kernel(sel(j-1)+2,3)],[kernel(sel(j-1)+1,4) kernel(sel(j-1)+2,4)])
    end
    %Internal kernel points - 1st unit process
    for i=2:N_comp-j-1
        kernel(sel(j)+i,:)=unit_process1(kernel(sel(j-1)+i+1,:),kernel(sel(j)+i-1,:));
        line([kernel(sel(j-1)+i,3) kernel(sel(j-1)+i+1,3)],[kernel(sel(j-1)+i,4) kernel(sel(j-1)+i+1,4)])
        line([kernel(sel(j-1)+i+1,3) kernel(sel(j)+i,3)],[kernel(sel(j-1)+1+i,4) kernel(sel(j)+i,4)])
    end
end
% Linear region
for j=0:steps-2
    if(j<=N_comp)
        ind_ext=sel(j)+N_comp-j;
    else
        ind_ext=sel(j);
    end
    for i=ind_ext:sel(j+1)-1
        k=i+1-sel(j);
        switch k
            case 1
                %UNIT PROCESS 3
                kernel(sel(j),:)=unit_process3(kernel(sel(j-1)+1,:));
                line([kernel(sel(j-1),3) kernel(sel(j-1)+1,3)],[kernel(sel(j-1),4) kernel(sel(j-1)+1,4)])
                line([kernel(sel(j),3) kernel(sel(j-1)+1,3)],[kernel(sel(j),4) kernel(sel(j-1)+1,4)])
            case 2
                %UNIT PROCESS 2
                kernel(sel(j)+1,:)=unit_process2(kernel(sel(j-1)+2,:),kernel(sel(j),:));
                line([kernel(sel(j-1)+1,3) kernel(sel(j-1)+2,3)],[kernel(sel(j-1)+1,4) kernel(sel(j-1)+2,4)])
                line([kernel(sel(j-1)+2,3) kernel(sel(j)+1,3)],[kernel(sel(j-1)+2,4) kernel(sel(j)+1,4)])
            otherwise
                %UNIT PROCESS 1
                ind_inn=sel(j-1)+(i-sel(j)+1);
                kernel(i,:)=unit_process1(kernel(ind_inn,:),kernel(i-1,:));
                line([kernel(ind_inn-1,3) kernel(ind_inn,3)],[kernel(ind_inn-1,4) kernel(ind_inn,4)])
                line([kernel(ind_inn,3) kernel(i,3)],[kernel(ind_inn,4) kernel(i,4)])
        end
    end
end

%% TRANSITION REGION
transition=zeros(2*steps*(steps+1),4);
XY=zeros(2*steps+1,2);
% Compute mean distance on the last right-running characteristic
S_ave=mean(sqrt(diff(kernel(sel(0:steps-1)-1,3)).^2+diff(kernel(sel(0:steps-1)-1,4)).^2));
S_plus=AR*S_ave;
mu_exit=asin(1/kernel(end,2));
%x_exit=sqrt(area_mach_nozzle(mach_exit,gamma))/tan(mu_exit)+kernel(end,3);

% Compute first transition point
transition(1,2)=kernel(end,2);
transition(1,3)=S_plus*cos(mu_exit)+kernel(end,3);
transition(1,4)=S_plus*sin(mu_exit)+kernel(end,4);
line([kernel(end,3) transition(1,3)],[kernel(end,4) transition(1,4)],'Color','black')
% Compute first transition line
for j=steps-2:-1:0
    transition(steps-j,:)=unit_process1(transition(steps-j-1,:),kernel(sel(j)-1,:)); 
    line([transition(steps-j,3) transition(steps-j-1,3)],[transition(steps-j,4) transition(steps-j-1,4)],'Color','black')
    line([transition(steps-j,3) kernel(sel(j)-1,3)],[transition(steps-j,4) kernel(sel(j)-1,4)],'Color','black')
end
transition(steps-j+1,:)=unit_process1(transition(steps-j,:),[transition(steps-j,1) transition(steps-j,2) 0 1]);
line([transition(steps-j+1,3) transition(steps-j,3)],[transition(steps-j+1,4) transition(steps-j,4)],'Color','black')
line([0 transition(steps-j+1,3)],[1 transition(steps-j+1,4)],'Color','black')

% Compute initial contour slope (theta_wmax)
m_lim=(transition(steps,4)-1)/transition(steps,3);
NEXT=false;
cc=steps; %identify column on the transition web
bb=0; %identify row on the transition web
XY(1,:)=[0 1]; %initial contour point (throat region)
k=2; %contour counter
while(NEXT==false)
    if(tan(theta_wmax)>m_lim)
        %fprintf('Intersects C- %d\n',k)
        NEXT=true;
        point_new=wall_minus(theta_wmax,XY(k-1,:),transition(bb*(steps+1)+cc,:),transition(bb*(steps+1)+cc+1,:)); %theta,r,x
        XY(k,:)=[point_new(3),point_new(2)]; %SAVE CONTOUR POINT
        line([0 XY(k,1)],[1 XY(k,2)],'Color','black','LineWidth',4)
        k=k+1;
        bb=bb+1;
    else
        %fprintf('Intersects C+ %d\n',k)
        NEXT=false;
        point_new=wall_plus(theta_wmax,XY(k-1,:),transition((bb-1)*(steps+1)+cc,:),transition(bb*(steps+1)+cc+1,:)); %theta,r,x
        XY(k,:)=[point_new(3),point_new(2)]; %SAVE CONTOUR POINT
        line([0 XY(k,1)],[1 XY(k,2)],'Color','black','LineWidth',4)
        k=k+1;
        cc=cc-1;
    end
end

i=0;
OUTLET=false;
while(OUTLET==false)
    i=i+1;
    transition(i*(steps+1)+1,2)=kernel(end,2);
    transition(i*(steps+1)+1,3)=S_plus*cos(mu_exit)+transition((i-1)*(steps+1)+1,3);
    transition(i*(steps+1)+1,4)=S_plus*sin(mu_exit)+transition((i-1)*(steps+1)+1,4);
    line([transition(i*(steps+1)+1,3) transition((i-1)*(steps+1)+1,3)],[transition(i*(steps+1)+1,4) transition((i-1)*(steps+1)+1,4)],'Color','black')
    for j=steps-2:-1:-1+(steps-cc)
        transition(i*(steps+1)+steps-j,:)=unit_process1(transition(i*(steps+1)+steps-j-1,:),transition((i-1)*(steps+1)+steps-j,:)); 
        line([transition(i*(steps+1)+steps-j,3) transition(i*(steps+1)+steps-j-1,3)],[transition(i*(steps+1)+steps-j,4) transition(i*(steps+1)+steps-j-1,4)],'Color','black')
        line([transition(i*(steps+1)+steps-j,3) transition((i-1)*(steps+1)+steps-j,3)],[transition(i*(steps+1)+steps-j,4) transition((i-1)*(steps+1)+steps-j,4)],'Color','black')
    end
    
    NEXT=false;
    while(NEXT==false && OUTLET==false)
        theta0=point_new(1);
        m_lim=(transition(bb*(steps+1)+cc,4)-XY(k-1,2))/(transition(bb*(steps+1)+cc,3)-XY(k-1,1)); %constant
        point_new=wall_contour(m_lim,theta0,XY(k-1,:),...
            transition((bb-1)*(steps+1)+cc,:),transition(bb*(steps+1)+cc,:),transition(bb*(steps+1)+cc+1,:));
        XY(k,:)=[point_new(3),point_new(2)]; %SAVE CONTOUR POINT
        line([XY(k-1,1) XY(k,1)],[XY(k-1,2) XY(k,2)],'Color','black','LineWidth',4)
        k=k+1;
        if(NEXT)
            bb=bb+1;
        else
            cc=cc-1;
            if(cc<=0) 
                OUTLET=true;
            end
        end
    end
end
title('Axisymmetric minimum-length nozzle contour','FontSize',20,'FontWeight','bold')
xlabel('x/r_{throat}','FontSize',16,'FontWeight','bold')
ylabel('y/r_{throat}','FontSize',16,'FontWeight','bold')

%% EXPORT
export=zeros(size(XY,1),3);
export(:,1)=1;
export(:,2:3)=XY(:,:);

format short
writematrix(export,'nozzle_axi.txt','Delimiter','tab')

end