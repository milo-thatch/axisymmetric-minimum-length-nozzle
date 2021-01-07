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

function [results]=unit_process1(point1,point2)
    global th gamma %options
    
    theta1=point1(1);
    theta13=point1(1);
    theta2=point2(1);
    theta23=point2(1);
    mach1=point1(2);
    mach13=point1(2);
    mach2=point2(2);
    mach23=point2(2);
    r1=point1(4);
    r13=point1(4);
    x1=point1(3);
    %x13=point1(3);
    r2=point2(4);
    r23=point2(4);
    x2=point2(3);
    %x23=point2(3);

    err=1e3;
    mach3_i=1;
    while(abs(err)>th)
        A13=sqrt(mach13^2-1)/mach13/(1+(gamma-1)/2*mach13^2);
        A23=sqrt(mach23^2-1)/mach23/(1+(gamma-1)/2*mach23^2);
        B13=-tan(theta13)/(sqrt(mach13^2-1)+tan(theta13))/r13;
        B23=-tan(theta23)/(sqrt(mach23^2-1)*tan(theta23)+1)/r23;
        D1=A13*mach1+theta1+B13*x1;
        D2=A23*mach2-theta2+B23*r2;
        C13=tan(theta13-asin(1/mach13));
        C23=tan(theta23+asin(1/mach23));
        E1=r1-C13*x1;
        E2=r2-C23*x2;
        x3=(E1-E2)/(C23-C13);
        r3=E2+C23*x3;
        mach3=(D1+D2-B23*r3-B13*x3)/(A13+A23);
        theta3=D1-A13*mach3-B13*x3;
    
        %compute averages
        mach13=0.5*(mach1+mach3);
        mach23=0.5*(mach2+mach3);
        theta13=0.5*(theta1+theta3);
        theta23=0.5*(theta2+theta3);
        r13=0.5*(r1+r3);
    
        %compute error
        err=mach3/mach3_i-1;
        mach3_i=mach3;
    end
    results(1)=theta3;
    results(2)=mach3;
    results(3)=x3;
    results(4)=r3;
end