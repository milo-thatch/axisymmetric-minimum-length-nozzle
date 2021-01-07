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

function [results]=unit_process3(point1)
    global th gamma 
    
    theta1=point1(1);
    theta13=point1(1);
    mach1=point1(2);
    mach13=point1(2);
    r1=point1(4);
    r13=point1(4);
    x1=point1(3);

    r3=0;
    err=1e3;
    mach3_i=1;
    while(abs(err)>th)
        A13=sqrt(mach13^2-1)/mach13/(1+(gamma-1)/2*mach13^2);
        B13=-tan(theta13)/(sqrt(mach13^2-1)+tan(theta13))/r13;
        C13=tan(theta13-asin(1/mach13));
        x3=(r3-r1)/C13+x1;
        mach3=mach1+(theta1-B13*(x3-x1))/A13;
    
        %compute averages
        mach13=0.5*(mach1+mach3);
        r13=0.5*(r1+r3);
    
        %compute error
        err=mach3/mach3_i-1;
        mach3_i=mach3;
    end
    results(1)=0;
    results(2)=mach3;
    results(3)=x3;
    results(4)=r3;
end