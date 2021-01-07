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

function results=wall_contour(m_lim,theta0,point0,pointA,pointB,pointC)
    global th NEXT
        
    tol=1e3;
    theta_new=theta0;
    point_new=[point0(2);point0(1)]; %YX
    while(tol>th)
        if(tan(0.5*(theta0+theta_new))>=m_lim) %important condition
            %intersects C-
            NEXT=true;
            mCB=(pointC(4)-pointB(4))/(pointC(3)-pointB(3));
            A=@(INC) [1,-mCB;1,-tan(0.5*(theta0+INC))]; %2x2 matrix
            B=@(INC) [pointC(4)-mCB*pointC(3);point0(2)-tan(0.5*(theta0+INC))*point0(1)];
            
            old_point=point_new; %YX
            old_theta=theta_new; %YX
            point_new=A(theta_new)\B(theta_new); %array gets YX
            theta_new=pointC(1)+sqrt(((point_new(2)-pointB(3))^2+(point_new(1)-pointB(4))^2)/...
                ((pointC(3)-pointB(3))^2+(pointC(4)-pointB(4))^2))*(pointB(1)-pointC(1));
            tol=max(abs(point_new(1)-old_point(1)),abs(point_new(2)-old_point(2)));
            tol=max(tol,abs(theta_new-old_theta));
        else
            %intersects C+
            NEXT=false;
            mAB=(pointB(4)-pointA(4))/(pointB(3)-pointA(3));
            A=@(INC) [1,-mAB;1,-tan(0.5*(theta0+INC))]; %2x2 matrix
            B=@(INC) [pointA(4)-mAB*pointA(3);point0(2)-tan(0.5*(theta0+INC))*point0(1)];
            
            old_point=point_new; %YX
            old_theta=theta_new; %YX
            point_new=A(theta_new)\B(theta_new); %array gets YX
            theta_new=pointA(1)+sqrt(((point_new(2)-pointA(3))^2+(point_new(1)-pointA(4))^2)/...
                ((pointB(3)-pointA(3))^2+(pointB(4)-pointA(4))^2))*(pointB(1)-pointA(1));
            tol=max(abs(point_new(1)-old_point(1)),abs(point_new(2)-old_point(2)));
            tol=max(tol,abs(theta_new-old_theta));
        end
    end
    results=[theta_new,point_new'];
end