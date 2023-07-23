m%Project 2
%Russell Keith
%3/20/23
clear;
close all;

%number of nodes
n=150;

%interaction range
k=1.2;
d=15;
r=k*d;
dprime=8;
rprime=k*dprime;

%dimensions
m=2;

%constants
%for sigma norm and dbeta
epsilon=0.1;
ralpha=sigmanorm(r,epsilon);
dbeta=sigmanorm(dprime,epsilon);
Delta_t=0.009;

%for algorithm 2 target position
%case 1
%qgamma(1,:)=[250,25];
%pgamma(1,:)=[0, 0];

%case 2
qgamma(1,:)=[40,25];
pgamma(1,:)=[0, 0];

c1a=30;
c2a=2*sqrt(c1a);
c1mt=1.1;
c2mt=2*sqrt(c1mt);
c1b=1500;
c2b=2*sqrt(c1a);

%assign a random x y location
for i = 1:n
    x = randi([0 70]);
    y = randi([0 70]);
    %store intial position and velocity
    Qi(i, :)=[x y];
    Pi(i, :)=[0 0];
end

%make scatter plot for starting nodes
scatter(Qi(:, 1),Qi(:, 2),10,"blue","filled",">");

%add radius circles to graph
for i=1:n
    %viscircles([Qi(i,1) Qi(i,2)],r);
end

%obstacles
yk=[100 25; 150 80; 200 10];
Rk=[15 25 30];

%more contants
t1 = 1;
x=0;
K=0;
z=0;
iterations=5000;

while (iterations <= 5000)
    edgecount=0;
    energy=0;

%=======================Show current state=======================%
    %connect neighbors by drawing a line
    %plot parameters
    f = figure(2);
    movegui(f,[200,200]);

    %disp points
    scatter(Qi(:, 1),Qi(:, 2),20,"blue","filled",">");
    hold on
    %disp target
    plot(qgamma(:,1),qgamma(:,2),'Marker','x','Color','r','LineWidth',2,LineStyle='-');

    for obstacle=1:width(Rk)
        %need to download Image Proccessing Toolbox 
        viscircles(yk(obstacle,:),Rk(1,obstacle));
    end
    

%===================Calculate formula parameters===============%
    %reset adjacency matrix
    A=zeros(n,n);

    %the good stuff.....
    for node=1:n
        %get the current state
        qi=Qi(node,:);
        pi=Pi(node,:);

        %get neighbor array of node
        neighborA= findNeighbors(n,Qi,qi,r);

        %for divation energy 
        edgecount=edgecount+width(neighborA);
        %calculate current energy in system
        energy = energy + DivationEnergy(neighborA,d,qi,Qi);

        %reset sums
        SumFig=0;
        SumFid=0;

        %iterate over neighbors
        for j=1:numel(neighborA)

            %get neighbor index
            neighborIndex=neighborA(1,j);
            
            %check if node has no neighbors
            if neighborIndex ~= 0
                
                %show current state draw line connecting neighbor to node
                neighborIndex=neighborA(1,j);
                nodex=[Qi(node,1) Qi(neighborIndex,1)];
                nodey=[Qi(node,2) Qi(neighborIndex,2)];
                line(nodex, nodey);

                %update adjacency matrix
                A(node,neighborIndex)=1;
                A(neighborIndex,node)=1;

                %get position and velocity of neighbor
                qj=Qi(neighborIndex,:);
                pj=Pi(neighborIndex,:);
                 
                %calculate variables alpha
                aij=BumpFunction(sigmanorm((qj-qi),epsilon)/ralpha);
                dalpha=sigmanorm(d,epsilon);
                phialpha=ActionFunction(sigmanorm((qj-qi),epsilon),ralpha,dalpha);
                nij=(qj-qi)/(sqrt(1+epsilon*(norm(qj-qi))^2));
                %total sums
                SumFig=SumFig + phialpha*nij;
                SumFid=SumFid + aij*(pj - pi);

            end
            
        end

        %find obstacles shadows for qi
        NiBeta = findObstacle(width(Rk),yk,Rk,rprime,qi);
        BetaSumLeft=0;
        BetaSumRight=0;

        if NiBeta~=0
            for obindex=1:width(NiBeta)
                obstacle=NiBeta(1,obindex);
                obRk=Rk(1,obstacle);
                obyk=yk(obstacle,:);
                qik=BetaPosition(obyk',qi',obRk);
                pik=BetaVelocity(obyk',qi',pi',obRk);
                nik=(qik-qi')/(sqrt(1+epsilon*(norm(qik-qi'))^2));
                bik=BetaActionFunction(qik,qi,dbeta,epsilon);
                phiBeta=RepulsiveActionFunction(sigmanorm((qik-qi'),epsilon),dbeta);
                BetaSumLeft=BetaSumLeft+phiBeta*nik;
                BetaSumRight=BetaSumRight+bik*(pik-pi');
            end
        end

        %for algorithm 3
        Uigamma(1,:)=-c1mt*(Qi(node,:)-qgamma(t1,:))-c2mt*(Pi(node,:)-pgamma(t1,:));

        Uialpha=c1a*SumFig+c2a*SumFid;

        UiBeta=c1b*BetaSumLeft'+c2b*BetaSumRight';

        %get acceleration of node
        Ui=Uialpha+Uigamma+UiBeta;

        %calculating velocity from acceleration for each node
        Pi(node,:)=Pi(node,:)+Ui*Delta_t;
        %calculating position from velocity for each node
        Qi(node,:)=Qi(node,:)+Pi(node,:)*Delta_t;
    end

    hold off
    
    %keep track of velocity and trajectory
    for m=1:n
        PlotPi(m,t1)=norm(Pi(m,1),Pi(m,2));
        t1Array(1,t1)=t1;
        plotQix(m,t1)=[Qi(m,1)];
        plotQiy(m,t1)=[Qi(m,2)];
    end

%=========================Updating tagets====================%
    %update counter
    t1 = t1 + 1;

    %center of mass
    COM(t1,:)=mean(Qi);
    
    %case 1
    %for static case
    %qgamma(t1,:) = [qgamma(1,1) qgamma(1,2)];
    %pgamma(t1,:) = [0 0];
    
    %case 2
    %for moving target with sine wave tragectory
    
    if qgamma(t1-1,1) < 250
        qgamma(t1,:) = [40+50*(t1*Delta_t) (25-50*sin(t1*Delta_t))];
        pgamma(t1,:) = (qgamma(t1,:) - qgamma(t1-1,:))/Delta_t;
    else
        qgamma(t1,:) = [250 50];
        pgamma(t1,:) = [0 0];
    end
    
%=================Plotting, Connectivity, Deviation, and Energy========% 
   
    %plotting connectivity
    if mod(t1,5) == 0 || t1 == 1
        x=x+1;
        hold on 
        con = figure(20);
        movegui(con,[1000,-400]);
        Connectivity(x,:) = [t1, ((1/n) * rank(A'))];
        p = plot(Connectivity(:,1),Connectivity(:,2));
        p.Color = '#000000';
        p.LineStyle = '-';
        %p.Marker='o';
        %axis([0 t1 0 1])
        drawnow;
        
        hold off
    end
    
    %plotting Deviation energy
    if mod(t1,5) == 0 || t1 == 1
        K=K+1;
        hold on 
        Eng = figure(21);
        movegui(Eng,[1200,-600]);
        energy=energy/((edgecount/2)+1);
        disp(energy);
        Energy(K,:) = [t1, energy];
        p = plot(Energy(:,1),Energy(:,2));
        p.Color = '#000000';
        p.LineStyle = '-';
       
        drawnow;
        hold off
    end
    
    %plotting COM
    if mod(t1,20) == 0 || t1 == 1
        %disp(datetime("now"));
        com = figure(25);
        movegui(com,[1000,200]);
        plot(COM(:,1),COM(:,2),'Color','b','LineWidth',2,'LineStyle','-');
        hold on 
        plot(qgamma(:,1),qgamma(:,2),'Marker','x','Color','r','LineWidth',2,LineStyle='-');

        %for case 1 and 2 case
        axis([0 400 -100 100])
       
        drawnow;
        hold off
    end

    %plotting velocity of every node
    if mod(t1,20) == 0 || t1 == 1 
        vel = figure(23);
        movegui(vel,[1000,-200]);
        plot(t1Array,PlotPi(1,:),'Color',rand(1,3));
        hold on 
        for i=2:n
            plot(t1Array,PlotPi(i,:),'Color',rand(1,3));
        end
        %for sine wave case
        axis([50*(t1*Delta_t)-200 50*(t1*Delta_t)+200 0 400])
        hold off
    end

    %plotting trajectory of every node
    if mod(t1,20) == 0 || t1 == 1 
        trajectory = figure(24);
        movegui(trajectory,[1000,200]);
        plot(plotQix(1,:),plotQiy(1,:));
        hold on 
        for i=2:n
            plot(plotQix(i,:),plotQiy(i,:),'Color',rand(1,3));
        end
        %for case 1 and 2
        axis([0 400 -200 200])
        hold off
    end
   
end

function sn = sigmanorm(z,e)
    sn=(1/e)*(sqrt(1+(e*((norm(z))^2)))-1);
end

function euclidean = euclideannorm(i, j)
    euclidean = sqrt((j(1,1) - i(1,1))^2 + (j(1,2) - i(1,2))^2);
end

%function for finding neighbors and updating neighbor array for each node
function Ni = findNeighbors(n,Qi,qi,r)
    %returns this if it has no neighbors
    Ni(1,1)=0;
    %find nieghbors for each node
    numN=1;
    for j=1:n
        qj=Qi(j,:);
        if (qi~=qj)
            %determine if node is within range
            if norm(qj-qi) < r
                %just want the index of neighbor
                Ni(1,numN)=j;
                numN=numN+1;
            end
        end
    end
end

function NiBeta = findObstacle(k,yk,Rk,rprime,qi)
    %yk and Rk are arrays
    %return if node is not near a obstacle
    NiBeta(1,1)=0;
    numK=0;
    for obstacle=1:k
        %get position of beta agent for every obstacle
        qik=BetaPosition(yk(obstacle,:)',qi',Rk(1,obstacle));
        %determine if shadow is within range
        if norm(qik - qi') < rprime
            numK=numK + 1;
            NiBeta(1,numK)=obstacle;
        end
    end
end

%Bump function
function bf = BumpFunction(z)
    h=0.2;
    if z >= 0 && z < h
        bf = 1;
    elseif z >= h && z <= 1
        bf = (1/2)*(1+cos(pi*((z-h)/(1-h))));
    else
        bf = 0;
    end
end

%Pair wise potential function with action Function
function phiAlpha = ActionFunction(z,ralpha,dalpha)
    phiAlpha = BumpFunction(z/ralpha)*pairwise(z-dalpha);
end

function bik = BetaActionFunction(qik,qi,dbeta,epsilon)
    bik = BumpFunction(sigmanorm(qik-qi,epsilon)/dbeta);
end

function phibeta = RepulsiveActionFunction(z,dbeta)
    phibeta=BumpFunction(z/dbeta)*(sigma1(z-dbeta)-1);
end

function phi = pairwise(z)
    %constants used in simulations    
    a=5;
    b=5;
    c=abs(a-b)/(sqrt(4*a*b));
    phi=(1/2)*((a+b)*sigma1(z+c)+(a-b));
end

function sig1 = sigma1(z)
    sig1 = z/(sqrt(1+(z^2)));
end

function Eq = DivationEnergy(Ni,d,qi,Qi)
    edgecount = 0;    
    phi = 0;
    %list of neighbors is Ni
    %check to see if node has neighbors
    if Ni(1,1)~=0
        for j=1:numel(Ni)
            %keep track of edges
            edgecount=edgecount+1;
            %get qj
            neighborIndex=Ni(1,j);
            qj=Qi(neighborIndex,:);
            %calculate phi
            phi = phi + ((norm(qj-qi)-d)^2);
        end
    end
    %calculate partial energy still need to be divided by edges after
    %summing all of the engery for each node
    Eq=phi;
end

function qik = BetaPosition(yk,qi,Rk)
    %this is for a single obstacle and node 
    %Rk and yk are not arrays
    Mu=Rk/(norm(qi-yk));
    qik=(Mu*qi)+((1-Mu)*yk);
end

function pik = BetaVelocity(yk,qi,pi,Rk)
    %this is for a single obstacle and node 
    %Rk and yk are not arrays
    Mu=Rk/norm(qi-yk);
    ak=(qi-yk)/norm(qi-yk);
    I=[1 0; 0 1];
    P=I-ak*ak';
    pik=Mu*P*pi;
end
