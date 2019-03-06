function [tablee3] = TRUSS

%{
The MATLAB code prepared for solving any truss system with or without 
temperature changes and fabrication errors. Hence, it needs user input. 
Please input a simple truss system which has proper element lengths 
like 10 m for better visual figure and calculation.
%}

format long
%taking EA inputs
prompt = {'Enter E kN/m2:','Enter A m2:'};
title = 'E and A';
dims = [1 70];
definput = {'200000000','0.003249'};
ema = inputdlg(prompt,title,dims,definput);
ea= str2double(ema{1})*str2double(ema{2});

%taking node numbers
prompt = {'Number of nodes:'};
title = 'Nodes';
dims = [1 70];
definput = {'2'};
nn = inputdlg(prompt,title,dims,definput);
node=str2double(nn{1});
nodedata=zeros(node,9);

for i=1:node

%takes input as node coordinates and support conditions
prompt = {'Enter node x dir m:','Enter node y dir m:','is support 1 or 0','Force at x kN','Force at y kN'};
title = ['Node: ',int2str(i)];
dims = [1 70];definput = {'0','0','1','0','0'};
nxy = inputdlg(prompt,title,dims,definput);

%saving node data
nodedata(i,:)=[str2double(nxy{1}),str2double(nxy{2}),str2double(nxy{3}),0,0,str2double(nxy{4}),str2double(nxy{5}),0,0];
end

for i=1:node
    if nodedata(i,3) == 1
        %taking support fixities and saves nodedata
        prompt = {'Enter x fixity 1 or 0:','Enter y fixity 1 or 0:'};
        title = ['Support at Node ',int2str(i),':'];
        dims = [1 70];
        definput = {'1','1'};
        ss = inputdlg(prompt,title,dims,definput);
        nodedata(i,4)=str2double(ss{1});
        nodedata(i,5)=str2double(ss{2});
    end
end

%calculates dof matrix 
dofdata=zeros(node,2);
sk = 0;
for i=1:node
    for j=4:5
        if nodedata(i,j) == 0
            sk = sk+1;
            dofdata(i,j-3)=sk;
        end
    end
end

%Kglobal matrix with using dof number for structure
KG=zeros(sk,sk);

%taking member numbers
prompt = {'Number of members: '};
title = 'Members';
dims = [1 70];
definput = {'1'};
mm = inputdlg(prompt,title,dims,definput);
member=str2double(mm{1});
memberdata=zeros(member,17);

for i=1:member

%taking member nodes and calculates lenght, cos, sin of the member
prompt = {'Enter first node','Enter second node','Enter member EA(if different):'};
title = ['Member: ',int2str(i)];
dims = [1 70];
definput = {'1','2',num2str(ea)};
mmb = inputdlg(prompt,title,dims,definput);
n1=str2double(mmb{1});
n2=str2double(mmb{2});
x1=nodedata(n1,1);
x2=nodedata(n2,1);
y1=nodedata(n1,2);
y2=nodedata(n2,2);
x=x2-x1;
y=y2-y1;
l=sqrt(x*x+y*y);
c=x/l;
s=y/l;
cs=c*s;
c2=c*c;
s2=s*s;

memberdata(i,:)=[x1,y1,x2,y2,c2,s2,cs,str2double(mmb{3}),n1,n2,0,0,0,0,l,c,s];   
end


%if there is temperature change in the problem it takes data
prompt = {'Is there any Temperature Change ? (1=yes or 0=no)'};
title = 'Temperature Change';
dims = [1 70];
definput = {'0'};
dt = inputdlg(prompt,title,dims,definput);

tempfdata=zeros(member,1);
if str2double(dt) == 1
    prompt = {'Input Thermal Coefficient'};
    title = 'Thermal Coefficient';
    dims = [1 70];
    definput = {'0.0000099'};
    tc = inputdlg(prompt,title,dims,definput);
    thc = str2double(tc{1});
    
        for i=1:member
            prompt = {'Temperature Change'};
            title = ['For Member ', int2str(i),' x: (', int2str(memberdata(i,1)),',',int2str(memberdata(i,2)),') y: (',int2str(memberdata(i,3)),',',int2str(memberdata(i,4)),')'];
            dims = [1 70];
            definput = {'0'};
            dtt = inputdlg(prompt,title,dims,definput);
            tempfdata(i,1)=thc*memberdata(i,8)*str2double(dtt{1});
        end
end

%if there is fabrication error in the problem it takes data
prompt = {'Is there any Fabrication Errors ? (1=yes or 0=no)'};
title = 'Fabrication Error (m)';
dims = [1 70];
definput = {'0'};
df = inputdlg(prompt,title,dims,definput);

errorfdata=zeros(member,1);
if str2double(df) == 1
    
        for i=1:member
            prompt = {'Fabrication Error'};
            title = ['For Member ', int2str(i),' x: (', int2str(memberdata(i,1)),',',int2str(memberdata(i,2)),') y: (',int2str(memberdata(i,3)),',',int2str(memberdata(i,4)),')'];
            dims = [1 70];
            definput = {'0'};
            dft = inputdlg(prompt,title,dims,definput);
            errorfdata(i,1)=memberdata(i,8)*str2double(dft{1})/memberdata(i,15);
        end
end

%calculates Kglobal of every member
KL= zeros(4,4,member);
for i=1:member
    KL(1,1,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,5);
    KL(1,2,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(1,3,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,5);
    KL(1,4,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(2,1,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(2,2,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,6);
    KL(2,3,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(2,4,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,6);
    KL(3,1,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,5);
    KL(3,2,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(3,3,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,5);
    KL(3,4,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(4,1,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(4,2,i)=-memberdata(i,8)/memberdata(i,15)*memberdata(i,6);
    KL(4,3,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,7);
    KL(4,4,i)=memberdata(i,8)/memberdata(i,15)*memberdata(i,6);
end


%pairs the members with dofs
for i=1:member
        memberdata(i,11)=dofdata(memberdata(i,9),1);
        memberdata(i,12)=dofdata(memberdata(i,9),2);
        memberdata(i,13)=dofdata(memberdata(i,10),1);
        memberdata(i,14)=dofdata(memberdata(i,10),2);
end


%enters Kglobal, parts of Kmember's
for i=1:member
    for x=1:4
        for y=1:4
            if memberdata(i,10+x) ~= 0 && memberdata(i,10+y) ~= 0
                KG(memberdata(i,10+x),memberdata(i,10+y))=KL(x,y,i)+KG(memberdata(i,10+x),memberdata(i,10+y));
            end
        end
    end
end


for i=1:node
    nodedata(i,8)=dofdata(i,1);
    nodedata(i,9)=dofdata(i,2);
end

%creates force matrix
forcedata=zeros(sk,1);
for i=1:node
    for j=1:2
        if dofdata(i,j) ~= 0 
            if j==1
            forcedata(dofdata(i,j))=forcedata(dofdata(i,j))+nodedata(i,6);
            elseif j==2
            forcedata(dofdata(i,j))=forcedata(dofdata(i,j))+nodedata(i,7);
            end
        end
    end
end


%adds temperature and fabrication errors as FEM
for i=1:member
    if tempfdata(i,1) ~= 0 || errorfdata(i,1) ~= 0
        for j=11:14
            if memberdata(i,j) ~= 0
                if j == 11
                forcedata(memberdata(i,j),1) = forcedata(memberdata(i,j),1) + -tempfdata(i,1)*memberdata(i,16)+ -errorfdata(i,1)*memberdata(i,16);
                elseif j ==12
                forcedata(memberdata(i,j),1) = forcedata(memberdata(i,j),1) + -tempfdata(i,1)*memberdata(i,17)+ -errorfdata(i,1)*memberdata(i,17);    
                elseif j ==13
                forcedata(memberdata(i,j),1) = forcedata(memberdata(i,j),1) + tempfdata(i,1)*memberdata(i,16)+ errorfdata(i,1)*memberdata(i,16);
                elseif j== 14
                forcedata(memberdata(i,j),1) = forcedata(memberdata(i,j),1) + tempfdata(i,1)*memberdata(i,17)+ errorfdata(i,1)*memberdata(i,17);
                end
            end
        end
    end
end

%solves for u
u=linsolve(KG,forcedata);

%arranges deflection matrices for every member
KF=zeros(member,4);
KLD=zeros(4,1,member);
for i=1:member
    for j=11:14
        if memberdata(i,j) ~= 0
            KLD(j-10,1,i)= u(memberdata(i,j));
        end
    end
end

%calculates member forces using partitioning
for i=1:member
    
    mult=[-memberdata(i,16),-memberdata(i,17),memberdata(i,16),memberdata(i,17)] * KLD(:,:,i);
    KF(i,1)=i;
    KF(i,2)=memberdata(i,8)/memberdata(i,15)*mult;
    KF(i,3)=-tempfdata(i,1);
    KF(i,4)=-errorfdata(i,1);
    KF(i,5)=KF(i,2)+KF(i,3)+KF(i,4);
end

tablee = table(KF(:,1),KF(:,2),KF(:,3),KF(:,4),KF(:,5));
tablee.Properties.VariableNames = {'Member_No','Force_kN','Force_Temp_kN','Force_F_Error_kN','Total_kN'};

tableeee = table(KF(:,1),KF(:,5),memberdata(:,15));
tableeee.Properties.VariableNames = {'Member_No','Force_kN','Member_Length_m'};

KS=zeros(node,1);
for i=1:node
    KS(i,1)=i;
end

%calculates support reactions using member forces
KSF=zeros(node,2);
for i=1:member
    for j=9:10;
        if j==9
        KSF(memberdata(i,j),:)=[memberdata(i,16)*-KF(i,5),memberdata(i,17)*-KF(i,5)]+KSF(memberdata(i,j),:);
        elseif j==10
        KSF(memberdata(i,j),:)=[memberdata(i,16)*KF(i,5),memberdata(i,17)*KF(i,5)]+KSF(memberdata(i,j),:);
        end
    end
end


j=1;
for i=1:node
    if nodedata(i,3) == 0
        KSF(j,:)=[];
        KS(j,:)=[];
        j=j-1;
    end
    j=j+1;
end

tablee2 = table(KS(:,1),KSF(:,1),KSF(:,2));
tablee2.Properties.VariableNames = {'Node_No','X_dir_Force_kN','Y_dir_Force_kN'};

DS=zeros(sk,2);

cnt=1;
for i=1:node
    for j=1:2
        if dofdata(i,j) ~= 0
            DS(cnt,1)=i;
            DS(cnt,2)=dofdata(i,j);
            cnt=cnt+1;
        end
    end
end


tablee3 = table(DS(:,1),DS(:,2),u(:,1));
tablee3.Properties.VariableNames = {'Node_No','Dof_No','Displacement_m'};

tablee4 = table(DS(:,1),DS(:,2),forcedata(:,1));
tablee4.Properties.VariableNames = {'Node_No','Dof_No','Force_kN'};


%outputs
fprintf('\nMember Global Matrices\n');
disp(KL);
fprintf('\nK(structure)\n');
disp(KG);
fprintf('\nForce Matrix\n');
disp(tablee4);
fprintf('\nSupport Reactions\n');
disp(tablee2);
fprintf('\nMember Forces\n');
disp(tableeee);

hold on
%calculates plot area
gx=0;
for i=1:member
    if gx < memberdata(i,3)
        gx=memberdata(i,3);
    end
    if gx < memberdata(i,1)
        gx=memberdata(i,1);
    end
end

gy=0;
for i=1:member
    if gy < memberdata(i,4)
        gy=memberdata(i,4);
    end
    if gy < memberdata(i,2)
        gy=memberdata(i,2);
    end
end

xlim([-2,gx+2]);
ylim([-2,gy+2]);

%plots members
for i=1:member
    line([memberdata(i,1),memberdata(i,3)],[memberdata(i,2),memberdata(i,4)]);
end

%plots dofs
for i=1:node
    for j=1:2
        if dofdata(i,j) ~= 0
            if j==1
                quiver(nodedata(i,1)-0.5,nodedata(i,2)+0.5,0.5,0,'MaxHeadSize',1.5,'Color','r','LineWidth',2);
                text(nodedata(i,1)+0.02,nodedata(i,2)+0.5,int2str(dofdata(i,j)));
            else
                quiver(nodedata(i,1)-0.5,nodedata(i,2)+0.5,0,0.5,'MaxHeadSize',1.5,'Color','r','LineWidth',2);
                text(nodedata(i,1)-0.5,nodedata(i,2)+1.02,int2str(dofdata(i,j)));
            end
        end
    end
end

end