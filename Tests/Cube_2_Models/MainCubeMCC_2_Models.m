close all;
clear;
%% Test: CUBE Elastic vs ElastoPlastic

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

% -------------------------- SET THE PHYSICS -------------------------

model = ModelType("Poromechanics_FEM");


% ----------------------- SIMULATION PARAMETERS ----------------------

fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);


% ------------------------------  MESH -------------------------------

% Create the Mesh object
topology = Mesh();

% Set the input file name
meshName = 'Mesh/CubeTetra5.msh';

% Import the mesh data into the Mesh object
topology.importGMSHmesh(meshName);

% MCCresults=load('camclay.dat');
% F=MCCresults(:,2);
% vectimeF=linspace(0,20,length(F));
% vecF=-(F-1.e5)/1.e3;
% plot(vectimeF,vecF);
vecF=[0 -50 -50.001 -50.002 -100 -100.001 -100.002 -150 -150.001 -150.002];
% vecF=[0 -10 -20 -30 -40 -50 -60 -70 -80 -90 -100 -110 -120 -130 -120 -110 -100 -90 -80 -70 -60 -50 -60 -70 -80 -90 -100 -110 -120 -130 -140 -150 -160];
vectimeF=linspace(0,simParam.tMax,length(vecF));
%% Case: Elastic


%----------------------------- MATERIALS -----------------------------

% Set the input file name
%
fileName = 'materialsList_Elastic.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);


%------------------------------ ELEMENTS -----------------------------

% Create object handling gauss point integration
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager 
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
Elastic_linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
Elastic_state = Elastic_linSyst.setState();

% Create and set the print utility
Elastic_printUtils = OutState(model,topology,'outTime.dat','folderName','Output_CubeElastic');


%------------------------ BOUNDARY CONDITIONS ------------------------
    
% UTILITY TO WRITE BC FILES 
writeBCfiles('BCs/dirPoroBottom','NodeBC','Dir',{'Poro','x','y','z'},'bottom_fix',0,0,topology,2);
writeBCfiles('BCs/neuPoroZ','SurfBC','Neu',{'Poro','z'},'Top_LoadZ',vectimeF,vecF,topology,3);


% Collect BC input file in a list
fileName = ["BCs/dirPoroBottom.dat","BCs/neuPoroZ.dat"];



% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);


%------------------------ INITIAL STATE ------------------------
% GeoParameters
nu=mat.db(1).ConstLaw.nu;
E=mat.db(1).ConstLaw.E;
depth=-5; % m
gamma=20; % kN/m^3 
M1=1;
M2=1; 
% M1=nu/(1-nu);
% M2=nu/(1-nu); 

% Cells Initial State
cellslist=dofmanager.getFieldCells("Poromechanics");
for i=1:length(cellslist)
    cID=cellslist(i);
    z=depth+elems.cellCentroid(cID,3); %Real Depth

    % Initial State
    sigmaZ=gamma*z;
    sigmaX=M1*sigmaZ;
    sigmaY=M2*sigmaZ;
    Elastic_state.iniStress(cID,1)=sigmaX;
    Elastic_state.iniStress(cID,2)=sigmaY;
    Elastic_state.iniStress(cID,3)=sigmaZ;
    p(i)=1/3*(sigmaX+sigmaY+sigmaZ);

    % Average Preconsolidation Stress pc
end

if elems.mesh.cellVTKType(1)==10 % element is tetra
    nptGauss=1;
elseif elems.mesh.cellVTKType(1)==12 % element is exa
    nptGauss=8;
end

if isa(mat.db(1).ConstLaw, 'ModifiedCamClay')
    pc=repmat(pc',1,nptGauss);
    newmap=mat.db(1);
    newmap.ConstLaw.pc=-pc;
    mat.db(1)=newmap;
elseif isa(mat.db(1).ConstLaw, 'Elastic')
%     matfileMCC=importdata('Materials\PorousMediaMCC.dat');
%     data=str2num(matfileMCC.textdata{2, 1}(1:70));
%     e=data(1);
%     k=data(6);
%     E=(1+e)/k*p;
%     newmap=mat.db(1);
%     newmap.ConstLaw.E=-E; 
end

Elastic_state.curr.stress = Elastic_state.iniStress;
Elastic_state.conv.stress = Elastic_state.iniStress;


% Print model initial state
Elastic_printUtils.printState(Elastic_linSyst,Elastic_state);

% ---------------------------- SOLUTION -------------------------------

% Create the object handling the (nonlinear) solution of the problem
solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,Elastic_printUtils,Elastic_state,Elastic_linSyst,GaussPts);
[Elastic_simState, Elastic_endState] = solver.NonLinearLoop();


% Finalize the print utility
Elastic_printUtils.finalize()
Elastic_Disp=Elastic_printUtils.results.expDispl;



%% Case Elastoplastic (MCC)
fid=fopen('HN_p_c.txt','w');
fprintf(fid,'\n');
fid=fopen('axis_sigma_z.txt','w');
fprintf(fid,'\n');
fid=fopen('axis_eps_z.txt','w');
fprintf(fid,'\n');
fid=fopen('HN_p.txt','w');
fprintf(fid,'\n');
fid=fopen('HN_q.txt','w');
fprintf(fid,'\n');
fid=fopen('HN_sigma_z.txt','w');
fprintf(fid,'\n');
%----------------------------- MATERIALS -----------------------------

% Set the input file name
%
fileName = 'materialsList.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);


%------------------------------ ELEMENTS -----------------------------

% Create object handling gauss point integration
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager 
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_CubeMCC');


%------------------------ BOUNDARY CONDITIONS ------------------------

% UTILITY TO WRITE BC FILES 
writeBCfiles('BCs/dirPoroBottom','NodeBC','Dir',{'Poro','x','y','z'},'bottom_fix',0,0,topology,2);
writeBCfiles('BCs/neuPoroZ','SurfBC','Neu',{'Poro','z'},'Top_LoadZ',vectimeF,vecF,topology,3);

% Collect BC input file in a list
fileName = ["BCs/dirPoroBottom.dat","BCs/neuPoroZ.dat"];

% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);


%------------------------ INITIAL STATE ------------------------
% GeoParameters
nu=mat.db(1).ConstLaw.nu;
depth=-5; % m
gamma=20; % kN/m^3  

% Cells Initial State
cellslist=dofmanager.getFieldCells("Poromechanics");
for i=1:length(cellslist)
    cID=cellslist(i);
    z=depth+elems.cellCentroid(cID,3); %Real Depth

    % Initial State
    sigmaZ=gamma*z;
    sigmaX=M1*sigmaZ;
    sigmaY=M2*sigmaZ;
    state.iniStress(cID,1)=sigmaX;
    state.iniStress(cID,2)=sigmaY;
    state.iniStress(cID,3)=sigmaZ;

    % Average Preconsolidation Stress pc
    OCR=1.5;
    pc(i)=OCR*1/3*(sigmaX+sigmaY+sigmaZ);
end

if elems.mesh.cellVTKType(1)==10 % element is tetra
    nptGauss=1;
elseif elems.mesh.cellVTKType(1)==12 % element is exa
    nptGauss=8;
end

if isa(mat.db(1).ConstLaw, 'ModifiedCamClay')
        pc=repmat(pc',1,nptGauss);
        newmap=mat.db(1);
        newmap.ConstLaw.pc=-pc;
        mat.db(1)=newmap;
else
end
state.conv.stress = state.iniStress;
state.curr.stress = state.iniStress;
fid=fopen('HN_p_c.txt','a');
fprintf(fid,'\n%.6f',mat.db(1).ConstLaw.pc(349));

% Print model initial state
printUtils.printState(linSyst,state);

% ---------------------------- SOLUTION -------------------------------
fprintf('\n MCC Model\n\n')

% Create the object handling the (nonlinear) solution of the problem
solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
[simState, endState] = solver.NonLinearLoop();


% Finalize the print utility
printUtils.finalize()

%% ------------------------ POST PROCESSING ----------------------------

% Recognizing Axis (0 0 0) + <0 0 1> Elements
if strcmp(meshName,"Mesh/CubeTetra5.msh")
    elem_x_OK=find(elems.cellCentroid(:,1)==0.0625);
    elem_y_OK=find(elems.cellCentroid(:,2)==0.0625);
elseif strcmp(meshName,"Mesh/CubeTetra7.msh")
    elem_x_OK=find(round(elems.cellCentroid(:,1),4)==0.0417);
    elem_y_OK=find(round(elems.cellCentroid(:,2),4)==0.0417);
elseif strcmp(meshName,"Mesh/CubeTetra3.msh")
    elem_x_OK=find(elems.cellCentroid(:,1)==0.1250);
    elem_y_OK=find(elems.cellCentroid(:,2)==0.1250);
elseif strcmp(meshName,"Mesh/CubeExa5.msh") 
    elem_x_OK=find(round(elems.cellCentroid(:,1),4)==0.125);
    elem_y_OK=find(round(elems.cellCentroid(:,2),4)==0.125);
elseif strcmp(meshName,"Mesh/CubeExa7.msh")
    elem_x_OK=find(round(elems.cellCentroid(:,1),4)==0.0833);
    elem_y_OK=find(round(elems.cellCentroid(:,2),4)==0.0833);
end
    elem_axis=intersect(elem_y_OK,elem_x_OK);
    highest_elem=intersect(find(elems.cellCentroid(:,3)==max(elems.cellCentroid(elem_axis,3))),elem_axis);

% Sorting Elements along axis
Z_coords=elems.cellCentroid(elem_axis,3);
[sortZ_coords,id]=sort(Z_coords);
% Sorting SigmaZ
sortelem_axis=elem_axis(id);
% sigmaz_sorted_t0=state.conv.stress(sortelem_axis,3);
% sigmaz_sorted_tend=endState.conv.stress(sortelem_axis,3);
axis_sigmaz=load('axis_sigma_z.txt');
axis_epsz=load('axis_eps_z.txt');

figure (7)
tiledlayout(1,2)
nexttile % Plot Sigma
for i=1:4:simParam.tMax
    plot(-axis_sigmaz(i*length(sortelem_axis)-(length(sortelem_axis)-1):i*length(sortelem_axis)),sortZ_coords,'*-',LineWidth=1.5)
    hold on
end
xlabel('σ_z [kPa]','FontSize',18)
ylabel('local z [m]','FontSize',18)
title('σ_z lungo asse verticale','FontSize',25)
legend('t = 0','t = 5','t = 10','t = 15','t = 20','FontSize',13)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on
hold off

nexttile
for i=1:4:simParam.tMax
    plot(-axis_epsz(i*length(sortelem_axis)-(length(sortelem_axis)-1):i*length(sortelem_axis)),sortZ_coords,'*-',LineWidth=1.5)
    hold on
end
xlabel('ε_z [kPa]','FontSize',18)
ylabel('local z [m]','FontSize',18)
title('ε_z lungo asse verticale','FontSize',25)
legend('t = 0','t = 5','t = 10','t = 15','t = 20','FontSize',13)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on
hold off

% % Plot Uz over time
node_x_OK=find(elems.mesh.coordinates(:,1)==0);
node_y_OK=find(elems.mesh.coordinates(:,2)==0);
node_z_OK=find(elems.mesh.coordinates(:,3)==0.5);
node_xy_OK=intersect(node_x_OK,node_y_OK);
central_Node=intersect(node_xy_OK,node_z_OK);
Uz_centralNode=printUtils.results.expDispl(central_Node*3,:);
Elastic_Uz_centralNode=Elastic_Disp(central_Node*3,:);
Time=0:simParam.tMax;

figure (2)
tiledlayout(1,2)
nexttile
plot(Time,Uz_centralNode,'r',LineWidth=1.5)
hold on
plot(Time,Elastic_Uz_centralNode,'b--',LineWidth=1.5)
grid on
xlabel('Tempo [s]','FontSize',18)
ylabel('u_z [m]','FontSize',18)
title('Spostamenti u_z nel tempo per il nodo [0, 0, 0.5]','FontSize',18)
legend('MCC','Elastic','FontSize',18)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
hold off

% Recognizing Axis (0 0 0) + <0 0 1> Nodes

node_x_OK=find(abs(0-elems.mesh.coordinates(:,1))<1.e-5);
node_y_OK=find(abs(0-elems.mesh.coordinates(:,2))<1.e-5);
node_axis=intersect(node_y_OK,node_x_OK);

% Sorting Nodes along axis
Z_coords=elems.mesh.coordinates(node_axis,3);
[sortZ_coords,id]=sort(Z_coords);
% Sorting UZ
sortnode_axis=node_axis(id);
U=printUtils.results.expDispl;
i=0;
for j=1:5:size(U,2)
    i=i+1;

    ux_sorted_t(:,i)=U(3*sortnode_axis-2,j);
    uy_sorted_t(:,i)=U(3*sortnode_axis-1,j);
    uz_sorted_t(:,i)=U(3*sortnode_axis,j);  
end

figure (3)
tiledlayout(1,2)
nexttile
 % Plot U_x
for i=1:size(uz_sorted_t,2)
    plot(ux_sorted_t(:,i),sortZ_coords,'-*',LineWidth=1.5)
    hold on
end
title('Spostamenti u_x lungo asse verticale','FontSize',22)
ylabel('z locale [m]','FontSize',18)
xlabel('u_x [m]','FontSize',18)
legend('t = 0','t = 5','t = 10','t = 15','t = 20','FontSize',13)
grid on
hold off

nexttile % Plot U_y
for i=1:size(uz_sorted_t,2)
    plot(uy_sorted_t(:,i),sortZ_coords,'-*',LineWidth=1.5)
    hold on
end
title('Spostamenti u_y lungo asse verticale','FontSize',22)
ylabel('z locale [m]','FontSize',18)
xlabel('u_y [m]','FontSize',18)
legend('t = 0','t = 5','t = 10','t = 15','t = 20','FontSize',13)
grid on
hold off

figure (2)
nexttile(2) % Plot U_z
for i=1:size(uz_sorted_t,2)
    plot(uz_sorted_t(:,i),sortZ_coords,'-*',LineWidth=1.5)
    hold on
end
title('Spostamenti u_z lungo asse verticale','FontSize',18)
ylabel('z locale [m]','FontSize',18)
xlabel('u_z [m]','FontSize',18)
legend('t = 0','t = 5','t = 10','t = 15','t = 20','FontSize',13)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on
hold off

if isa(mat.db(1).ConstLaw, 'ModifiedCamClay')
    figure (4)
    t=tiledlayout(1,2);
    nexttile
    tetramesh(elems.mesh.cells, elems.mesh.coordinates,-pc, 'EdgeColor', [0 0 0],'FaceAlpha',0.6);
    colormap(flipud(parula(128)));
    colorbar;
    hold on
    axis equal tight;
    xlabel('coordinate x'); ylabel('coordinate y'); zlabel('coordinate z');
    title('Pressione di preconsolidazione p_c (Inizio simulazione)','FontSize',18);
    clim([-0.8*min(pc) 0.92*max(mat.db(1).ConstLaw.pc)]);
    camproj('perspective');
    hold off

    nexttile
    tetramesh(elems.mesh.cells, elems.mesh.coordinates,mat.db(1).ConstLaw.pc, 'EdgeColor', [0 0 0],'FaceAlpha',0.6);
    colormap(flipud(parula(128)));
    colorbar;
    hold on
    axis equal tight;
    xlabel('coordinate x'); ylabel('coordinate y'); zlabel('coordinate z');
    title('Pressione di preconsolidazione p_c (Fine simulazione)','FontSize',18);
    clim([-0.8*min(pc) 0.92*max(mat.db(1).ConstLaw.pc)]);
    camproj('perspective');
    hold off
end

sigmaz_highelem=load('HN_sigma_z.txt');
epsilonz_highelem=load('HN_eps_z.txt');
Elastic_sigmaz_highelem=load('Elastic_HN_sigma_z.txt');
Elastic_epsilonz_highelem=load('Elastic_HN_eps_z.txt');

figure (5)
tiledlayout(2,2)
nexttile(1)   
plot(0:simParam.tMax,-sigmaz_highelem(end-simParam.tMax:end),'r',LineWidth=1.5)
hold on
plot(0:simParam.tMax,-Elastic_sigmaz_highelem(end-simParam.tMax:end),'b--',LineWidth=1.5)
ylabel('σ_z [kPa]','FontSize',18)
xlabel('Tempo [s]','FontSize',18)
title('σ_z nel tempo per il FE in sommità','FontSize',25)
legend('MCC','Elastic','FontSize',15)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on

nexttile(3)
plot(0:simParam.tMax,-epsilonz_highelem(end-simParam.tMax:end),'r',LineWidth=1.5)
hold on
plot(0:simParam.tMax,-Elastic_epsilonz_highelem(end-simParam.tMax:end),'b--',LineWidth=1.5)
ylabel('ε_z [-]','FontSize',18)
xlabel('Tempo [s]','FontSize',18)
title('ε_z nel tempo per il FE in sommità','FontSize',25)
legend('MCC','Elastic','FontSize',15)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on

nexttile([2 1])
plot(-epsilonz_highelem(end-simParam.tMax:end),-sigmaz_highelem(end-simParam.tMax:end),'r',LineWidth=1.5)
hold on
plot(-Elastic_epsilonz_highelem(end-simParam.tMax:end),-Elastic_sigmaz_highelem(end-simParam.tMax:end),'b--',LineWidth=1.5)
ylabel('σ_z [kPa]','FontSize',18)
xlabel('ε_z [-]','FontSize',18)
title('σ_z - ε_z per il FE in sommità','FontSize',25)
legend('MCC','Elastic','FontSize',15)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on

pc_HN=load('HN_p_c.txt');
figure (6)
plot(0:length(pc_HN)-1,pc_HN,'r',LineWidth=1.5)
hold on
ylim([0.8*min(pc_HN) 1.2*max(pc_HN)])
xlabel('Tempo [s]','FontSize',18)
ylabel('p_c [kPa]','FontSize',18)
title('Tensione di preconsolidazione p_c nel tempo','FontSize',25)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on 
hold off

figure (8)
P=load("HN_p.txt");
Q=load("HN_q.txt");
plot(-P,Q,'r-',LineWidth=1.5)
hold on
xlabel('Tensione volumetrica p [kPa]','FontSize',18)
ylabel('Tensione deviatorica q [kPa]','FontSize',18)
title('Diagramma p-q','FontSize',25)
grid on 
%     % CSL critical state line
    pp=linspace(0,pc_HN(end)/2,250);
    qq=@(pp) pp.*mat.db(1).ConstLaw.M;
    plot(pp,qq(pp),'--k')
        % yielding surface
        pp=linspace(0,pc_HN(1),250);    
        qq=@(pp) sqrt(-pp.*(pp-pc_HN(1)).*mat.db(1).ConstLaw.M^2);
        plot(pp,qq(pp),'b',LineWidth=1)
for j=1:length(pc_HN)
    if mod(j,simParam.tMax+1)==0
        pcyield=pc_HN(j);
        pp=linspace(0,pcyield,250);    
        qq=@(pp) sqrt(-pp.*(pp-pcyield).*mat.db(1).ConstLaw.M^2);
        plot(pp,qq(pp),'b--',LineWidth=0.5)
    else
    end
end
legend('Curva p-q','CSL','F=0 iniziale','F=0 finale')

figure (1)
plot(vectimeF,-vecF,'r',LineWidth=1.5)
hold on
xlabel('Tempo [s]','FontSize',18);
ylabel('σ_z [kPa]','FontSize',18)
title('Entità del carico applicato q nel tempo','FontSize',25)
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on