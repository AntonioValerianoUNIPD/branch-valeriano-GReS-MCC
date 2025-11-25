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

F=-10;
timeF=5;
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
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_CubeElastic');


%------------------------ BOUNDARY CONDITIONS ------------------------
    
% UTILITY TO WRITE BC FILES 
writeBCfiles('BCs/dirPoroBottom','NodeBC','Dir',{'Poro','x','y','z'},'bottom_fix',0,0,topology,2);
writeBCfiles('BCs/neuPoroZ','SurfBC','Neu',{'Poro','z'},'Top_LoadZ',[0 timeF],[0 F],topology,3);
% writeBCfiles('BCs/neuPoroY','SurfBC','Neu',{'Poromechanics','y'},'Lat_LoadY',[0 0.1],[0 -100],topology,2);
% writeBCfiles('BCs/neuPoroX','SurfBC','Neu',{'Poromechanics','x'},'Lat_LoadX',[0 0.1],[0 -100],topology,2);

% Collect BC input file in a list
fileName = ["BCs/dirPoroBottom.dat","BCs/neuPoroZ.dat"];
% fileName = ["BCs/dirPoroBottom.dat","BCs/neuPoroZ.dat","BCs/neuPoroY.dat","BCs/neuPoroX.dat"];

% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);


%------------------------ INITIAL STATE ------------------------
% GeoParameters
nu=mat.db(1).ConstLaw.nu;
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

% Print model initial state
printUtils.printState(linSyst,state);

% ---------------------------- SOLUTION -------------------------------

fprintf('\n Elastic Model\n\n')
% Create the object handling the (nonlinear) solution of the problem
solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
[Elastic_simState, Elastic_endState] = solver.NonLinearLoop();


% Finalize the print utility
Elastic_printUtils.finalize()
Elastic_Disp=Elastic_printUtils.results.expDispl;
