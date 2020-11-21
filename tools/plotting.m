function [] = plotting()

  clc; more off;

  delimiterIn   = ' ';
  headerlinesIn = 1;
    
  files_mesh      = 'MESH_*.dat';
  files_solution  = 'SOLUTION_*.dat';
  
  %------------------------------% 
  % Importing the MESH DATA
  %------------------------------%
  import_mesh = dir(files_mesh);
  nfiles_mesh = length(import_mesh);
  for i=1:nfiles_mesh
    if (i==1)
      filename = import_mesh(i).name;
      mydata_mesh = importdata(filename,delimiterIn,headerlinesIn);
      x = mydata_mesh.data(:,1);
      xIni = min(x);
      xEnd = max(x);
    end
    if (i==2)
      filename = import_mesh(i).name;
      mydata_mesh = importdata(filename,delimiterIn,headerlinesIn);
      y = mydata_mesh.data(:,1);
      yIni = min(y);
      yEnd = max(y);
    end
  end
  [X,Y] = meshgrid(x,y);

  %------------------------------% 
  % Importing the SOLUTION DATA
  %------------------------------%
  import_solution = dir(files_solution);
  nfiles_solution = length(import_solution);
  for i=1:nfiles_solution
    filename = import_solution(i).name;
    mydata_solution = importdata(filename,delimiterIn,headerlinesIn);
    u = mydata_solution.data(:,1);
    uMin = min(u);
    uMax = max(u);
    
    aFile = sprintf('%s',filename);
    fprintf('Plotting file %s\n',aFile);
    U = reshape(u,[size(X,2),size(Y,1)]);
    U = transpose(U);
    contourf(X,Y,U);

    xlim([xIni,xEnd]);
    ylim([yIni,yEnd]);
    pbaspect([1 1 1]);
    daspect([1 1 1]);
    pause;
  end

end
