function domainMesh = read_mesh(filename)

  mesh_file = fopen(filename);  

  % if it fails create new mesh
  if mesh_file == -1
        error("Mesh file not found");
  end

  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domainMesh.type = sscanf(fgets(mesh_file),'%s');

  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.num_nodes=nodes_number; 

  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end

  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.num_elements=polygons_number;  

  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end

  % read indices of nodes on the boundary
  line = fgets(mesh_file); % read commented line
  bound = fgets(mesh_file); % read commented line
  domainMesh.boundary_nodes = sscanf(bound,'%d'); %Left

  fclose(mesh_file);
  
end

