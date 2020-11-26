# Grid Refinement Functions

Make_Grid <- function(x_coords, y_coords, cell_size, buffer){
  #check for same length in x and y 
  if (length(x_coords) != length(y_coords)){
    return('The x and y coordinates do not have the same number of values')
  }
  # calculate the origin
  x_origin = floor(min(x_coords)) - (cell_size * buffer)
  y_origin = floor(min(y_coords)) - (cell_size * buffer)
  # calculate the farthest distance out
  x_dist = ceiling(max(x_coords)) + (cell_size * buffer)
  y_dist = ceiling(max(y_coords)) + (cell_size * buffer)
  # grid check, if last value in the sequence is not greater than distance in x and y need to add one more cell size to vector
  x_seq <- seq(x_origin ,x_dist, cell_size)
  if (max(x_seq) < (x_dist - cell_size/2)){
    x_seq <- c(x_seq, max(x_seq) + cell_size)
  }
  y_seq <- seq(y_origin ,y_dist, cell_size)
  if (max(y_seq) < (y_dist - cell_size/2)){
    y_seq <- c(y_seq, max(y_seq) + cell_size)
  }
  # create grid data set
  grid_data <- CJ(x = x_seq, y = y_seq)
  # add cell information
  grid_data$cell_size = cell_size
  grid_data$cell_id = 1:nrow(grid_data)
  grid_data
}

Ghost_Nodes <- function(Grid, x_coords, y_coords, ref_method){
  # method selection
  if (ref_method == 'Single') {
    Ghosts <- lapply(1:length(x_coords), function(a){
      temp_point_x <- x_coords[a]
      temp_point_y <- y_coords[a]
      temp_cell_size <- min(Grid$cell_size)
      # locate cell ids
      Grid_Nodes <- Grid[between(x = x, lower = temp_point_x-temp_cell_size/2, upper = temp_point_x+temp_cell_size/2, incbounds = T)][between(x = y, lower = temp_point_y-temp_cell_size/2, upper = temp_point_y+temp_cell_size/2)]
      
      Grid_Nodes <- Grid_Nodes[cell_size == min(Grid_Nodes$cell_size)][1]
      
      Grid_Nodes
    })
  }else{
    
    if(ref_method != 'Radial'){return("No method found with that name. Please check spelling of the reference method")}
    # run for each point pair
    Ghosts <- lapply(1:length(x_coords), function(a){
      temp_point_x <- x_coords[a]
      temp_point_y <- y_coords[a]
      temp_cell_size <- min(Grid$cell_size)
      
      #find max/min x and y by multiplying the point location by 1.5 the current minimum distance in the grid
      min_x = temp_point_x - temp_cell_size*1.5
      min_y = temp_point_y - temp_cell_size*1.5
      max_x = temp_point_x + temp_cell_size*1.5
      max_y = temp_point_y + temp_cell_size*1.5
      
      # extract cell ids that will be needed
      Grid_Nodes <- Grid[between(x, min_x, max_x)][between(y, min_y, max_y)]
      Grid_Nodes
    })
  }
  Grid_Nodes <- rbindlist(Ghosts)
  Grid_Nodes
}

Quad_Tree <- function(Grid, ghost_nodes){
  # this is the outer loop that controls the number of refinements to do
  # using ghost nodes refine grid down number of steps desired
  new_cells <- lapply(ghost_nodes$cell_id, function(b) {
    temp_grid <- Grid[cell_id == b]
    # reduce each side by 2
    new_cell1 = data.table(x = temp_grid$x + temp_grid$cell_size/4, y = temp_grid$y + temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell2 = data.table(x = temp_grid$x + temp_grid$cell_size/4, y = temp_grid$y - temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell3 = data.table(x = temp_grid$x - temp_grid$cell_size/4, y = temp_grid$y - temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cell4 = data.table(x = temp_grid$x - temp_grid$cell_size/4, y = temp_grid$y + temp_grid$cell_size/4, cell_size = temp_grid$cell_size/2, cell_id = 0)
    new_cells <- rbindlist(list(new_cell1,new_cell2, new_cell3, new_cell4))
  })
  # remove old cells and add in new cells, order by x and then y
  new_cells <- rbindlist(new_cells)
  Grid <- Grid[!cell_id %in% ghost_nodes$cell_id]
  Grid <- rbind(Grid, new_cells)
  Grid <- Grid[order(x,y)]
  Grid$cell_id <- 1:nrow(Grid)
  Grid
}

Grid_Refinement_Function <- function(x_coords, y_coords, cell_size, buffer, num_ref, ref_method) {
  temp_grid <- Make_Grid(x_coords = x_coords, y_coords = y_coords, cell_size = cell_size, buffer = buffer)
  # check if refinement number is zero and return normal grid
  if (num_ref > 0){
    # now we enter into the refine grid loop
    for (i in 1:num_ref){
      # first find the ghost nodes
      ghosts <- Ghost_Nodes(Grid = temp_grid, x_coords = x_coords, y_coords = y_coords, ref_method = ref_method)
      # then run quad tree refinement
      temp_grid <- Quad_Tree(Grid = temp_grid, ghost_nodes = ghosts)
    }
  }
  # if no more refinements return the updated grid
  temp_grid
}