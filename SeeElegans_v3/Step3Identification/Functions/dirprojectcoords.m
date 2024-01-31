function dircoords=dirprojectcoords(anterdir,dorsaldir,leftdir,midpoint,coords)
dircoords(1)=dirproject(anterdir,midpoint,coords);    
dircoords(2)=dirproject(dorsaldir,midpoint,coords); 
dircoords(3)=dirproject(leftdir,midpoint,coords);


