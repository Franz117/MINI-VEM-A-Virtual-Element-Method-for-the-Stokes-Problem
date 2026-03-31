function [dof_U1,dof_U2,dof_P] = el_to_dof(ne,nv,connect,ie)
    
    dof_U1 = [connect; nv+ie];
    dof_U2 = [connect+nv+ne; 2*nv+ ne + ie];
    dof_P  = connect + 2*nv + 2*ne;

end

