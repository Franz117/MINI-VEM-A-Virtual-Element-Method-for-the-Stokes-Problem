function [dof_U1,dof_U2,dof_P] = node_to_dof(node,ne,nv)
    dof_U1 = node;
    dof_U2 = node+ne+nv;
    dof_P = node+2*ne+2*nv;

end

