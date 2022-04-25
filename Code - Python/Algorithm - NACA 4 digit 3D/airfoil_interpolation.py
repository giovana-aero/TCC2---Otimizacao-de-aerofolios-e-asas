import numpy as np
from airfoil_rotation import airfoil_rotation
from numpy.matlib import repmat

def airfoil_interpolation(coo_r,coo_t,op,d,tw):
    # Encontrar o bordo de ataque de ambos os aerofólios
    delta_r = np.zeros(coo_r.shape[0]-1)
    for i in range(coo_r.shape[0]-1):
        delta_r[i] = coo_r[i+1,0] - coo_r[i,0]
    
    
    delta_b_r = np.array([int(x) for x in delta_r >= 0])  #int(delta_r >= 0)
    #for i in range(len(delta_b_r)):
     #   if delta_b_r[i] == 1: LE = i;break
    LE = np.argwhere(delta_b_r == 1)[0,0]
    
    # Separar superfícies (ordenadas)
    ex_r = coo_r[0:LE+1,:]
    in_r = coo_r[LE:,:]
    ex_t = coo_t[0:LE+1,:]
    in_t = coo_t[LE:,:]
    
    # Fazer interpolações
    # 0 -> raiz
    # 1 -> ponta
    ex_intp = (ex_r*(1-op) + ex_t*op)
    in_intp = (in_r*(1-op) + in_t*op)
    
    # Montar coordenadas
    coo_intp = np.vstack((ex_intp,in_intp[1:,:]))
    
    # Fazer rotação do perfil
    if tw != 0:coo_intp = airfoil_rotation(coo_intp,tw)
    
    if d != 0: # Passar as ordenadas pro eixo z e aplicar um deslocamento ao longo do eixo y
    #    n = size(coo_intp,1);
        coo_intp = np.hstack((np.array([coo_intp[:,0]]).transpose(),repmat(d,coo_intp.shape[0],1),np.array([coo_intp[:,1]]).transpose()))
        
    #NODE(1:sec_N,:) = [coo_t(:,1),repmat(-b/2,sec_N,1),coo_t(:,2)];
    # else # Retornar as coordenadas no formato 2D comum
    

    return coo_intp





