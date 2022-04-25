function aero = run_apame(pop,dat)

% Esta função escreve o arquivo de entrada do APAME, faz a simulação da asa
% e retorna os coeficientes aerodinâmicos relevantes

% Dados da simulação
cases = dat.cases;
v_ref = dat.v_ref;
rho = dat.rho;
p_atm = dat.p_atm;
aoa = dat.aoa;
mach = dat.mach;
% Dados da geometria
b = pop.b;
mac = pop.mac;
NODE = pop.NODE;
PANEL = pop.PANEL;
sec_N = pop.sec_N;
S = pop.S;


aero = zeros(cases,4);
for i = 1:cases
    
    % Imprimir o arquivo de entrada ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fid = fopen('apame_input.inp','w');
    fprintf(fid,'APAME input file\nVERSION 3.1\n\n\n');
    % Parâmetros do escoamento
    fprintf(fid,'# FLOW PARAMETERS\n\n');
    fprintf(fid,'# Airspeed [m/s]\n');
    fprintf(fid,'AIRSPEED %f\n',v_ref(i));
    fprintf(fid,'# Air density [kg/m^3\n');
    fprintf(fid,'DENSITY %.20f\n',rho(i));
    fprintf(fid,'# Atmospheric pressure [Pa]\n');
    fprintf(fid,'PRESSURE %.10f\n',p_atm(i));
    fprintf(fid,'# Prandtl-glauert corretion:\n#  0-no correction\n#  *-Mach number\n');
    fprintf(fid,'MACH %f\n',mach(i));
    fprintf(fid,'# Number of cases\n# Angle of attack [degrees]\n# sideslip angle [degrees]\n');
    fprintf(fid,'CASE_NUM 1\n')
    fprintf(fid,'%f\n0\n\n\n',aoa(i));
    % Valores de referência
    fprintf(fid,'# REFERENCE VALUES\n\n');
    fprintf(fid,'# Wingspan [m]\n');
    fprintf(fid,'WINGSPAN %f\n',b);
    fprintf(fid,'# Mean aerodynamic chord [m]\n');
    fprintf(fid,'MAC %f\n',mac);
    fprintf(fid,'# Wing surface area [m^2]\n');
    fprintf(fid,'SURFACE %f\n',S);
    fprintf(fid,'# Reference point [m]\n');
    fprintf(fid,'ORIGIN *\n0 0 0\n\n\n');
    % Parâmetros do solucionador
    fprintf(fid,'# SOLVER PARAMETERS\n\n');
    fprintf(fid,'# Singularity method:\n#  0-constant source/doublet\n#  1-constant doublet\n');
    fprintf(fid,'METHOD 0\n');
    fprintf(fid,'# Error\n');
    fprintf(fid,'ERROR 0.0000001\n');
    fprintf(fid,'# Collocation point depth\n');
    fprintf(fid,'COLLDIST 0.0000001\n');
    fprintf(fid,'# "Far field" coefficient\n');
    fprintf(fid,'FARFIELD 5\n');
    fprintf(fid,'# Collocation point calculation:\n#  0-approximative\n#  1-accurate\n');
    fprintf(fid,'COLLCALC 0\n');
    fprintf(fid,'# Interpolation method/order for velocity calculations:\n#  0-nodal\n#  1-first\n#  2-second\n');
    fprintf(fid,'VELORDER 1\n\n\n');
    % Resultados
    fprintf(fid,'# RESULT REQUESTS\n#  0-no\n# 1-yes\n');
    fprintf(fid,'RESULTS 1\n');
    fprintf(fid,'#  1 - coefiicients\n');
    fprintf(fid,'RES_COEF 1\n');
    fprintf(fid,'#  2 - forces\n');
    fprintf(fid,'RES_FORC 0\n');
    fprintf(fid,'#  3 - geometry\n');
    fprintf(fid,'RES_GEOM 0\n');
    fprintf(fid,'#  4 - velocity\n');
    fprintf(fid,'RES_VELO 0\n');
    fprintf(fid,'#  5 - pressure\n');
    fprintf(fid,'RES_PRES 0\n');
    fprintf(fid,'#  6 - center points\n');
    fprintf(fid,'RES_CENT 0\n');
    fprintf(fid,'#  7 - dipole values\n');
    fprintf(fid,'RES_DOUB 0\n');
    fprintf(fid,'#  8 - source values\n');
    fprintf(fid,'RES_SORC 0\n');
    fprintf(fid,'#  9 - velocity components\n');
    fprintf(fid,'RES_VELC 0\n');
    fprintf(fid,'# 10 - mesh characteristics\n');
    fprintf(fid,'RES_MESH 0\n');
    fprintf(fid,'# 11 - static pressure\n');
    fprintf(fid,'RES_STAT 0\n');
    fprintf(fid,'# 12 - dynamic pressure\n');
    fprintf(fid,'RES_DYNA 0\n');
    fprintf(fid,'# 13 - manometer pressure\n');
    fprintf(fid,'RES_MANO 0\n\n\n');
    % Geometria
    fprintf(fid,'# GEOMETRY\n\n');
    fprintf(fid,'NODES %d\n',size(NODE,1));
    fprintf(fid,'%.10f %.10f %.10f\n',NODE');
    fprintf(fid,'\nPANELS %d\n',size(PANEL,1));
    fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',PANEL(1:end-sec_N+1,:)');
    fprintf(fid,'%d %d %d %d %d %d %d\n',PANEL(end-sec_N+2:end,1:7)');
    fclose(fid);

    % Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    % Nota: para fazer a simulação, é necessário abrir o executável e inserir o nome
    % do arquivo a ser analisado sem a extensão .inp. Devido à natureza desse processo,
    % é necessário o uso de um arquivo permanente contendo o nome do arquivo de input
    % do apame. A omissão da exetnsão do arquivo permanente foi realizada com o 
    % intuito de identificar melhor esse arquivo.
    system('apame_win64.exe < apame_cmd_input');clc

    % Apagar arquivos desnecessários
    delete('fort.2');
    delete('apame_input.res');

    % Ler os resultados ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fidpolar = fopen('apame_input.log');
    dataBuffer1 = textscan(fidpolar,'%f %f %f %f %f %f','HeaderLines',51,... % CX, CY, CZ, CL, CM, CN (body reference frame)
        'CollectOutput',1,...
        'Delimiter',''); fclose(fidpolar);
    fidpolar = fopen('apame_input.log');
    dataBuffer2 = textscan(fidpolar,'%f %f %f','HeaderLines',55,... % CD, CK, CL (aerodynamic reference frame)
        'CollectOutput',1,...
        'Delimiter',''); fclose(fidpolar);

    CL = dataBuffer2{1,1}(1,3);
    CD = dataBuffer2{1,1}(1,1);
    CM = dataBuffer1{1,1}(1,5);
        
    if isnan(CL) || CD <= 0
        aero = 'n'; break
    else
        aero(i,:) = [CL,CD,CL/CD,CM];
    end
    
end

%delete('apame_input.log');

end