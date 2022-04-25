# Objetivo: gerar uma asa e simular
# Descrição desta versão: tradução do script Wing.vspscript pro python

import openvsp as vsp


name = 'wing_test_sim1'

# Reiniciar o modelo 
vsp.ClearVSPModel()


# //==== Create A Multi Section Wing and Change Some Parameters ====//

#//==== Add Wing ====//
wid = vsp.AddGeom( "WING", "" )

#//===== Insert A Couple More Sections =====//
vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )
vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )
vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )

#//===== Cut The Original Section =====//
vsp.CutXSec( wid, 1 )

#//===== Change Driver =====//
# isto funciona da mesma forma que no gui: o projetista determina três parâmetros
# entre os 8 possíveis pra manipular a planta da asa
vsp.SetDriverGroup( wid, 1, vsp.AREA_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )

vsp.SetParmVal( wid, "RotateAirfoilMatchDideralFlag", "WingGeom", 1.0 )

#//===== Change Some Parameters 1st Section ====//
vsp.SetParmVal( wid, "Root_Chord", "XSec_1", 7.0 )
vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", 3.0 )
vsp.SetParmVal( wid, "Area", "XSec_1", 45.0 )
vsp.SetParmVal( wid, "Sweep", "XSec_1", 40.0 )

#//==== Because Sections Are Connected Change One Section At A Time Then Update ====//
vsp.Update()

#//===== Change Some Parameters 2nd Section ====//
vsp.SetParmVal( wid, "Tip_Chord", "XSec_2", 2.0 )
vsp.SetParmVal( wid, "Sweep", "XSec_2", 60.0 )
vsp.SetParmVal( wid, "Dihedral", "XSec_2", 30.0 )
vsp.Update()

#//===== Change Some Parameters 3rd Section ====//
vsp.SetParmVal( wid, "Sweep", "XSec_3", 60.0 )
vsp.SetParmVal( wid, "Dihedral", "XSec_3", 80.0 )
vsp.Update()

#//==== Change Airfoil ====//
vsp.SetParmVal( wid, "Camber", "XSecCurve_0", 0.02 )
vsp.Update()


# //==== Check For API Errors ====//
# while ( GetNumTotalErrors() > 0 )
# {
#     ErrorObj err = PopLastError();
#     Print( err.GetErrorString() );
#    }
errorMgr = vsp.ErrorMgrSingleton_getInstance()
num_err = errorMgr.GetNumTotalErrors()
for i in range(0, num_err):
    err = errorMgr.PopLastError()
    print("error = ", err.m_ErrorString)

# Exportar arquivo do openvsp
vsp.WriteVSPFile( name + ".vsp3")