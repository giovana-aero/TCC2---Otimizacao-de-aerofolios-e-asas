# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:58:37 2021

@author: Guga Weffort
"""

import openvsp as vsp

# Reiniciar o modelo 
vsp.ClearVSPModel()

print( "Begin Airfoil Export Script\n" ) 
    
print( "--> Creating Wing Geometry\n" ) 

# //==== Add Wing Geometry ====//
wing_id = vsp.AddGeom( "WING" )

# //==== Turn Symmetry OFF ====//
vsp.SetParmVal( vsp.FindParm( wing_id, "Sym_Planar_Flag", "Sym"), 0.0 )  #// Note: bool input not supported in SetParmVal

# // Note: The shape of this wing was selected for demonstration purposes only

# //===== Add XSec ====//
vsp.InsertXSec( wing_id, 1, vsp.XS_SIX_SERIES )

# //===== Change Some Parameters 1st Section ====//
vsp.SetParmVal( wing_id, "Twist", "XSec_1", 20.0 )
vsp.SetParmVal( wing_id, "InLEMode", "XSec_1", vsp.BLEND_MATCH_OUT_LE_TRAP )
vsp.SetParmVal( wing_id, "InTEMode", "XSec_1", vsp.BLEND_MATCH_OUT_TE_TRAP )

vsp.SetParmVal( wing_id, "ThickChord", "XSecCurve_0", 0.2 )

# //===== Define Number of Airfoils to Interpolate in Section 1 ====//
vsp.SetParmVal( wing_id, "SectTess_U", "XSec_1", 9 )

vsp.Update()

# //===== Change Some Parameters 2nd Section ====//
vsp.SetParmVal( wing_id, "Span", "XSec_2", 3.0 )
vsp.SetParmVal( wing_id, "Root_Chord", "XSec_2", 1.5 )
vsp.SetParmVal( wing_id, "Tip_Chord", "XSec_2", 0.5 )
vsp.SetParmVal( wing_id, "Sweep", "XSec_2", 50.0 )
vsp.SetParmVal( wing_id, "Twist", "XSec_2", 20.0 )
vsp.SetParmVal( wing_id, "Dihedral", "XSec_2", 20.0 )

vsp.Update()

#//===== Set Airfoil Export Parms in Vehicle Parm Container ====//
veh_id = vsp.FindContainer( "Vehicle", 0 )

vsp.SetParmVal( vsp.FindParm( veh_id, "AFWTessFactor", "AirfoilExport" ), 2.0 ) #// double W tesselation for export
vsp.SetParmVal( vsp.FindParm( veh_id, "AFAppendGeomIDFlag", "AirfoilExport" ), 1.0 ) #// Note: bool input not supported in SetParmVal

vsp.Update()

# //===== Export All Airfoils and Metadata *.csv ====//
print(  "--> Exporting Airfoil Files\n" ) 

vsp.ExportFile( "Airfoil_Metadata.csv", vsp.SET_ALL, vsp.EXPORT_SELIG_AIRFOIL )

#//==== Check For API Errors ====//
# while vsp.GetNumTotalErrors() > 0:
#     err = vsp.PopLastError()
#     print( err.GetErrorString() )


# Exportar arquivo do openvsp
vsp.WriteVSPFile("AirfoilExport.vsp3")

print("COMPLETE")