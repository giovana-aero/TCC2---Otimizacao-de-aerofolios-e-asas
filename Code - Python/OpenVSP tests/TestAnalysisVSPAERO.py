
import openvsp as vsp
# import degen_geom as DG
# import utilities
# import charm

vsp.ClearVSPModel()

#//==== Create some test geometries ====// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(  "--> Generating Geometries" ) 
# Print( string( "" ) );

pod_id = vsp.AddGeom( "POD", "" )
wing_id = vsp.AddGeom( "WING", "" )

vsp.SetParmVal( wing_id, "X_Rel_Location", "XForm", 2.5 )
vsp.SetParmVal( wing_id, "TotalArea", "WingGeom", 25 )

subsurf_id = vsp.AddSubSurf( wing_id, vsp.SS_CONTROL, 0 )

vsp.Update()

# //==== Setup export filenames ===//
fname_vspaerotests = "TestVSPAero.vsp3"

# //==== Save Vehicle to File ====//
print( "-->Saving vehicle file to: ")
print( fname_vspaerotests)
# Print( "" );
vsp.WriteVSPFile( fname_vspaerotests, vsp.SET_ALL )
print( "COMPLETE" )
vsp.Update();



# fazer simulação em método de painéis 3D

print( "-> Begin TestVSPAeroComputeGeomPanel" ) 
# Print( string( "" ) );

# //open the file created in GenerateGeom
fname_vspaerotests = "TestVSPAero.vsp3"
vsp.ReadVSPFile( fname_vspaerotests )  #// Sets VSP3 file name

# //==== Analysis: VSPAero Compute Geometry ====//
analysis_name = "VSPAEROComputeGeometry"
print( analysis_name )

# // Set defaults
vsp.SetAnalysisInputDefaults( analysis_name )

# // Set to panel method
# analysis_method = vsp.GetIntAnalysisInput( analysis_name, "AnalysisMethod" )
analysis_method = [vsp.PANEL] 
vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method )

# // list inputs, type, and current values
vsp.PrintAnalysisInputs( analysis_name )
print( "" )

# // Execute
print( "Executing..." )
rid = vsp.ExecAnalysis( analysis_name )
print( "COMPLETE" )

# // Get & Display Results
vsp.PrintResults( rid )