# https://groups.google.com/g/openvsp/c/4Jh-T6AXWqE/m/AMKlv38HAQAJ

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:45:55 2021

@author: Hemant Sriram
"""

import SUAVE
from SUAVE.Core import Units, Data

import xlwt
import string
from openpyxl import Workbook

try:

    # import vsp_g as vsp
    import openvsp as vsp
    vsp.VSPCheckSetup()
    vsp.VSPRenew()
    
except ImportError:
    # This allows SUAVE to build without OpenVSP
    pass
import numpy as np

import time
import pdb

## @ingroup Input_Output-OpenVSP
def analysis(tag):
    t = time.time()   
    try:
        
        vsp.ClearVSPModel()
        
    except NameError:
        print('VSP import failed')
        
    vsp.VSPCheckSetup()
    vsp.VSPRenew()


    #open the file created in vsp_write

    vsp.ReadVSPFile(tag)
    
    #==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    compgeom_name = "VSPAEROComputeGeometry";

    # Set defaults
    vsp.SetAnalysisInputDefaults(compgeom_name);

    # list inputs, type, and current values
    vsp.PrintAnalysisInputs(compgeom_name);

    # Execute
    compgeom_resid=vsp.ExecAnalysis(compgeom_name);

    # Get & Display Results
    vsp.PrintResults( compgeom_resid );
    

    #==== Analysis: VSPAero Compute Geometry ====//

    analysis_name="VSPAeroSweep"

    #Set defaults

    vsp.SetAnalysisInputDefaults(analysis_name)

    #Change some input values
    #    Analysis method

    analysis_method = [vsp.VORTEX_LATTICE]


    vsp.SetIntAnalysisInput( analysis_name, "AnalysisMethod", analysis_method, 0 )
    
    #Reference geometry set
    
    geom_set=[0]
    vsp.SetIntAnalysisInput( analysis_name, "GeomSet", geom_set );
                           
    #Reference areas, lengths
    
    #Modified by Hemant Sriram 4 October 2021
    ids = vsp.FindGeoms()
    wing_id =[]
    
    for i in range(len(ids)):
        if (vsp.GetGeomTypeName(ids[i]) == 'Wing'):
            wing_id.append(ids[i])

    
    # wing_id=vsp.FindGeomsWithName("MainWing")
    
    sref=[float(0)]*1
    bref=[float(0)]*1
    cref=[float(0)]*1
    
    
    # pdb.set_trace()
    sref[0]=(vsp.GetParmVal(wing_id[0],"TotalArea", "WingGeom"))
    bref[0]=(vsp.GetParmVal(wing_id[0],"TotalSpan", "WingGeom"))
    cref[0]=(vsp.GetParmVal(wing_id[0],"TotalChord", "WingGeom"))
    
    ref_flag=[3]
    
    
    vsp.SetDoubleAnalysisInput( analysis_name, 'Sref', sref )
    vsp.SetDoubleAnalysisInput( analysis_name, 'bref', bref )
    vsp.SetDoubleAnalysisInput( analysis_name, 'cref', cref )
    vsp.SetIntAnalysisInput( analysis_name, "RefFlag", ref_flag )
    
                             
    #Freestream parameters
    #Alpha
    # pdb.set_trace()
    alpha_start=[float(1)]
    alpha_end=[float(1)]
    alpha_npts=[1]
    
    vsp.SetDoubleAnalysisInput( analysis_name, "AlphaStart", alpha_start )
    vsp.SetDoubleAnalysisInput( analysis_name, "AlphaEnd", alpha_end )
    vsp.SetIntAnalysisInput( analysis_name, "AlphaNpts", alpha_npts )
    # pdb.set_trace()
    #Beta
    beta_start=[float(0)]
    beta_end=[float(0)]
    beta_npts=[1]
    
    
    vsp.SetDoubleAnalysisInput( analysis_name, "BetaStart", beta_start )
    vsp.SetDoubleAnalysisInput( analysis_name, "BetaEnd", beta_end )
    vsp.SetIntAnalysisInput( analysis_name, "BetaNpts", beta_npts );
                           
    #Mach
    # pdb.set_trace()
    mach_start=[float(0.147)]
    mach_end=[float(0.147)]
    mach_npts=[1]
    
    vsp.SetDoubleAnalysisInput( analysis_name, "MachStart", mach_start )
    vsp.SetDoubleAnalysisInput( analysis_name, "MachEnd", mach_end )
    vsp.SetIntAnalysisInput( analysis_name, "MachNpts", mach_npts )
    
    # pdb.set_trace()                        
    vsp.Update()
              
    ##Case Setup
    # pdb.set_trace()
    wakeNumIter=[5]
    wakeSkipUntilIter=[6]
    batch_mode_flag=[1]
    
    vsp.SetIntAnalysisInput( analysis_name, "WakeNumIter", wakeNumIter )
    vsp.SetIntAnalysisInput( analysis_name, "WakeSkipUntilIter", wakeSkipUntilIter )
    vsp.SetIntAnalysisInput( analysis_name, "BatchModeFlag", batch_mode_flag )
    # pdb.set_trace()                        
    vsp.Update()


    #list inputs, type, and current values

    vsp.PrintAnalysisInputs(analysis_name)
    print("Analysis has be written to VSP")
    print("Starting / Executing Analysis")

    #Execute
    # pdb.set_trace()
    rid = vsp.ExecAnalysis(analysis_name)
    
    print("Analysis complete")

    #Get & Display Results

    vsp.PrintResults(rid)
                    
    #Write in CSV
    
    csvname='Simulation'
    
    vsp.WriteResultsCSVFile(rid,csvname+'.csv')

    # Check for errors

    errorMgr = vsp.ErrorMgrSingleton_getInstance()
    num_err = errorMgr.GetNumTotalErrors()
    for i in range(0, num_err):
        err = errorMgr.PopLastError()
        print("error = ", err.m_ErrorString)
        
    print ('FINISHED VSPAERO SIMULATION')
    
    data = []
    with open(csvname+'.csv') as f:
        for line in f:
            data.append([word for word in line.split(",") if word])
            
    wb = Workbook()
    sheet = wb.active
    for row_index in range(len(data)):
        for col_index, letter in zip(range(len(data[row_index])), string.ascii_uppercase):
            sheet[letter+str(row_index+1)]= data[row_index][col_index]

    wb.save(csvname+'.xlsx')
    elapsed = time.time() - t
    print("Elapsed Time is: {}".format(elapsed))
    return 