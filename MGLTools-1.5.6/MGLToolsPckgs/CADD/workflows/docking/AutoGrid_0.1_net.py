#!/bin/ksh ~/.mgltools/pythonsh
########################################################################
#
#    Vision Network - Python source code - file generated by vision
#    Wednesday 23 February 2011 18:15:15 
#    
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler, Michel Sanner and TSRI
#   
# revision: Guillaume Vareille
#  
#########################################################################
#
# $Header: /opt/cvs/CADD/workflows/docking/AutoGrid_0.1_net.py,v 1.1 2011/04/15 16:16:48 nadya Exp $
#
# $Id: AutoGrid_0.1_net.py,v 1.1 2011/04/15 16:16:48 nadya Exp $
#


if __name__=='__main__':
    from sys import argv
    if '--help' in argv or '-h' in argv or '-w' in argv: # run without Vision
        withoutVision = True
        from Vision.VPE import NoGuiExec
        ed = NoGuiExec()
        from NetworkEditor.net import Network
        import os
        masterNet = Network("process-"+str(os.getpid()))
        ed.addNetwork(masterNet)
    else: # run as a stand alone application while vision is hidden
        withoutVision = False
        from Vision import launchVisionToRunNetworkAsApplication, mainLoopVisionToRunNetworkAsApplication
	if '-noSplash' in argv:
	    splash = False
	else:
	    splash = True
        masterNet = launchVisionToRunNetworkAsApplication(splash=splash)
        import os
        masterNet.filename = os.path.abspath(__file__)
from traceback import print_exc
## loading libraries ##
from AutoDockTools.VisionInterface.Adt import Adt
from WebServices.VisionInterface.WSNodes import wslib
from Vision.StandardNodes import stdlib
try:
    masterNet
except (NameError, AttributeError): # we run the network outside Vision
    from NetworkEditor.net import Network
    masterNet = Network()

masterNet.getEditor().addLibraryInstance(Adt,"AutoDockTools.VisionInterface.Adt", "Adt")

masterNet.getEditor().addLibraryInstance(wslib,"WebServices.VisionInterface.WSNodes", "wslib")

masterNet.getEditor().addLibraryInstance(stdlib,"Vision.StandardNodes", "stdlib")

from WebServices.VisionInterface.WSNodes import addOpalServerAsCategory
try:
    addOpalServerAsCategory("http://ws.nbcr.net/opal2", replace=False)
except:
    pass
try:
    addOpalServerAsCategory("http://kryptonite.nbcr.net/opal2", replace=False)
except:
    pass
try:
    ## saving node ComputeGrids ##
    from Adt.Macro.ComputeGrids import ComputeGrids
    ComputeGrids_0 = ComputeGrids(constrkw={}, name='ComputeGrids', library=Adt)
    masterNet.addNode(ComputeGrids_0,391,290)
    apply(ComputeGrids_0.configure, (), {'paramPanelImmediate': 1, 'expanded': False})
    prepareGPF_kryptonite_nbcr_net_4 = ComputeGrids_0.macroNetwork.nodes[3]
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['singlelib'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['r_url'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['zpoints'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['filter_file_url'].widget.set(r"", run=False)
    apply(prepareGPF_kryptonite_nbcr_net_4.inputPortByName['lib'].widget.configure, (), {'choices': ('sample', 'NCIDS_SC', 'NCI_DS1', 'NCI_DS2', 'oldNCI', 'human_metabolome', 'chembridge_building_blocks', 'drugbank_nutraceutics', 'drugbank_smallmol', 'fda_approved')})
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['lib'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['ypoints'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['xcenter'].widget.set(r"auto", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['p'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['o'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['zcenter'].widget.set(r"auto", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['v'].widget.set(0, run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['userlib'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['xpoints'].widget.set(r"", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['localRun'].widget.set(0, run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['ycenter'].widget.set(r"auto", run=False)
    prepareGPF_kryptonite_nbcr_net_4.inputPortByName['execPath'].widget.set(r"", run=False)
    autogrid_kryptonite_nbcr_net_5 = ComputeGrids_0.macroNetwork.nodes[4]
    autogrid_kryptonite_nbcr_net_5.inputPortByName['infile_url'].widget.set(r"", run=False)
    autogrid_kryptonite_nbcr_net_5.inputPortByName['l'].widget.set(r"output.glg", run=False)
    autogrid_kryptonite_nbcr_net_5.inputPortByName['o'].widget.set(0, run=False)
    autogrid_kryptonite_nbcr_net_5.inputPortByName['p'].widget.set(r"", run=False)
    autogrid_kryptonite_nbcr_net_5.inputPortByName['localRun'].widget.set(0, run=False)
    autogrid_kryptonite_nbcr_net_5.inputPortByName['execPath'].widget.set(r"", run=False)
    GetURLFromList_8 = ComputeGrids_0.macroNetwork.nodes[7]
    GetURLFromList_8.inputPortByName['ext'].widget.set(r"gpf", run=False)

    ## saving connections for network ComputeGrids ##
    ComputeGrids_0.macroNetwork.freeze()
    ComputeGrids_0.macroNetwork.unfreeze()

    ## modifying MacroInputNode dynamic ports
    input_Ports_1 = ComputeGrids_0.macroNetwork.ipNode
    input_Ports_1.outputPorts[1].configure(name='GetComputeGridsInputs_ligands')
    input_Ports_1.outputPorts[2].configure(name='GetComputeGridsInputs_receptor_pdbqt')
    input_Ports_1.outputPorts[3].configure(name='GetComputeGridsInputs_gpf_obj')

    ## modifying MacroOutputNode dynamic ports
    output_Ports_2 = ComputeGrids_0.macroNetwork.opNode
    output_Ports_2.inputPorts[1].configure(singleConnection='auto')
    output_Ports_2.inputPorts[2].configure(singleConnection='auto')
    output_Ports_2.inputPorts[1].configure(name='MakeAutogridResultObj_autogrid_result_obj')
    output_Ports_2.inputPorts[2].configure(name='GetMainURLFromList_newurl')
    ComputeGrids_0.inputPorts[0].configure(name='GetComputeGridsInputs_ligands')
    ComputeGrids_0.inputPorts[0].configure(datatype='LigandDB')
    ComputeGrids_0.inputPorts[1].configure(name='GetComputeGridsInputs_receptor_pdbqt')
    ComputeGrids_0.inputPorts[1].configure(datatype='receptor_prepared')
    ComputeGrids_0.inputPorts[2].configure(name='GetComputeGridsInputs_gpf_obj')
    ComputeGrids_0.inputPorts[2].configure(datatype='gpf_template')
    ## configure MacroNode input ports
    ComputeGrids_0.outputPorts[0].configure(name='MakeAutogridResultObj_autogrid_result_obj')
    ComputeGrids_0.outputPorts[0].configure(datatype='autogrid_results')
    ComputeGrids_0.outputPorts[1].configure(name='GetMainURLFromList_newurl')
    ComputeGrids_0.outputPorts[1].configure(datatype='string')
    ## configure MacroNode output ports
    ComputeGrids_0.shrink()
    apply(ComputeGrids_0.configure, (), {'paramPanelImmediate': 1, 'expanded': False})
except:
    print "WARNING: failed to restore ComputeGrids named ComputeGrids in network masterNet"
    print_exc()
    ComputeGrids_0=None

try:
    ## saving node WebBrowser ##
    from WebServices.VisionInterface.WSNodes import WebBrowserNode
    WebBrowser_9 = WebBrowserNode(constrkw={}, name='WebBrowser', library=wslib)
    masterNet.addNode(WebBrowser_9,310,412)
    apply(WebBrowser_9.inputPortByName['url'].configure, (), {'defaultValue': None})
    WebBrowser_9.inputPortByName['url'].rebindWidget()
    WebBrowser_9.inputPortByName['url'].widget.set(r"", run=False)
    WebBrowser_9.inputPortByName['url'].unbindWidget()
    apply(WebBrowser_9.configure, (), {'paramPanelImmediate': 1})
except:
    print "WARNING: failed to restore WebBrowserNode named WebBrowser in network masterNet"
    print_exc()
    WebBrowser_9=None

try:
    ## saving node PublicServerLigandDB ##
    from Adt.Input.PublicServerLigandDB import PublicServerLigandDB
    PublicServerLigandDB_10 = PublicServerLigandDB(constrkw={}, name='PublicServerLigandDB', library=Adt)
    masterNet.addNode(PublicServerLigandDB_10,16,40)
    PublicServerLigandDB_10.inputPortByName['server_lib'].widget.set(r"sample", run=False)
except:
    print "WARNING: failed to restore PublicServerLigandDB named PublicServerLigandDB in network masterNet"
    print_exc()
    PublicServerLigandDB_10=None

try:
    ## saving node PrepareReceptor ##
    from Adt.Macro.PrepareReceptor import PrepareReceptor
    PrepareReceptor_11 = PrepareReceptor(constrkw={}, name='PrepareReceptor', library=Adt)
    masterNet.addNode(PrepareReceptor_11,408,153)
    apply(PrepareReceptor_11.configure, (), {'paramPanelImmediate': 1, 'expanded': False})
    Pdb2pqrWS_14 = PrepareReceptor_11.macroNetwork.nodes[2]
    Pdb2pqrOpalService_ws_nbcr_net_18 = Pdb2pqrWS_14.macroNetwork.nodes[3]
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['noopt'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['phi'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['psi'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['verbose'].widget.set(1, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['chain'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['nodebump'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['chi'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['ligand'].widget.set(r"", run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['hbond'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['with_ph'].widget.set(r"", run=False)
    apply(Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['forcefield'].widget.configure, (), {'choices': ('AMBER', 'CHARMM', 'PARSE', 'TYL06')})
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['forcefield'].widget.set(r"AMBER", run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['clean'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['inId'].widget.set(r"", run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['apbs_input'].widget.set(0, run=False)
    apply(Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['ffout'].widget.configure, (), {'choices': ('AMBER', 'CHARMM', 'PARSE', 'TYL06')})
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['ffout'].widget.set(r"", run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['localRun'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['rama'].widget.set(0, run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['execPath'].widget.set(r"", run=False)
    Pdb2pqrOpalService_ws_nbcr_net_18.inputPortByName['assign_only'].widget.set(0, run=False)
    GetURLFromList_19 = Pdb2pqrWS_14.macroNetwork.nodes[4]
    GetURLFromList_19.inputPortByName['ext'].widget.set(r"pqr", run=False)

    ## saving connections for network Pdb2pqrWS ##
    Pdb2pqrWS_14.macroNetwork.freeze()
    Pdb2pqrWS_14.macroNetwork.unfreeze()

    ## modifying MacroInputNode dynamic ports
    input_Ports_15 = Pdb2pqrWS_14.macroNetwork.ipNode
    input_Ports_15.outputPorts[1].configure(name='CheckFileFormat_value')

    ## modifying MacroOutputNode dynamic ports
    output_Ports_16 = Pdb2pqrWS_14.macroNetwork.opNode
    output_Ports_16.inputPorts[1].configure(singleConnection='auto')
    output_Ports_16.inputPorts[2].configure(singleConnection='auto')
    output_Ports_16.inputPorts[1].configure(name='UpdateReceptor_receptor_obj')
    output_Ports_16.inputPorts[2].configure(name='UpdateReceptor_pdb2pqr_result')
    Pdb2pqrWS_14.inputPorts[0].configure(name='CheckFileFormat_value')
    Pdb2pqrWS_14.inputPorts[0].configure(datatype='receptor')
    ## configure MacroNode input ports
    Pdb2pqrWS_14.outputPorts[0].configure(name='UpdateReceptor_receptor_obj')
    Pdb2pqrWS_14.outputPorts[0].configure(datatype='receptor')
    Pdb2pqrWS_14.outputPorts[1].configure(name='UpdateReceptor_pdb2pqr_result')
    Pdb2pqrWS_14.outputPorts[1].configure(datatype='string')
    ## configure MacroNode output ports
    Pdb2pqrWS_14.shrink()
    PrepareReceptorWS_21 = PrepareReceptor_11.macroNetwork.nodes[3]
    PrepareReceptorOpalService_ws_nbcr_net_25 = PrepareReceptorWS_21.macroNetwork.nodes[3]
    PrepareReceptorOpalService_ws_nbcr_net_25.inputPortByName['o'].widget.set(r"", run=False)
    PrepareReceptorOpalService_ws_nbcr_net_25.inputPortByName['v'].widget.set(0, run=False)
    PrepareReceptorOpalService_ws_nbcr_net_25.inputPortByName['localRun'].widget.set(0, run=False)
    PrepareReceptorOpalService_ws_nbcr_net_25.inputPortByName['execPath'].widget.set(r"", run=False)
    GetURLFromList_26 = PrepareReceptorWS_21.macroNetwork.nodes[4]
    GetURLFromList_26.inputPortByName['ext'].widget.set(r"pdbqt", run=False)
    DownloadToFile_27 = PrepareReceptorWS_21.macroNetwork.nodes[5]
    DownloadToFile_27.inputPortByName['overwrite'].widget.set(1, run=False)

    ## saving connections for network PrepareReceptorWS ##
    PrepareReceptorWS_21.macroNetwork.freeze()
    PrepareReceptorWS_21.macroNetwork.unfreeze()

    ## modifying MacroInputNode dynamic ports
    input_Ports_22 = PrepareReceptorWS_21.macroNetwork.ipNode
    input_Ports_22.outputPorts[1].configure(name='CheckFileFormat_value')
    input_Ports_22.outputPorts[2].configure(name='PrepareReceptorOpalService_ws_nbcr_net_C')

    ## modifying MacroOutputNode dynamic ports
    output_Ports_23 = PrepareReceptorWS_21.macroNetwork.opNode
    output_Ports_23.inputPorts[1].configure(singleConnection='auto')
    output_Ports_23.inputPorts[2].configure(singleConnection='auto')
    output_Ports_23.inputPorts[1].configure(name='UpdateReceptor_receptor_prepared_obj')
    output_Ports_23.inputPorts[2].configure(name='UpdateReceptor_receptor_result')
    PrepareReceptorWS_21.inputPorts[0].configure(name='CheckFileFormat_value')
    PrepareReceptorWS_21.inputPorts[0].configure(datatype='receptor')
    PrepareReceptorWS_21.inputPorts[1].configure(name='PrepareReceptorOpalService_ws_nbcr_net_C')
    PrepareReceptorWS_21.inputPorts[1].configure(datatype='boolean')
    ## configure MacroNode input ports
    PrepareReceptorWS_21.outputPorts[0].configure(name='UpdateReceptor_receptor_prepared_obj')
    PrepareReceptorWS_21.outputPorts[0].configure(datatype='receptor_prepared')
    PrepareReceptorWS_21.outputPorts[1].configure(name='UpdateReceptor_receptor_result')
    PrepareReceptorWS_21.outputPorts[1].configure(datatype='string')
    ## configure MacroNode output ports
    PrepareReceptorWS_21.shrink()

    ## saving connections for network PrepareReceptor ##
    PrepareReceptor_11.macroNetwork.freeze()
    PrepareReceptor_11.macroNetwork.unfreeze()

    ## modifying MacroInputNode dynamic ports
    input_Ports_12 = PrepareReceptor_11.macroNetwork.ipNode
    input_Ports_12.outputPorts[1].configure(name='Pdb2pqrWS_CheckFileFormat_value')
    input_Ports_12.outputPorts[2].configure(name='PrepareReceptorWS_PrepareReceptorOpalService_ws_nbcr_net_C')

    ## modifying MacroOutputNode dynamic ports
    output_Ports_13 = PrepareReceptor_11.macroNetwork.opNode
    output_Ports_13.inputPorts[1].configure(singleConnection='auto')
    output_Ports_13.inputPorts[2].configure(singleConnection='auto')
    output_Ports_13.inputPorts[1].configure(name='PrepareReceptorWS_UpdateReceptor_receptor_prepared_obj')
    output_Ports_13.inputPorts[2].configure(name='PrepareReceptorWS_UpdateReceptor_receptor_result')
    PrepareReceptor_11.inputPorts[0].configure(name='Pdb2pqrWS_CheckFileFormat_value')
    PrepareReceptor_11.inputPorts[0].configure(datatype='receptor')
    PrepareReceptor_11.inputPorts[1].configure(name='PrepareReceptorWS_PrepareReceptorOpalService_ws_nbcr_net_C')
    PrepareReceptor_11.inputPorts[1].configure(datatype='boolean')
    ## configure MacroNode input ports
    PrepareReceptor_11.outputPorts[0].configure(name='PrepareReceptorWS_UpdateReceptor_receptor_prepared_obj')
    PrepareReceptor_11.outputPorts[0].configure(datatype='receptor_prepared')
    PrepareReceptor_11.outputPorts[1].configure(name='PrepareReceptorWS_UpdateReceptor_receptor_result')
    PrepareReceptor_11.outputPorts[1].configure(datatype='string')
    ## configure MacroNode output ports
    PrepareReceptor_11.shrink()
    apply(PrepareReceptor_11.configure, (), {'paramPanelImmediate': 1, 'expanded': False})
except:
    print "WARNING: failed to restore PrepareReceptor named PrepareReceptor in network masterNet"
    print_exc()
    PrepareReceptor_11=None

try:
    ## saving node StructureBrowser ##
    from Adt.Input.StructureBrowser import StructureBrowser
    StructureBrowser_29 = StructureBrowser(constrkw={}, name='StructureBrowser', library=Adt)
    masterNet.addNode(StructureBrowser_29,408,35)
    StructureBrowser_29.inputPortByName['receptor_file'].widget.set(r"AutoGrid_0.1_input/2HTY_A.pdbqt", run=False)
    apply(StructureBrowser_29.configure, (), {'paramPanelImmediate': 1})
except:
    print "WARNING: failed to restore StructureBrowser named StructureBrowser in network masterNet"
    print_exc()
    StructureBrowser_29=None

try:
    ## saving node GPFTemplateBrowser ##
    from Adt.Input.GPFTemplateBrowser import GPFTemplateBrowser
    GPFTemplateBrowser_30 = GPFTemplateBrowser(constrkw={}, name='GPFTemplateBrowser', library=Adt)
    masterNet.addNode(GPFTemplateBrowser_30,826,34)
    GPFTemplateBrowser_30.inputPortByName['gpf_template_file'].widget.set(r"AutoGrid_0.1_input/2HTY_A.gpf", run=False)
    apply(GPFTemplateBrowser_30.configure, (), {'paramPanelImmediate': 1})
except:
    print "WARNING: failed to restore GPFTemplateBrowser named GPFTemplateBrowser in network masterNet"
    print_exc()
    GPFTemplateBrowser_30=None

try:
    ## saving node DownloadSaveDir ##
    from WebServices.VisionInterface.WSNodes import DownloadSaveDirNode
    DownloadSaveDir_31 = DownloadSaveDirNode(constrkw={}, name='DownloadSaveDir', library=wslib)
    masterNet.addNode(DownloadSaveDir_31,502,397)
    apply(DownloadSaveDir_31.inputPortByName['url'].configure, (), {'defaultValue': None})
    DownloadSaveDir_31.inputPortByName['url'].rebindWidget()
    DownloadSaveDir_31.inputPortByName['url'].widget.set(r"", run=False)
    DownloadSaveDir_31.inputPortByName['url'].unbindWidget()
    apply(DownloadSaveDir_31.configure, (), {'paramPanelImmediate': 1})
except:
    print "WARNING: failed to restore DownloadSaveDirNode named DownloadSaveDir in network masterNet"
    print_exc()
    DownloadSaveDir_31=None

try:
    ## saving node Output Directory Browser ##
    from Vision.StandardNodes import DirBrowserNE
    Output_Directory_Browser_32 = DirBrowserNE(constrkw={}, name='Output Directory Browser', library=stdlib)
    masterNet.addNode(Output_Directory_Browser_32,905,178)
    Output_Directory_Browser_32.inputPortByName['directory'].widget.set(r"AutoGrid_0.1_output", run=False)
    apply(Output_Directory_Browser_32.configure, (), {'paramPanelImmediate': 1})
except:
    print "WARNING: failed to restore DirBrowserNE named Output Directory Browser in network masterNet"
    print_exc()
    Output_Directory_Browser_32=None

#masterNet.run()
masterNet.freeze()

## saving connections for network AutoGrid_0.1 ##
if PublicServerLigandDB_10 is not None and ComputeGrids_0 is not None:
    try:
        masterNet.connectNodes(
            PublicServerLigandDB_10, ComputeGrids_0, "ligDB", "GetComputeGridsInputs_ligands", blocking=True
            , splitratio=[0.71040241021275885, 0.32432210892569779])
    except:
        print "WARNING: failed to restore connection between PublicServerLigandDB_10 and ComputeGrids_0 in network masterNet"
if PrepareReceptor_11 is not None and ComputeGrids_0 is not None:
    try:
        masterNet.connectNodes(
            PrepareReceptor_11, ComputeGrids_0, "PrepareReceptorWS_UpdateReceptor_receptor_prepared_obj", "GetComputeGridsInputs_receptor_pdbqt", blocking=True
            , splitratio=[0.37771324748037094, 0.29247414694141388])
    except:
        print "WARNING: failed to restore connection between PrepareReceptor_11 and ComputeGrids_0 in network masterNet"
if StructureBrowser_29 is not None and PrepareReceptor_11 is not None:
    try:
        masterNet.connectNodes(
            StructureBrowser_29, PrepareReceptor_11, "receptor_obj", "Pdb2pqrWS_CheckFileFormat_value", blocking=True
            , splitratio=[0.62789739082139784, 0.69687814536090076])
    except:
        print "WARNING: failed to restore connection between StructureBrowser_29 and PrepareReceptor_11 in network masterNet"
if GPFTemplateBrowser_30 is not None and ComputeGrids_0 is not None:
    try:
        masterNet.connectNodes(
            GPFTemplateBrowser_30, ComputeGrids_0, "gpf_template", "GetComputeGridsInputs_gpf_obj", blocking=True
            , splitratio=[0.39300490327970311, 0.4601676294958682])
    except:
        print "WARNING: failed to restore connection between GPFTemplateBrowser_30 and ComputeGrids_0 in network masterNet"
if ComputeGrids_0 is not None and WebBrowser_9 is not None:
    try:
        masterNet.connectNodes(
            ComputeGrids_0, WebBrowser_9, "GetMainURLFromList_newurl", "url", blocking=True
            , splitratio=[0.30896252552604309, 0.21760077157963498])
    except:
        print "WARNING: failed to restore connection between ComputeGrids_0 and WebBrowser_9 in network masterNet"
if ComputeGrids_0 is not None and DownloadSaveDir_31 is not None:
    try:
        masterNet.connectNodes(
            ComputeGrids_0, DownloadSaveDir_31, "GetMainURLFromList_newurl", "url", blocking=True
            , splitratio=[0.66446304698655212, 0.30837982537508668])
    except:
        print "WARNING: failed to restore connection between ComputeGrids_0 and DownloadSaveDir_31 in network masterNet"
if Output_Directory_Browser_32 is not None and DownloadSaveDir_31 is not None:
    try:
        masterNet.connectNodes(
            Output_Directory_Browser_32, DownloadSaveDir_31, "directory", "path", blocking=True
            , splitratio=[0.34510585020439977, 0.73683619505012499])
    except:
        print "WARNING: failed to restore connection between Output_Directory_Browser_32 and DownloadSaveDir_31 in network masterNet"
masterNet.runOnNewData.value = False

if __name__=='__main__':
    from sys import argv
    lNodePortValues = []
    if (len(argv) > 1) and argv[1].startswith('-'):
        lArgIndex = 2
    else:
        lArgIndex = 1
    while lArgIndex < len(argv) and argv[lArgIndex][-3:]!='.py':
        lNodePortValues.append(argv[lArgIndex])
        lArgIndex += 1
    masterNet.setNodePortValues(lNodePortValues)
    if '--help' in argv or '-h' in argv: # show help
        masterNet.helpForNetworkAsApplication()
    elif '-w' in argv: # run without Vision and exit
         # create communicator
        from NetworkEditor.net import Communicator
        masterNet.communicator = Communicator(masterNet)
        print 'Communicator listening on port:', masterNet.communicator.port

        import socket
        f = open(argv[0]+'.sock', 'w')
        f.write("%s %i"%(socket.gethostbyname(socket.gethostname()),
                         masterNet.communicator.port))
        f.close()

        masterNet.run()

    else: # stand alone application while vision is hidden
        if '-e' in argv: # run and exit
            masterNet.run()
        elif '-r' in argv or len(masterNet.userPanels) == 0: # no user panel => run
            masterNet.run()
            mainLoopVisionToRunNetworkAsApplication(masterNet.editor)
        else: # user panel
            mainLoopVisionToRunNetworkAsApplication(masterNet.editor)

