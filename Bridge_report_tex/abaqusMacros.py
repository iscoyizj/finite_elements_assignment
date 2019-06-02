# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def Macro1():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.ModelFromInputFile(name='Bridge-1', 
        inputFileName='D:/Abacus/FEM Course/Bridge/Bridge-1.inp')
    a = mdb.models['Bridge-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    mdb.Job(name='Job-1', model='Bridge-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)

def Getmaxdisplacement():
    import math
	# Get the current working directory
    cwd=os.getcwd()
    # Open the .odb
    o=session.openOdb(name='{0}//Job-{1}.odb'.format(cwd,1));
    RS=o.steps['Step-1'].frames[-1].fieldOutputs['U'].values;
    maxdisp=0
    for i in RS :
        l=math.sqrt(i.data[0]*i.data[0]+i.data[1]*i.data[1]+i.data[2]*i.data[2])
        if l>maxdisp:
            maxdisp=l
    return maxdisp
	

def Getp():
    import math
    # Get the current working directory
    cwd=os.getcwd()
    # Open the .odb
    o=session.openOdb(name='{0}//Job-{1}.odb'.format(cwd,1))	
    RF=o.steps['Step-1'].frames[-1].fieldOutputs['S'].values
    with open( 'D:\Abacus\FEM Course\Bridge\Scomponent.txt','w') as f1:
     for v in RF:
        f1.writelines( str( v.baseElementType) + ' ')
        f1.writelines( str( v.elementLabel) + ' ')
        f1.writelines( str( v.mises) + ' \n')
     else:
        print 'ok'
        f1.close
    with open( 'D:\Abacus\FEM Course\Bridge\elmementofset.txt','w') as f2:
     for v in RF:
        elementnode = o.rootAssembly.instances.elements[v.elementLabel-1].connectivity
        f1.writelines( str( v.baseElementType) + ' ')
        f2.writelines( str( v.elementLabel) + ' ')
        f2.writelines( str( elementnode) + ' \n')
     else:
        f2.close
        print 'ok2'
    nd = mdb.models['Brige-1'].parts.nodes
    with open( 'D:\Abacus\FEM Course\Bridge\nodeofset.txt','w') as f3:
     v = 1
     while v <= len(nd) :
        nodescoors = mdb.models['Model-1'].parts.sets.nodes[v - 1].coordinates
        nodeslabel = mdb.models['Model-1'].parts.sets.nodes[v - 1].label
        f3.writelines( str( nodeslabel) + ' ')
        f3.writelines( str( nodescoors) + ' \n')
        v = v + 1
     else:
        f3.close
        print 'ok3'	
	
	
		
def Ex4Post():
    ALLIE2=[]
    E=Getmaxdisplacement();   
    ALLIE2.append(E);
    Getp();		
    return ALLIE2