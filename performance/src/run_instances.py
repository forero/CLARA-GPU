import os.path
from math import *
import subprocess as sp
import datetime

def write_report(report_name, command="ls",
                 TotalPhotons=1e4, BusyPhotons=1e5, BlockSize=128):

    a = sp.Popen("time " +command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    lines = a.stderr.readlines()
    timeline = lines[0].rstrip('\n')
    out = open(report_name, "a")
    out.write("%d %d %d %s %s\n"%(TotalPhotons, BusyPhotons, BlockSize, timeline, command))
    out.close()

def write_config(filename, InputDir="./", OutputDir="./", BaseName="toto",
                 TotalPhotons=1e4, BusyPhotons=1e5, BlockSize=128):
    out = open(filename, "w")
    out.write("BlockSize                          %d\n"%(BlockSize))
    out.write("BusyPhotons                        %d\n"%(BusyPhotons))
    out.write("TotalPhotons                       %d\n"%(TotalPhotons))
    out.write("InputDir                           %s\n"%(InputDir))
    out.write("CubeName                           grid_object_100\n")
    out.write("OutputDir                          %s\n"%(OutputDir))
    out.write("OutputFile                         %s\n"%(BaseName))
    out.write("NeufeldSlab                        0\n")
    out.write("NeufeldCube                        0\n")
    out.write("ExpandingSphere                    1\n")
    out.write("TestParallelVel                    0\n")
    out.write("TestParallelVelFast                0\n")
    out.write("TestFirstScatter                   0\n")
    out.write("TestRND                            0\n")
    out.write("TestPerpVel                        0\n")
    out.write("SimulationCube                     0\n")
    out.write("HomogeneousInit                    1\n")
    out.write("Test_a                             0.00014\n")
    out.write("Test_x                             2\n")
    out.write("OutputTestFile                     test_b\n")
    out.write("Temperature                        10000\n")
    out.write("Tau                                1000\n")
    out.write("NumberDensityHI                    100000\n")
    out.write("VmaxSphere                         0\n")
    out.write("GrainSize                          1e-06\n")
    out.write("TauDust                            0\n")
    out.write("DustAbsorptionProb                 0.5\n")
    out.write("InputFrequency                     0\n")
    out.write("LuminosityPerPackage               1\n")
    out.write("EffectiveEta                       0.71\n")
    out.write("EffectiveDust                      0.001\n")
    out.write("ThresholdCube                      10000\n")
    out.write("UseDust                            0\n")
    out.write("UseVelocities                      0\n")
    out.write("UseAtomSpeedUp                     0\n")
    out.write("OutputInitList                     1\n")
    out.write("OutputFinalList                    1\n")
    out.write("OutputBinary                       0\n")
    out.write("UnitMass_in_g                      1.989e+43\n")
    out.write("UnitVelocity_in_cm_per_s           100000\n")
    out.write("UnitLength_in_cm                   3.08568e+21\n")
    out.write("UnitLymanLuminosity                1e+41\n")
    out.write("RandomSeed                         14923\n")
    out.close()


TotalPhotonsList = [1e3, 1e4, 1e5]
BusyPhotonsList = [1,10,100]
BlockSizeList = [128, 256, 512, 1024]
RUNID=range(5)

OutputPath="/home/forero/Dropbox/CLARA-GPU/performance/output/"
ConfigPath="/home/forero/Dropbox/CLARA-GPU/performance/config/"
ReportPath="/home/forero/Dropbox/CLARA-GPU/performance/report/"
ExecPath="/home/forero/Dropbox/CLARA-GPU/src_kernel_3/"
ExecName="./clara"



now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d_%H-%M")
now = str(now)

ReportName="clara_GPU_report_"+now+".dat"


for TotalPhoton in TotalPhotonsList:
    for BusyPhoton in BusyPhotonsList:
        for BlockSize in BlockSizeList:
            for i in RUNID:
                RootName="RUN_%d_BusyPh_%d_TotalPh_%d_BlockSize_%d"%\
                (i, log10(BusyPhoton*TotalPhoton), log10(TotalPhoton), BlockSize)

                ParamName=ConfigPath+RootName+".input"
                OutputFile=OutputPath+RootName+"_out.ascii"
                InputFile=OutputPath+RootName+"_out.ascii"

                print i, BusyPhoton*TotalPhoton, TotalPhoton, BlockSize

                if (not os.path.exists(OutputFile)):
                    write_config(ParamName, InputDir=ConfigPath, 
                                 OutputDir=OutputPath, BaseName=RootName,
                                 TotalPhotons=TotalPhoton, 
                                 BusyPhotons=BusyPhoton*TotalPhoton, 
                                 BlockSize=BlockSize)

                    run = ExecPath+ExecName+" "+ParamName
                    print run
                    write_report(ReportPath+ReportName, command=run,
                                 TotalPhotons=TotalPhoton, 
                                 BusyPhotons=BusyPhoton*TotalPhoton, 
                                 BlockSize=BlockSize)


