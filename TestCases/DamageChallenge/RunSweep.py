import json
import numpy
import os
import subprocess

DoSims = True

def PrepareInOutFiles(SimName, InputJSON):
    folderName = "Results/"+SimName
    fileName = folderName + "/inputs.json"

    InputJSON["Restart"]["Saveloc"] = folderName + "/Results/Backup"
    InputJSON["Outputs"]["SaveFolder"] = folderName + "/Results"

    #json_formatted_str = json.dumps(RefData, indent=2)
    #print(json_formatted_str)

    if (os.path.isdir(folderName)==False):
        os.makedirs(folderName)

    if (os.path.isfile(fileName)==True):
        os.remove(fileName)

    with open(fileName, "w") as outfile:
        json_object = json.dumps(InputJSON, indent=4)
        outfile.write(json_object)
    
    if (os.path.isdir(folderName+"/Results")==False):
        os.makedirs(folderName+"/Results")
    

def RunSim(SimName):
    subprocess.run(['mpiexec',"-np", "8", "./IceCode","Results/"+SimName+"/inputs.json"])

def PrepareBatch(SeriesName, SimName):
    with open(SeriesName, "a") as outfile:
        outfile.write(SimName+"\n")


###### MAIN CODE
with open("DamageChallenge.json", "r") as read_file:
    RefData = json.load(read_file)

SeriesName = "ToRun.txt"
if (os.path.isfile(SeriesName)==True):
    os.remove(SeriesName)

for split in ["SpecStrains","VolStrains"]:
#    for mesh in ["quadratic_025mm", "cubic_025mm", "triangleLinear_025mm", "triangleQuadratic_025mm", "quadratic_0125mm", "cubic_0125mm", "triangleLinear_0125mm", "triangleQuadratic_0125mm"]:
    for mesh in ["quadratic_0125mm", "cubic_0125mm", "triangleLinear_0125mm", "triangleQuadratic_0125mm"]:
        for dist in ["Linear","Quadratic","HO_Linear","HO_Quadratic"]:
            for case in ["I","II"]:
                if ((dist=="HO_Linear" or dist=="HO_Quadratic") and (mesh=="triangleLinear_025mm" or mesh=="triangleQuadratic_025mm" or mesh=="triangleLinear_0125mm" or mesh=="triangleQuadratic_0125mm")):
                    pass
                else:
                    SimName = case+"_"+split+"_"+dist+"_"+mesh
                    print(SimName)

                    CurData = RefData
                    CurData["properties"]["PhaseField"]["Split"] = split
                    CurData["properties"]["PhaseField"]["DistributionFunction"] = dist
                    CurData["mesh"]["file"] = mesh+".h5"
                    if (case=="I"):
                        CurData["nonlinSolver"]["Initialization"]["Notch_x"] = 38.1e-3
                    else:
                        CurData["nonlinSolver"]["Initialization"]["Notch_x"] = 47.63e-3

                    PrepareInOutFiles(SimName, CurData)

                    if DoSims:
                        RunSim(SimName)
                    else:
                        PrepareBatch(SeriesName, SimName)
