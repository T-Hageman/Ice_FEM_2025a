#include "Plot_BoundaryFill.h"
#ifdef LOCALENVIRONMENT

#include <iostream>
#include <set>

#include "../physics.h"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellType.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkChartXYZ.h>
#include <vtkPlotSurface.h>
#include <vtkTable.h>
#include <vtkChartXYZ.h>
#include <vtkChartXY.h>
#include <vtkContextMouseEvent.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkPlotSurface.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vtkVersion.h>
#include <vtkCubeAxesActor.h>
#include <vtkColorSeries.h>
#include <vtkVariantArray.h>
#include <vtkAxis.h>
#include <vtkChartLegend.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>


Plot_BoundaryFill::Plot_BoundaryFill(std::string name, inputData &inputs, Physics *phys): Vis_Plot(name, inputs, phys) {
    std::string ElemGroupIndex_Name; inputs.GetRequired(ElemGroupIndex_Name, {"Visualisation",Name,"ElementGroup"});
    ElemGroupIndex_Data = physics->mesh->GetElementGroupIdx(ElemGroupIndex_Name);

    Deformations = false;
    inputs.GetOptional(Deformations, {"Visualisation",Name,"Deformed"});

    if(Deformations){
        std::string UGroup; inputs.GetRequired(UGroup, {"Visualisation","ElementGroup_u"});
        ElemGroupIndex_u = physics->mesh->GetElementGroupIdx(UGroup);

        DeformationScale = 1000.0;
        inputs.GetOptional(DeformationScale, {"Visualisation",Name,"DeformedScale"});
        dofTypes_u.resize(physics->mesh->dim);
        for (size_t i = 0; i < physics->mesh->dim; i++){
            physics->dofspace->getDofTypesSteps(DofNames_u[i], dofTypes_u[i], dofStep_u);
        }
    }

    inputs.GetRequired(DataName, {"Visualisation",Name,"Dof"});
    size_t UNUSED;
    phys->dofspace->getDofTypesSteps(DataName, UNUSED, DofStep_Data);

    inputs.GetRequired(GrowDir, {"Visualisation",Name,"GrowDir"});
    DataScale = 1.0;
    inputs.GetOptional(DataScale, {"Visualisation",Name,"DataScale"});

    LogScale = false;
    inputs.GetOptional(LogScale, {"Visualisation",Name,"LogScale"});

    if (rank == 0){
        // Set the background color.
        colors->SetColor("BkgColor", colors->GetColor4d("Snow").GetData());
        renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
        renderWindow->AddRenderer(renderer);
        renderer->ResetCamera();
    }
}

Plot_BoundaryFill::~Plot_BoundaryFill(){

}

void Plot_BoundaryFill::Plot(){
    std::vector<std::vector<double>> XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords;

    size_t nNodes_per_Elem = physics->mesh->GetExportMeshCoords(ElemGroupIndex_Data, XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords);
    size_t nElems = XCoords.size();

    std::vector<std::vector<double>> NodeData, deformX, deformY, NotNeeded;
    NodeData.resize(nElems); 
    for (size_t i = 0; i < nElems; i++) NodeData[i].resize(nNodes_per_Elem);

    physics->GetNodalDataToSave(DataName, ElemGroupIndex_Data, NodeData);
    if (LogScale){
        GetLog10(NodeData);
    }

    if (Deformations){
        deformX.resize(nElems); 
        for (size_t i = 0; i < nElems; i++) deformX[i].resize(nNodes_per_Elem);

        deformY.resize(nElems); 
        for (size_t i = 0; i < nElems; i++) deformY[i].resize(nNodes_per_Elem);

        physics->GetNodalDataToSave(DofNames_u[0], ElemGroupIndex_u, deformX);
        physics->GetNodalDataToSave(DofNames_u[1], ElemGroupIndex_u, deformY);
    }

    std::vector<double> DRange(2); DRange[0] = min(NodeData); DRange[1] = max(NodeData);
    if(size>1){
        MPI_Allreduce(MPI_IN_PLACE, &DRange[0], 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &DRange[1], 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    }

    vtkNew<vtkLookupTable> colorLookupTable; 
    colorLookupTable->SetHueRange(0.66667, 0);
    double offset = 0.0;
    if (DRange[0] == DRange[1]){
        offset = 1e-6;
    }
    colorLookupTable->SetTableRange(DRange[0]-offset, DRange[1]+offset);
    colorLookupTable->Build();
    
    MPI_Barrier(PETSC_COMM_WORLD);

    if (nElems>0){
        nNodes_per_Elem = XCoords[0].size();
    } else {
        nNodes_per_Elem = 0;
    }
    if (size>1){
        MPI_Allreduce(MPI_IN_PLACE, &nNodes_per_Elem, 1, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    }


    std::vector<size_t> order=physics->mesh->ElementGroups[ElemGroupIndex_Data].BaseElem->PlottingOrder;
    size_t nPoints = 2*order.size()*nElems;

    std::vector<double> x0(nPoints), x1(nPoints), x2(nPoints), c0(nPoints);
    std::vector<size_t> connects_i(nPoints);
    double c;

    size_t ptCounter = 0;
    for (size_t i = 0; i < nElems; i++){
        for (size_t j = 0; j < order.size(); j++){
            size_t pt = order[j];
            if (Deformations){
                x0[ptCounter] = XCoords[i][pt]+DeformationScale*deformX[i][pt];
                x1[ptCounter] = YCoords[i][pt]+DeformationScale*deformY[i][pt];
            } else {
                x0[ptCounter] = XCoords[i][pt];
                x1[ptCounter] = YCoords[i][pt];
            }
            x2[ptCounter] = 0.0;
            c0[ptCounter] = NodeData[i][pt];

            connects_i[ptCounter] = i;

            ptCounter += 1;
        }
        for (int j = order.size()-1; j >= 0; j--){
            size_t pt = order[j];

            if (GrowDir=="+y"){
                if (Deformations){
                    x0[ptCounter] = XCoords[i][pt]+DeformationScale*deformX[i][pt];
                    x1[ptCounter] = YCoords[i][pt]+DeformationScale*deformY[i][pt]+NodeData[i][pt]*DataScale;
                } else {
                    x0[ptCounter] = XCoords[i][pt];
                    x1[ptCounter] = YCoords[i][pt]+NodeData[i][pt]*DataScale;
                }
            } else if (GrowDir=="-y"){
                if (Deformations){
                    x0[ptCounter] = XCoords[i][pt]+DeformationScale*deformX[i][pt];
                    x1[ptCounter] = YCoords[i][pt]+DeformationScale*deformY[i][pt]-NodeData[i][pt]*DataScale;
                } else {
                    x0[ptCounter] = XCoords[i][pt];
                    x1[ptCounter] = YCoords[i][pt]-NodeData[i][pt]*DataScale;
                }
            } else if (GrowDir=="+x"){
                if (Deformations){
                    x0[ptCounter] = XCoords[i][pt]+DeformationScale*deformX[i][pt]+NodeData[i][pt]*DataScale;
                    x1[ptCounter] = YCoords[i][pt]+DeformationScale*deformY[i][pt];
                } else {
                    x0[ptCounter] = XCoords[i][pt]+NodeData[i][pt]*DataScale;
                    x1[ptCounter] = YCoords[i][pt];
                }
            } else {
                if (Deformations){
                    x0[ptCounter] = XCoords[i][pt]+DeformationScale*deformX[i][pt]-NodeData[i][pt]*DataScale;
                    x1[ptCounter] = YCoords[i][pt]+DeformationScale*deformY[i][pt];
                } else {
                    x0[ptCounter] = XCoords[i][pt]-NodeData[i][pt]*DataScale;
                    x1[ptCounter] = YCoords[i][pt];
                }
            }
            
            x2[ptCounter] = 0.0;
            c0[ptCounter] = NodeData[i][pt];

            connects_i[ptCounter] = i;

            ptCounter += 1;
        }
    }

    std::vector<size_t> NPoints(size), NElems(size);
    MPI_Gather(&nPoints, 1, my_MPI_SIZE_T, NPoints.data(), 1, my_MPI_SIZE_T, 0, PETSC_COMM_WORLD);
    MPI_Gather(&nElems, 1, my_MPI_SIZE_T, NElems.data(), 1, my_MPI_SIZE_T, 0, PETSC_COMM_WORLD);

    vtkIdType connects[2*order.size()]; 
    std::vector<MPI_Request> MPI_Reqs(5);
    std::vector<MPI_Status> MPI_Stats(5);

    if(rank==0){
        vtkNew<vtkUnsignedCharArray> colorData;
        colorData->SetNumberOfComponents(3);
        colorData->Allocate(std::reduce(NPoints.begin(),NPoints.end()));
        colorData->SetName("Colors");

        points->Allocate(std::reduce(NPoints.begin(),NPoints.end()));
        ugrid->Allocate(std::reduce(NElems.begin(),NElems.end()));

        double x[3];
        size_t iPrevious;
        size_t totalPoints = 0;

        for (size_t core = 0; core < size; core++){
            if (core!=0){
                x0.resize(NPoints[core]);
                x1.resize(NPoints[core]);
                x2.resize(NPoints[core]);
                c0.resize(NPoints[core]);

                connects_i.resize(NPoints[core]);

                if (NPoints[core]>0){
                    MPI_Irecv(x0.data(), NPoints[core], MPI_DOUBLE, core, core*10+0, PETSC_COMM_WORLD, &MPI_Reqs[0]);
                    MPI_Irecv(x1.data(), NPoints[core], MPI_DOUBLE, core, core*10+1, PETSC_COMM_WORLD, &MPI_Reqs[1]);
                    MPI_Irecv(x2.data(), NPoints[core], MPI_DOUBLE, core, core*10+2, PETSC_COMM_WORLD, &MPI_Reqs[2]);
                    MPI_Irecv(c0.data(), NPoints[core], MPI_DOUBLE, core, core*10+3, PETSC_COMM_WORLD, &MPI_Reqs[3]);

                    MPI_Irecv(connects_i.data(), NPoints[core], my_MPI_SIZE_T, core, core*10+4, PETSC_COMM_WORLD, &MPI_Reqs[4]);

                    MPI_Waitall(5, MPI_Reqs.data(), MPI_Stats.data());
                }
            }
            
            if (NPoints[core]>0){
                size_t idx_con = 0;
                iPrevious = connects_i[0];

                for (size_t i = 0; i < NPoints[core]; i++){
                    x[0] = x0[i];
                    x[1] = x1[i];
                    x[2] = x2[i];
                    c = c0[i];

                    points->InsertPoint(totalPoints, x);

                    if (connects_i[i]!=iPrevious){
                        ugrid->InsertNextCell(VTK_POLYGON, 2*order.size(), connects);
                        idx_con = 0;
                    } 
                    connects[idx_con] = totalPoints;
                    idx_con += 1;
                    iPrevious = connects_i[i];

                    double dcolor[3];
                    unsigned char ccolor[3];
                    colorLookupTable->GetColor(c, dcolor);
                    for (size_t k = 0; k < 3; k++){
                        ccolor[k] = static_cast<unsigned char>(255.0 * dcolor[k]);
                    }
                    colorData->InsertNextTypedTuple(ccolor);

                    totalPoints += 1;
                }
                //ugrid->InsertNextCell(VTK_POLYGON, 2*order.size(), connects);
            }
            
        }

        ugrid->SetPoints(points);
        ugrid->GetPointData()->SetScalars(colorData);
    } else {
        if (nPoints>0){
            MPI_Isend(x0.data(), nPoints, MPI_DOUBLE, 0, rank*10+0, PETSC_COMM_WORLD, &MPI_Reqs[0]);
            MPI_Isend(x1.data(), nPoints, MPI_DOUBLE, 0, rank*10+1, PETSC_COMM_WORLD, &MPI_Reqs[1]);
            MPI_Isend(x2.data(), nPoints, MPI_DOUBLE, 0, rank*10+2, PETSC_COMM_WORLD, &MPI_Reqs[2]);
            MPI_Isend(c0.data(), nPoints, MPI_DOUBLE, 0, rank*10+3, PETSC_COMM_WORLD, &MPI_Reqs[3]);

            MPI_Isend(connects_i.data(), nPoints, my_MPI_SIZE_T, 0, rank*10+4, PETSC_COMM_WORLD, &MPI_Reqs[4]);
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);


    if (rank==0){
        ColourBar->SetLookupTable(colorLookupTable);
        if (DoOnce){
            ugridMapper->SetInputData(ugrid);

            ugridActor->SetMapper(ugridMapper);
            ugridActor->GetProperty()->EdgeVisibilityOff();

            renderer->AddActor(ugridActor);

            std::vector<double> XRange(2);
            std::vector<double> YRange(2);
            XRange[0] = ugridMapper->GetBounds()[0];
            XRange[1] = ugridMapper->GetBounds()[1];
            YRange[0] = ugridMapper->GetBounds()[2];
            YRange[1] = ugridMapper->GetBounds()[3];

            renderer->GetActiveCamera()->ParallelProjectionOff();
            renderer->GetActiveCamera()->SetPosition(0.5*(XRange[0]+XRange[1]), 0.5*(YRange[0]+YRange[1]), 3*std::max(XRange[1]-XRange[0], YRange[1]-YRange[0]));
            renderer->GetActiveCamera()->SetFocalPoint(0.5*(XRange[0]+XRange[1]), 0.5*(YRange[0]+YRange[1]), 0);

            AddGrid(ugridMapper->GetBounds());
            AddColourBar();

            renderWindow->SetWindowName(DataName.c_str());
            // Set the icon for the window
            auto icon = LoadIconFromEmbeddedPNG();
            renderWindow->SetIcon(icon);

            textActor->SetPosition2(10, 10);
            textActor->GetTextProperty()->SetFontSize(18);
            textActor->GetTextProperty()->SetColor(colors->GetColor3d("Black").GetData());
            renderer->AddActor2D(textActor);

            renderWindow->SetWindowName(DataName.c_str());
            DoOnce = false;
        }
        std::string tString; tString="t="+std::to_string(physics->time+physics->timeScheme->dt);
        textActor->SetInput(tString.c_str());
        renderer->ResetCamera();
        renderWindow->Render();
    }
    MPI_Barrier(PETSC_COMM_WORLD);
}
#endif