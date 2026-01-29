#include "Plot_IPSurf_3D.h"
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
#include <vtkContourFilter.h>
#include <vtkPerlinNoise.h>
#include <vtkSampleFunction.h>
#include <vtkPointInterpolator.h>
#include <vtkPointGaussianMapper.h>
#include <vtkPointInterpolator.h>
#include <vtkGaussianKernel.h>
#include <vtkImageData.h>

Plot_IPSurf_3D::Plot_IPSurf_3D(std::string name, inputData &inputs, Physics *phys): Vis_Plot(name, inputs, phys) {
    std::string ElemGroupIndex_Name; inputs.GetRequired(ElemGroupIndex_Name, {"Visualisation",Name,"ElementGroup"});
    ElemGroupIndex_Data = physics->mesh->GetElementGroupIdx(ElemGroupIndex_Name);

    inputs.GetRequired(DataName, {"Visualisation",Name,"IP"});

    LogScale = false;
    inputs.GetOptional(LogScale, {"Visualisation",Name,"LogScale"});

    sliceDir = "x";
    inputs.GetOptional(sliceDir, {"Visualisation",Name,"Slice"});

    if (rank == 0){
        // Set the background color.
        colors->SetColor("BkgColor", colors->GetColor4d("Snow").GetData());
        renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
        renderWindow->AddRenderer(renderer);
        renderer->ResetCamera();
    }
}

Plot_IPSurf_3D::~Plot_IPSurf_3D(){

}

void Plot_IPSurf_3D::Plot(){
    std::vector<std::vector<double>> XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords;

    size_t nNodes_per_Elem = physics->mesh->GetExportMeshCoords(ElemGroupIndex_Data, XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords);
    size_t nElems = XCoords.size();
    size_t nIP = physics->mesh->ElementGroups[ElemGroupIndex_Data].BaseElem->ipcount;

    std::vector<std::vector<double>> IPData;
    IPData.resize(nElems); 
    for (size_t i = 0; i < nElems; i++) IPData[i].resize(nIP);

    physics->GetIPDataToSave(DataName, ElemGroupIndex_Data, IPData);
    if (LogScale){
        GetLog10(IPData);
    }

    std::vector<double> DRange(2); DRange[0] = min(IPData); DRange[1] = max(IPData);
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

    CreatePlot(XIPCoords, YIPCoords, ZIPCoords, IPData, colorLookupTable);

    MPI_Barrier(PETSC_COMM_WORLD);
    
    if (rank==0){
        ColourBar->SetLookupTable(colorLookupTable);
        if (DoOnce){
            if (true){
                Mapper2->SetInputData(ugrid);
                Mapper2->ScalarVisibilityOn();
                contourActor->SetMapper(Mapper2);
            } else {
                contourMapper->SetInputConnection(SliceFilter->GetOutputPort());
                contourMapper->ScalarVisibilityOn();
                contourActor->SetMapper(contourMapper);
            }

            contourActor->GetProperty()->EdgeVisibilityOn();

            renderer->AddActor(contourActor);

            std::vector<double> XRange(2);
            std::vector<double> YRange(2);
            std::vector<double> ZRange(2);
            XRange[0] = polyData->GetBounds()[0];
            XRange[1] = polyData->GetBounds()[1];
            YRange[0] = polyData->GetBounds()[2];
            YRange[1] = polyData->GetBounds()[3];
            ZRange[0] = polyData->GetBounds()[4];
            ZRange[1] = polyData->GetBounds()[5];

            renderer->GetActiveCamera()->ParallelProjectionOn();
            if (sliceDir=="x"){
                renderer->GetActiveCamera()->SetPosition(-1, 0, 0);
                renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
                renderer->GetActiveCamera()->SetViewUp(0,0,1);
            } else if (sliceDir=="y"){
                renderer->GetActiveCamera()->SetPosition(-1, 1, 0);
                renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
                renderer->GetActiveCamera()->SetViewUp(0,0,1);
            } else {
                renderer->GetActiveCamera()->SetPosition(0, 0, -1);
                renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
                renderer->GetActiveCamera()->SetViewUp(0,1,0);
            }

            AddGrid(polyData->GetBounds());
            AddColourBar();

            renderWindow->SetWindowName(DataName.c_str());
            // Set the icon for the window
            auto icon = LoadIconFromEmbeddedPNG();
            renderWindow->SetIcon(icon);

            DoOnce = false;
        }
        if (false){
            contourMapper->Update();
        } else {
            Mapper2->Update();
        }
        renderer->ResetCamera();
        renderWindow->Render();
    }

    MPI_Barrier(PETSC_COMM_WORLD);
}

void Plot_IPSurf_3D::CreatePlot(std::vector<std::vector<double>> &XCoords, std::vector<std::vector<double>> &YCoords, std::vector<std::vector<double>> &ZCoords, 
                          std::vector<std::vector<double>> &Data,
                          vtkNew<vtkLookupTable>& colorLookupTable){
    size_t nElems, nNodes_per_Elem;
    nElems = XCoords.size();
    if (nElems>0){
        nNodes_per_Elem = XCoords[0].size();
    } else {
        nNodes_per_Elem = 0;
    }
    size_t nPoints = nElems*nNodes_per_Elem;

    std::vector<double> x(nPoints), y(nPoints), z(nPoints), d(nPoints);
    size_t p=0;
    for (size_t i = 0; i < nElems; i++){
        for (size_t j = 0; j < nNodes_per_Elem; j++){
            x[p] = XCoords[i][j];
            y[p] = YCoords[i][j];
            z[p] = ZCoords[i][j];
            d[p] = Data[i][j];
            p += 1;
        }
    }

    std::vector<MPI_Request> MPI_Reqs(4);
    std::vector<MPI_Status> MPI_Stats(4);
    if (rank==0){
        if (size==1){
            scalars->SetNumberOfComponents(3);
            scalars->Allocate(nPoints);
            points->Allocate(nPoints);
            //scalars->SetName("phi");
            
            for (size_t i = 0; i < nPoints; i++){
                points->InsertPoint(i, x[i], y[i], z[i]);

                double dcolor[3];
                unsigned char ccolor[3];
                colorLookupTable->GetColor(d[i], dcolor);
                for (size_t k = 0; k < 3; k++){
                    ccolor[k] = static_cast<unsigned char>(255.0 * dcolor[k]);
                }
                scalars->InsertTuple3(i, ccolor[0], ccolor[1], ccolor[2]);
            }
        } else {
            std::vector<size_t> Points_per_core(size);
            MPI_Gather(&nPoints, 1, my_MPI_SIZE_T, Points_per_core.data(), 1, my_MPI_SIZE_T, 0, PETSC_COMM_WORLD);

            size_t TotalPointCount = sum(Points_per_core);
            scalars->SetNumberOfComponents(3);
            scalars->Allocate(TotalPointCount);
            points->Allocate(TotalPointCount);
            scalars->SetName("Phi"); 

            p = 0;
            for (size_t c = 0; c < size; c++){
                if (c!=0){
                    nPoints = Points_per_core[c];
                    x.resize(nPoints);
                    y.resize(nPoints);
                    z.resize(nPoints);
                    d.resize(nPoints);

                    MPI_Irecv(x.data(), nPoints, MPI_DOUBLE, c, c*10+1, PETSC_COMM_WORLD, &MPI_Reqs[0]);
                    MPI_Irecv(y.data(), nPoints, MPI_DOUBLE, c, c*10+2, PETSC_COMM_WORLD, &MPI_Reqs[1]);
                    MPI_Irecv(z.data(), nPoints, MPI_DOUBLE, c, c*10+3, PETSC_COMM_WORLD, &MPI_Reqs[2]);
                    MPI_Irecv(d.data(), nPoints, MPI_DOUBLE, c, c*10+4, PETSC_COMM_WORLD, &MPI_Reqs[3]);
                    MPI_Waitall(4, MPI_Reqs.data(), MPI_Stats.data());
                }
                for (size_t i = 0; i < nPoints; i++){
                    if (abs(x[i])<0.05e-2){
                        if (DoOnce){
                            points->InsertPoint(p, x[i], y[i], z[i]);
                        }

                        double dcolor[3];
                        unsigned char ccolor[3];
                        colorLookupTable->GetColor(d[i], dcolor);
                        for (size_t k = 0; k < 3; k++){
                            ccolor[k] = static_cast<unsigned char>(255.0 * dcolor[k]);
                        }
                        scalars->InsertTuple3(p, ccolor[0], ccolor[1], ccolor[2]);

                        p += 1;
                    }
                }
            }
        }

        // Step 3: Create a VTK PolyData to store points and scalars
        if (DoOnce){
            std::cout << "start DOONCE\n";
            polyData->SetPoints(points);
            polyData->GetPointData()->SetScalars(scalars);
            polyData->GetPointData()->Modified();

            delaunay->SetInputData(polyData);
            delaunay->Modified();
            delaunay->Update();

            ugrid->DeepCopy(delaunay->GetOutput());

            CutPlane->SetOrigin(0,0,0);
            if (sliceDir=="x"){
                CutPlane->SetNormal(0,1,0);
            }
            if (sliceDir == "y"){
                CutPlane->SetNormal(1,0,0);
            } else {
                CutPlane->SetNormal(0,0,1);
            }

            SliceFilter->SetCutFunction(CutPlane);
            SliceFilter->SetInputData(ugrid);
            SliceFilter->Update();
        } else {
            polyData->GetPointData()->SetScalars(scalars);
            ugrid->GetPointData()->SetScalars(cleanPolyData->GetOutput()->GetPointData()->GetScalars());
            SliceFilter->SetInputData(ugrid);
        }
        SliceFilter->Modified();
        SliceFilter->Update();

    } else {
        MPI_Gather(&nPoints, 1, my_MPI_SIZE_T, NULL, 1, my_MPI_SIZE_T, 0, PETSC_COMM_WORLD);

        MPI_Isend(x.data(), nPoints, MPI_DOUBLE, 0, rank*10+1, PETSC_COMM_WORLD, &MPI_Reqs[0]);
        MPI_Isend(y.data(), nPoints, MPI_DOUBLE, 0, rank*10+2, PETSC_COMM_WORLD, &MPI_Reqs[1]);
        MPI_Isend(z.data(), nPoints, MPI_DOUBLE, 0, rank*10+3, PETSC_COMM_WORLD, &MPI_Reqs[2]);
        MPI_Isend(d.data(), nPoints, MPI_DOUBLE, 0, rank*10+4, PETSC_COMM_WORLD, &MPI_Reqs[3]);
    }

    MPI_Barrier(PETSC_COMM_WORLD);
}
#endif