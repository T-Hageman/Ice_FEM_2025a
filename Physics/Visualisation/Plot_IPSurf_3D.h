#ifndef PLOT_IPSURF3D_H
#define PLOT_IPSURF3D_H

#ifdef LOCALENVIRONMENT
#include "Visualisation.h"

#include <vtkContourFilter.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkCleanPolyData.h>

class Plot_IPSurf_3D: public Vis_Plot{

    public:
        Plot_IPSurf_3D(std::string name, inputData& inputs, Physics* phys);
        ~Plot_IPSurf_3D();
        void Plot();

    private:
        std::string sliceDir;

        vtkNew<vtkFloatArray> scalars;
        void CreatePlot(std::vector<std::vector<double>> &XCoords, std::vector<std::vector<double>> &YCoords, std::vector<std::vector<double>> &ZCoords, 
                          std::vector<std::vector<double>> &Data,
                          vtkNew<vtkLookupTable>& colorLookupTable);

        vtkNew<vtkPolyData> polyData;
        vtkNew<vtkCleanPolyData> cleanPolyData;

        vtkNew<vtkCutter> SliceFilter;
        vtkNew<vtkPlane> CutPlane;

        vtkNew<vtkPolyDataMapper> contourMapper;
        vtkNew<vtkDataSetMapper> Mapper2;
        vtkNew<vtkActor> contourActor;

        vtkNew<vtkDelaunay3D> delaunay;
};

#endif
#endif