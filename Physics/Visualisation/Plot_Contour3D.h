#ifndef PLOT_CONTOUR3D_H
#define PLOT_CONTOUR3D_H

#ifdef LOCALENVIRONMENT
#include "Visualisation.h"

#include <vtkContourFilter.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>

class Plot_Contour3D: public Vis_Plot{

    public:
        Plot_Contour3D(std::string name, inputData& inputs, Physics* phys);
        ~Plot_Contour3D();
        void Plot();

    private:
        vtkNew<vtkFloatArray> scalars;
        void CreateContour(std::vector<std::vector<double>> &XCoords, std::vector<std::vector<double>> &YCoords, std::vector<std::vector<double>> &ZCoords, 
                          std::vector<std::vector<double>> &Data);

        vtkNew<vtkPolyData> polyData;
        vtkNew<vtkContourFilter> contourFilter;

        vtkNew<vtkPolyDataMapper> contourMapper;
        vtkNew<vtkActor> contourActor;

        vtkNew<vtkDelaunay3D> delaunay;
};

#endif
#endif