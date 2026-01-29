#ifndef PLOT_BoundaryFill_H
#define PLOT_BoundaryFill_H

#include "Visualisation.h"
#ifdef LOCALENVIRONMENT
class Plot_BoundaryFill: public Vis_Plot{

    public:
        Plot_BoundaryFill(std::string name, inputData& inputs, Physics* phys);
        ~Plot_BoundaryFill();
        void Plot();



    private:
        std::vector<std::string> DofNames_u = {"ux","uy","uz"};
        size_t dofStep_u;
        std::vector<size_t> dofTypes_u;

        std::string DataName;

        size_t ElemGroupIndex_u, ElemGroupIndex_Data;
        bool Deformations;
        double DeformationScale;

        std::string GrowDir;
        double DataScale;
        bool LogScale;

        //display mesh data
        vtkNew<vtkPoints> points;
        vtkNew<vtkUnstructuredGrid> ugrid;
        vtkNew<vtkDataSetMapper> ugridMapper;
        vtkNew<vtkActor> ugridActor;
        vtkNew<vtkScalarBarActor> ColourBar;
};

#endif
#endif