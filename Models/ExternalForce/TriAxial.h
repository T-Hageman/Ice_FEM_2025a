#ifndef TRIAXIALMODEL_H
#define TRIAXIALMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Fracture/PhaseField/PhaseFieldUtil.h"

class TriAxialModel: public BaseModel {
    public:
        std::string Name;
        std::vector<std::string> DofNames_u = {"ux","uy","uz"};      //names for displacement degrees of freedom
        std::string Dofname_PF = "phase";

        PhaseFieldUtility* PF_Util;

        TriAxialModel(Physics& My_Physics, std::string MyName);
        ~TriAxialModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

        size_t hasTimeData(std::vector<std::string>& DataNames);
        void GetTimeData(std::vector<double>& DataNames);
        void Commit(int CommitType);

    protected:

    private:
        size_t ElementGroupIndex_Side, ElementGroupIndex_Top;
        size_t dofStep, dofStep_Phase, dofType_phase;
        std::vector<size_t> dofTypes;

        bool DamageForces;

        double p0;
        double uRate;

        double F_Ext, F_ExtOld;
        double uz_Sum, uz_SumOld;
        double Area;
        double u0;
        double tbreak, maxLoad; 

        double kdummy;
        double t_postFail;
        bool MinLoadPassed;
};

void Register_TriAxialModel();
BaseModel* New_TriAxialModel(Physics& My_Physics, std::string MyNameIn);

#endif

