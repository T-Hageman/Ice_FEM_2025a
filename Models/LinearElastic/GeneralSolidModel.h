#ifndef GENERALSOLIDMODEL_H
#define GENERALSOLIDMODEL_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Materials/SolidMaterial.h"

class GeneralSolidModel: public BaseModel{
    public:
        std::vector<std::string> DofNames_u = {"ux","uy"};      //names for displacement degrees of freedom

        GeneralSolidModel(Physics& My_Physics, std::string MyName);
        ~GeneralSolidModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);
        void ResetStep();

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat &K, Vec &f, Constrainer *cons, size_t step);
        void Commit(int CommitType);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>> &Data);
        void SaveStressComponent(std::string DataName, std::vector<std::vector<double>> &Data, bool Damaged);
        void SaveStrainComponent(std::string DataName, std::vector<std::vector<double>> &Data);
    protected:

    private:
        void Assemble_U(Vec &f, Mat &K);

        std::string SolidMatName;   //name of solid material
        SolidMaterial* M;

        size_t ElemGroupIndex_u;  //index of the element groups
        std::vector<size_t> dofTypes_u;                 // dof indices for displacements
        size_t Step_u;                      // steps in which displacements and phasefield are resolved (NOT THE SAME STEP)

        Eigen::Vector3d g={0, -9.81, 0};   //gravity vector

        bool IncludeInertia;

        std::vector<std::vector<Vector6d>> Plastic_Strain, Plastic_StrainOld;    //plastic strain history
        std::vector<std::vector<std::vector<double>>> MatHist, MatHistOld;
};

void Register_GeneralSolidModel();
BaseModel* New_GeneralSolidModel(Physics& My_Physics, std::string MyNameIn);

#endif

