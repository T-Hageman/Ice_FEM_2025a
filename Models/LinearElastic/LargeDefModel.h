#ifndef LARGEDEFMODEL_H
#define LARGEDEFMODEL_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "../Materials/SolidMaterial.h"

class LargeDeformationModel: public BaseModel{
    public:
        std::vector<std::string> DofNames_u = {"ux","uy"};      //names for displacement degrees of freedom

        LargeDeformationModel(Physics& My_Physics, std::string MyName);
        ~LargeDeformationModel();

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
};

void Register_LargeDeformationModel();
BaseModel* New_LargeDeformationModel(Physics& My_Physics, std::string MyNameIn);

#endif

