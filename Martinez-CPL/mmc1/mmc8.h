#ifndef SOMEDEFS
#define SOMEDEFS

#include <math.h>
#include <stdio.h>

#include <vector>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Core>




typedef Eigen::SparseMatrix<double> MSparseMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> MTriplet;  //triplets to fill sparse matrix
typedef std::vector<MTriplet> Container;   //container for triplets. Sparse matrix is filled from Container
typedef Eigen::SparseLU<Eigen::SparseMatrix<double> > SolverLU;  //sparse solver



#define MAX_SPECIES 10
#define MAX_STEPS_C 10
#define MAX_POINTS_DERIV 7

#define EMECHANISM 0
#define ECMECHANISM 1

#define PLANE 0
#define CYLINDRICAL 1
#define SPHERICAL 2


//for first order kinetics

struct TChemKineticOrder1
{
    int TotalTerms;
    int Index[MAX_STEPS_C];
    double Terms[MAX_STEPS_C];

};

//Data to manage simulations
struct GeneralData
{
    double Gamma;
    double hh;
    double Lambda;
    double Cinfi[MAX_SPECIES];
    double XMax;
    int NPointsGridSpecies;
    int PointsDerivs;
    int PointsCurrent;
    int OneRow;
    int Mechanism;
    int Electrode;
    int NSpecies;
    int TimeSteps;
    double FactorG;
    double DeltaT;
    double R0;
    double KKinec[MAX_STEPS_C];


};






#endif // SOMEDEFS

