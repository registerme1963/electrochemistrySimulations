#ifndef SIMUL_H
#define SIMUL_H
#include "somedefs.h"
//objet to manage simulations
class Simul
{
private:
    struct GeneralData MyData;
    Container MContainer;
    MSparseMat MatrixA;
    SolverLU MSolver;
    Eigen::VectorXd Vectorb, Vectorx;
    double *XX;
    double *Gamma2i;
    double *Gammai;
    double *TermIndp;
    double CoefAlpha[MAX_POINTS_DERIV];
    double CoefBeta[MAX_POINTS_DERIV];
    double CoefBeta0[MAX_POINTS_DERIV];
    double *Concentration;
    TChemKineticOrder1 TKinO1[MAX_SPECIES];

public:
    Simul(struct GeneralData Data);
    ~Simul();
    void FillMatrixA(double Potential);
    void BulkConCen();
    void ObtainAlphaAndBeta();
    double CoeffMatrixN2(int ii, int PosRel);
    void ObtainKineticTerms();
    void AddKineticToMatrix();
    void FillIndTerm();
    double ObtainBN2(int ii,int IndexSpecies);
    void SolveLinearSistem(void);
    double ObtainFlux(int Species);
    void ChangeNernstPotential(int row, int column, double Potential);
    void AnalyzeMatrixPattern(void);
    void FactorizeMatrix(void);
    void PrintConcen(char *fichero);
    //void PrintB(char *fichero);
    void SurfaceConditions(double Potential);
    void Nernst(int Row,int Ox,int Red,double Potential);
    void Fluxes(int Row,int NEspecies,int *Index);
};

#endif // SIMUL_H
