#ifndef _LSHADE_HPP
#define _LSHADE_HPP

#include<stdio.h>
#include<omp.h>
#include<sys/time.h>
#include<random>
#include<vector>
#include<algorithm>
#include<thread>
#include"../fit/truss.hpp"

using namespace std;

typedef struct{
    vector<double> data;
    double fit;
    double weight;
    double pen;
    double cr;
    double f;
    int id;
    int mut;
} Individual;

typedef struct{
    vector<Individual*> ind;
    vector<Individual*> trial;
    vector<Individual*> aux;
    vector<double> sf;
    vector<double> scr;
    vector<double> fitDiff;
    vector<double> mf;
    vector<double> mcr;
    vector<double> mfcp;
    int np;
    int npi;
    int npf;
    int d;
    int aval;
    int k;
    int h;
    int iBest;
    double li;
    double ui;
    double p;
    double rArc;
    double c;
    double sumDiffM1;
    double sumDiffM2;
} Population;

//vec.reserve(n);

//Main functions
Population *initializePop(FILE *input);
void mainPopulationFitness(Population *pop, double (*p) (Individual*, int));
void trialPopulationFitness(Population *pop, double (*p) (Individual*, int));
void generateTrialPop(Population *pop);
void mutation(Population *pop, int ind_parent, int *ind_mut, void (*p) (Population*, int, int*));
void selectionMin(Population *pop);
void popAval(Population *pop, bool mainPop, int nThreads);
//Main functions end

//JADE/SHADE functions
vector<int> findPBest(Population *pop);
bool compareByFitness(Individual *a, Individual *b);
double wArithmeticMean(vector<double> &v, vector<double> &d, double sumDiff);
double wLehmerMean(vector<double> &v, vector<double> &d, double sumDiff);
void updateMeans(Population *pop);
void mutationCurrentToPBest(Population *pop, int ind_parent, int *ind_mut);
void mutationCurrentToOrdPBest(Population *pop, int ind_parent, int *ind_mut);
void copyIndividualToAux(Population *pop, int index);
void adjustSettings(Population *pop);
void linearPopulationReduction(Population *pop, int currentAval);
void alternativePopulationReduction(Population *pop, int currentAval);
//JADE functions end

//Crowding functions
double euclideanDistance(Individual *ind1, Individual *ind2);
double* nearestNeighbor(Population *pop, int trialIndex);
void selectionCrowdingMin(Population *pop);
void selectionNGMin(Population *pop, int g, double *avg);
//Crowding functions end

//Clean-up fuctions
void destroyIndividual(Individual *ind);
void destroyPopulation(Population *pop);
//Clean-up functions end

//Input function
void parametersInput(Population *pop, FILE *arq);
//Input function end

//Auxiliary Functions
double meanFitness(Population *pop);
double standardDeviationFitness(Population *pop, double mean);
double standardDeviationRuns(vector<vector<double>> &meanRuns, double mean, int runs, int gen);
double meanArray(double *arr, int size);
double standardDeviationArray(double *arr, int size);
double wtime();
void copyIndividual(Individual *original, Individual *copy);
//Auxiliary Functions End

//Output functions
double populationDiversity(Population *pop);
int bestFitnessIndex(Population *pop);
//Output functions end

//Visualization Functions
void printDataIndividual(Individual *ind);
void printDataPop(Population *pop);
void printDataTrial(Population *pop);
void printFitIndividual(Individual *ind);
void printFitPop(Population *pop);
//Visualization Functions end

#endif
