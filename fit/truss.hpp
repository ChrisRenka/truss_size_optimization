#ifndef _TRUSS_HPP
#define _TRUSS_HPP

#include<stdio.h>
#include<cmath>

#define PI 3.141592653598

//Finite element analysis and structure modelling were developed based on:
//https://github.com/labicom-udesc/ECO-framework/blob/master/ProblemFunc.c

void setupTruss(int nTruss, int nThreads, int *avalPerInd);

void setupTruss10(int nThreads, int *avalPerInd);
void setupTruss17(int nThreads, int *avalPerInd);
void setupTruss200(int nThreads, int *avalPerInd);

void setupTruss25(int nThreads, int *avalPerInd);
void setupTruss72(int nThreads, int *avalPerInd);
void setupTruss120(int nThreads, int *avalPerInd);

void setupTruss72case1(int nThreads, int *avalPerInd);
void setupTruss72case2(int nThreads, int *avalPerInd);

void distributeAreasEqual(double *sol, double *area, int id);
void distributeAreas200(double *sol, double *area, int id);
void distributeAreas25(double *sol, double *area, int id);
void distributeAreas72(double *sol, double *area, int id);
void distributeAreas120(double *sol, double *area, int id);

void destroyTruss();
double trussFit(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);
double trussFit2D(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);
double trussFit3D(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);

double fit25(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);
double fit72(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);
double fit120(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);
double fit200(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id);

double getMinArea();
double getMaxArea();
int getDim();

#endif