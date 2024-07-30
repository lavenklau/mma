#pragma once

#ifndef __MMAOPTER_HPP
#define __MMAOPTER_HPP


#ifdef _EXPORT_MMAOPT
#ifdef _WIN32
#define API_MMAOPT extern "C" __declspec(dllexport)
#else
#define API_MMAOPT __attribute__((visibility("default")))
#endif
#else
#ifdef _WIN32
#define API_MMAOPT extern "C" __declspec(dllimport)
#else
#define API_MMAOPT
#endif
#endif


API_MMAOPT void mmasub(int ncontrain, int nvar, int itn, double* xvar, double* xmin, double* xmax, double* xold1, double* xold2,
	double f0val, double* df0dx, double* gval, double* dgdx, double* low, double* upp,
	double a0, double* a, double* c, double* d, double move);


#endif
