#ifndef _Sound_and_LPC2_h_
#define _Sound_and_LPC2_h_
/* Sound_and_LPC2.h
 *
 * Copyright (C) 1994-2025 David Weenink
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 djmw 19971103
 djmw 20020812 GPL header
*/

#include "SoundFrameIntoSampledFrame2.h"
#include "Sound.h"
#include "LPC.h"
#include "SVD.h"

/*
	20240603:
	The output of the Sound_to_LPC_<x> might be a little bit different from previous outputs because:
	1. The sound frame now always has an odd number of samples, irrespective of the sampling frequency or the window shape.
		This means that also the window function will have an odd number of samples. Therefore the sample at the centre of the sound frame always has weight 1.0.
		Previously the number of samples could be even or odd, depending on how rounding turned out.
	2. The Gaussian window function was slightly improved.
	3. The precision of the autocorrelation and covariance method have been improved a little by using some `longdouble` accumulators.
*/

Thing_define (SoundFrameIntoLPCFrame2, SoundFrameIntoSampledFrame2) {
	
	mutableLPC outputLPC;
	integer order;
	integer currentOrder;  // TODO djmw 20250825 is this one necessary 
	integer orderp1; 	// convenience order+1
	autoVEC a;			// common work vector of dimension orderp1

	virtual void initBasicSoundFrameAndLPC (constSound input, mutableLPC outLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape);
	
	void initHeap ()
		override;
};

/*********************** Autocorrelation method *************************************************************/

Thing_define (SoundFrameIntoLPCFrame2Auto, SoundFrameIntoLPCFrame2) {
	
	autoVEC r;		// orderp1
	autoVEC rc;		// orderp1
	
	void initHeap ()
		override;
		
	bool inputFrameIntoOutputFrame (integer currentFrame)
		override;
};

autoSoundFrameIntoLPCFrame2Auto SoundFrameIntoLPCFrame2Auto_create (constSound input, mutableLPC outLPC,
	double effectiveAnalysisWidth, kSound_windowShape windowShape);


void Sound_into_LPC_auto2 (constSound me, mutableLPC thee, double analysisWidth);

autoLPC Sound_to_LPC_auto2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency);

/*********************** Covariance method *************************************************************/

Thing_define (SoundFrameIntoLPCFrame2Covar, SoundFrameIntoLPCFrame2) {
	integer order2;	// size: order * (order + 1) / 2
	autoVEC b;		// size: order2
	autoVEC grc;	// size: order
	autoVEC beta;	// size: order
	autoVEC cc;		// size: orderp1
	
	void initBasicSoundFrameAndLPC (constSound input, mutableLPC outLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape) override;

	void initHeap ()
		override;
	
	bool inputFrameIntoOutputFrame (integer currentFrame)
		override;
};

autoSoundFrameIntoLPCFrame2Covar SoundFrameIntoLPCFrame2Covar_create (constSound input, mutableLPC outLPC,
	double effectiveAnalysisWidth, kSound_windowShape windowShape);

void Sound_into_LPC2_covar (constSound me, mutableLPC thee, double analysisWidth);

autoLPC Sound_to_LPC2_covar (constSound me, mutableLPC thee, double analysisWidth);

/*********************** Burg method *************************************************************/

Thing_define (SoundFrameIntoLPCFrame2Burg, SoundFrameIntoLPCFrame2) {
	autoVEC b1;		// size: soundFrameSize
	autoVEC b2;		// size: soundFrameSize
	autoVEC aa;		// size: order

	void initHeap ()
		override;

	double burg (VEC const& a, constVEC const& x, integer& frameAnalysisInfo);

	bool inputFrameIntoOutputFrame (integer currentFrame)
		override;

};

autoSoundFrameIntoLPCFrame2Burg SoundFrameIntoLPCFrame2Burg_create (constSound input, mutableLPC outLPC,
	double effectiveAnalysisWidth, kSound_windowShape windowShape);

void Sound_into_LPC_burg (constSound me, mutableLPC thee, double analysisWidth);

autoLPC Sound_to_LPC_burg (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency);

/*********************** Marple method *************************************************************/

Thing_define (SoundFrameIntoLPCFrame2Marple, SoundFrameIntoLPCFrame2) {
	double tol1, tol2;
	autoVEC c;	// orderp1)
	autoVEC d;	// orderp1)
	autoVEC r;	// orderp1)
	
	void initBasicSoundFrameAndLPCMarple (constSound inputSound, mutableLPC outputLPC, 
		double effectiveAnalysisWidth, kSound_windowShape windowShape, double tol1, double tol2);
	
	void initHeap ()
		override;
	
	bool inputFrameIntoOutputFrame (integer iframe)
		override;
		
};

autoSoundFrameIntoLPCFrame2Marple SoundFrameIntoLPCFrame2Marple_create (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, double preEmphasisFrequency, double tol1, double tol2);

void Sound_into_LPC_marple2 (constSound me, mutableLPC thee, double analysisWidth, double tol1, double tol2);

autoLPC Sound_to_LPC_marple2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt, 
	double preEmphasisFrequency, double tol1, double tol2);

/*********************** Robust method (LPC & Sound) *************************************************************/

Thing_define (LPCAndSoundFramesIntoLPCFrameRobust2, SoundFrameIntoLPCFrame2) {
	constLPC inputLPC;
	mutableLPC outputLPC;
	constSound inputSound;
	
	//oo_STRUCT (LPC_Frame, otherInputLPCFrame)
	
	integer currentPredictionOrder;
	double k_stdev;
	integer iter;
	integer itermax;
	integer huber_iterations; // = 5;
	bool wantlocation;	//
	bool wantscale;	//
	double location = 0.0;
	double scale;
	double tol1;
	double tolSVD = 1e-10;
	autoVEC error; //  soundFrameSize)
	autoVEC sampleWeights;	// soundFrameSize
	autoVEC coefficients;	// inputLPC -> maxnCoefficients
	autoVEC covariancesw;	// inputLPC -> maxnCoefficients
	autoMAT covarmatrixw;	// inputLPC -> maxnCoefficients, inputLPC -> maxnCoefficients
	integer computedSVDworksize;
	autoSVD svd;
	autoVEC svdwork1;	// computedSVDworksize)
	autoVEC svdwork2;	// order)
	autoVEC filterMemory;	// order)
	autoVEC huberwork;	// soundFrameSize)

	void initHeap ()
		override;

	bool inputFrameIntoOutputFrame (integer iframe)
		override;

private:

	void resize ();
	void setSampleWeights ();
	void setCovariances ();
	void solvelpc ();
};

autoLPCAndSoundFramesIntoLPCFrameRobust2 LPCAndSoundFramesIntoLPCFrameRobust2_create (constLPC inputLPC,
	constSound inputSound, mutableLPC outputLPC, double effectiveAnalysisWidth, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation);

void LPC_and_Sound_into_LPC_robust2 (constLPC inputLPC, constSound inputSound, mutableLPC outpuLPC, double effectiveAnalysisWidth, 
	double preEmphasisFrequency, double k_stdev, integer itermax, double tol, bool wantlocation);

autoLPC LPC_and_Sound_to_LPC_robust2 (LPC inputLPC, constSound inputSound, double effectiveAnalysisWidth,
	double preEmphasisFrequency, double k_stdev, integer itermax, double tol, bool wantlocation);

/*********************** Robust method (Sound) *************************************************************/

Thing_define (SoundFrameIntoLPCFrameRobust2, SoundFrameIntoLPCFrame2Auto) {
	autoLPCAndSoundFramesIntoLPCFrameRobust2 lpcAndSoundFramesIntoLPCframe;
	
	bool inputFrameIntoOutputFrame (integer iframe)
		override;
	
	void initHeap ()
		override;
};

autoSoundFrameIntoLPCFrameRobust2 SoundFrameIntoLPCFrameRobust2_create (constSound inputSound, mutableLPC outputLPC, mutableLPC inputLPC,
	double effectiveAnalysisWidth, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation);

void Sound_into_LPC_robust2 (constSound me, mutableLPC thee, double analysisWidth, kSound_windowShape windowShape,
	double k_stdev,	integer itermax, double tol, bool wantlocation);

autoLPC Sound_to_LPC_robust2 (constSound me, int predictionOrder, double effectiveAnalysisWidth, double dt,
	double preEmphasisFrequency, kSound_windowShape windowShape, double k_stdev, integer itermax, double tol, bool wantlocation);


#if 0

/*
	Precondition:
		Sound has been resampled and pre-emphasized
*/

void Sound_into_LPC_robust (constSound me, mutableLPC thee, double analysisWidth,
	double k_stdev,	integer itermax, double tol, bool wantlocation);

/*
 * Function:
 *	Calculate linear prediction coefficients according to following model:
 *  Minimize E(m) = Sum(n=n0;n=n1; (x [n] + Sum(k=1;k=m; a [k]*x [n-k])))
 * Method:
 *  The minimization is carried out by solving the equations:
 *  Sum(i=1;i=m; a [i]*c [i] [k]) = -c [0] [k] for k=1,2,...,m
 *  where c [i] [k] = Sum(n=n0;n=n1;x [n-i]*x [n-k])
 *  1. Covariance:
 *		n0=m; n1 = N-1;
 *      c [i] [k] is symmetric, positive semi-definite matrix
 *  	Markel&Gray, LP of Speech, page 221;
 *  2. Autocorrelation
 *		signal is zero outside the interval;
 *      n0=-infinity; n1=infinity
 *      c [i] [k] symmetric, positive definite Toeplitz matrix
 *  	Markel&Gray, LP of Speech, page 219;
 * Preconditions:
 *	predictionOrder > 0;
 *  preEmphasisFrequency >= 0;
 *
 * Burg method: see Numerical recipes Chapter 13.
 *
 * Marple method: see Marple, L. (1980), A new autoregressive spectrum analysis
 *		algorithm, IEEE Trans. on ASSP 28, 441-453.
 *	tol1 : stop iteration when E(m) / E(0) < tol1
 *	tol2 : stop iteration when (E(m)-E(m-1)) / E(m-1) < tol2,
 */

autoSound LPC_Sound_filter (constLPC me, constSound thee, bool useGain);
/*
	E(z) = X(z)A(z),
	A(z) = 1 + Sum (k=1, k=m, a(k)z^-k);

	filter:
		given e & a, determine x;
		x(n) = e(n) - Sum (k=1, m, a(k)x(n-k))
	useGain determines whether the LPC-gain is used in the synthesis.
*/

void LPC_Sound_filterWithFilterAtTime_inplace (constLPC me, mutableSound thee, integer channel, double time);

autoSound LPC_Sound_filterWithFilterAtTime (constLPC me, constSound thee, integer channel, double time);

autoSound LPC_Sound_filterInverse (constLPC me, constSound thee);
/*
	E(z) = X(z)A(z),
	A(z) = 1 + Sum (k=1, k=m, a(k)z^-k);

	filter inverse:
		given x & a, determine e;
		e(n) = x(n) + Sum (k=1, m, a(k)x(n-k))
*/

autoSound LPC_Sound_filterInverseWithFilterAtTime (constLPC me, constSound thee, integer channel, double time);

void LPC_Sound_filterInverseWithFilterAtTime_inplace (constLPC me, mutableSound thee, integer channel, double time);

/*
	For all LPC analysis
*/
void checkLPCAnalysisParameters_e (double sound_dx, integer sound_nx, double physicalAnalysisWidth, integer predictionOrder);
#endif

#endif /* _Sound_and_LPC2_h_ */
