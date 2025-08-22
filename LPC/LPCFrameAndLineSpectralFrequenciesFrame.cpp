/* LPCFrameIntoLineSpectralFrequenciesFrame.cpp
 *
 * Copyright (C) 2024,2025 David Weenink
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#include "melder.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame.h"
#include "Roots_and_Formant.h"

#include "oo_DESTROY.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_COPY.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_EQUAL.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_CAN_WRITE_AS_ENCODING.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_WRITE_TEXT.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_WRITE_BINARY.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_READ_TEXT.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_READ_BINARY.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"
#include "oo_DESCRIPTION.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"

Thing_implement (LPCFrameIntoLineSpectralFrequenciesFrame, LPCFrameIntoSampledFrame, 0);

void structLPCFrameIntoLineSpectralFrequenciesFrame :: allocateOutputFrames () {
	for (integer iframe = 1; iframe <= output -> nx; iframe ++) {
		LineSpectralFrequencies_Frame lsff = & outputLSF -> d_frames [iframe];
		LineSpectralFrequencies_Frame_init (lsff, outputLSF -> maximumNumberOfFrequencies);
	}
}

/*
	Conversion from Y(w) to a polynomial in x (= 2 cos (w))
	From: Joseph Rothweiler (1999), "On Polynomial Reduction in the Computation of LSP Frequencies."
	IEEE Trans. on ASSP 7, 592--594.
*/
static void cos2x (VECVU const& g) {
	for (integer i = 3; i <= g.size; i ++) {
		for (integer j = g.size; j > i; j --)
			g [j - 2] -= g [j];
		g [i - 2] -= 2.0 * g [i];
	}
}

static void LPC_Frame_into_Polynomial_sum (LPC_Frame me, Polynomial psum) {
	/*
		Fs (z) = A(z) + z^-(p+1) A(1/z)
	*/
	const integer order = my nCoefficients, g_order = (order + 1) / 2; // orders
	psum -> coefficients [order + 2] = 1.0;
	for (integer i = 1; i <= order; i ++)
		psum -> coefficients [order + 2 - i] = my a [i] + my a [order + 1 - i];

	psum -> coefficients [1] = 1.0;
	psum -> numberOfCoefficients = order + 2;

	if (order % 2 == 0) // order even
		Polynomial_divide_firstOrderFactor (psum, -1.0, nullptr);
	/*
		Transform to cos(w) terms
	*/
	for (integer i = 1; i <= g_order + 1; i ++)
		psum ->  coefficients [i] = psum -> coefficients [g_order + i];

	psum -> numberOfCoefficients = g_order + 1;
	/*
		To Chebychev
	*/
	cos2x (psum -> coefficients.part(1, psum -> numberOfCoefficients));
}

static void LPC_Frame_into_Polynomial_dif (LPC_Frame me, Polynomial pdif) {
	/*
		Fa (z) = A(z) - z^-(p+1)A(1/z)
	*/
	const integer order = my nCoefficients;
	pdif -> coefficients [order + 2] = -1.0;
	for (integer i = 1; i <= order; i ++)
		pdif -> coefficients [order + 2 - i] = - my a [i] + my a [order + 1 - i];

	pdif -> coefficients [1] = 1.0;
	pdif -> numberOfCoefficients = order + 2;

	if (order % 2 == 0) {
		/*
			Fa(z)/(z-1)
		*/
		Polynomial_divide_firstOrderFactor (pdif, 1.0, nullptr);
	} else {
		/*
			Fa(z) / (z^2 - 1)
		*/
		Polynomial_divide_secondOrderFactor (pdif, 1.0);
	}
	/*
		Transform to cos(w) terms
	*/
	integer g_order = pdif -> numberOfCoefficients / 2;
	for (integer i = 1; i <= g_order + 1; i ++)
		pdif -> coefficients [i] = pdif -> coefficients [g_order + i];

	pdif -> numberOfCoefficients = g_order + 1;
	/*
		To Chebychev
	*/
	cos2x (pdif -> coefficients.part(1, pdif -> numberOfCoefficients));
}

static integer Polynomial_into_Roots_searchOnGrid (Polynomial me, Roots thee, double gridSize) {
	Melder_assert (thy numberOfRoots >= my numberOfCoefficients - 1);
	double xmin = my xmin;
	integer numberOfRootsFound = 0;
	while (xmin < my xmax && numberOfRootsFound < my numberOfCoefficients - 1) {
		double xmax = xmin + gridSize;
		xmax = xmax > my xmax ? my xmax : xmax;
		const double root = Polynomial_findOneSimpleRealRoot_ridders (me, xmin, xmax);
		if (isdefined (root) && (numberOfRootsFound == 0 || thy roots [numberOfRootsFound].real() != root)) {
			thy roots [++ numberOfRootsFound]. real (root); // root not at border of interval
			thy roots [numberOfRootsFound]. imag (0.0);
		}
		xmin = xmax;
	}
	return numberOfRootsFound;
}

bool structLPCFrameIntoLineSpectralFrequenciesFrame :: inputFrameToOutputFrame () {
	LPC_Frame inputFrame = & inputLPC -> d_frames [currentFrame];
	LineSpectralFrequencies_Frame outputFrame = & outputLSF -> d_frames [currentFrame];
	Melder_assert (inputFrame -> nCoefficients == inputFrame -> a.size); // check invariant
	const double maximumFrequency = outputLSF -> maximumFrequency;
	/*
		Construct Fs and Fa
		divide out the zeros
		transform to polynomial equations gsum and gdif of half the order
	*/
	LPC_Frame_into_Polynomial_sum (inputFrame, gsum.get());
	const integer half_order_gsum = gsum -> numberOfCoefficients - 1;
	LPC_Frame_into_Polynomial_dif (inputFrame, gdif.get());
	const integer half_order_gdif = gdif -> numberOfCoefficients - 1;
	double currentGridSize = gridSize;
	integer numberOfBisections = 0, numberOfRootsFound = 0;
	while (numberOfRootsFound  < half_order_gsum && numberOfBisections < 10) {
		numberOfRootsFound = Polynomial_into_Roots_searchOnGrid (gsum.get(), roots.get(), currentGridSize);
		currentGridSize *= 0.5;
		numberOfBisections++;
	}
	frameAnalysisInfo = 0;
	if (numberOfBisections >= 10)
		frameAnalysisInfo = 1; // too many bisections
	/*
		[gsum-> xmin, gsum -> xmax] <==> [nyquistFrequency, 0],
		i.e. highest root corresponds to lowest frequency
	*/
	for (integer i = 1; i <= half_order_gsum; i ++)
		outputFrame -> frequencies [2 * i - 1] = acos (roots -> roots [half_order_gsum + 1 - i].real() / 2.0) / NUMpi * maximumFrequency;
	/*
		The roots of gdif lie inbetween the roots of gsum
	*/
	for (integer i = 1; i <= half_order_gdif; i ++) {
		const double xmax = roots -> roots [half_order_gsum + 1 - i].real();
		const double xmin = ( i == half_order_gsum ? gsum -> xmin : roots -> roots [half_order_gsum - i].real() );
		const double root = Polynomial_findOneSimpleRealRoot_ridders (gdif.get(), xmin, xmax);
		if (isdefined (root))
			outputFrame -> frequencies [2 * i] = acos (root / 2.0) / NUMpi * maximumFrequency;
		else
			outputFrame -> numberOfFrequencies --;
	}
	outputFrame -> frequencies.resize (outputFrame -> numberOfFrequencies); // maintain invariant
	return true;
}

#if 0
/* g [0]+g [1]x + ... g [m]*x^ m = 0 ; m should be even
 * Semenov, Kalyuzhny, Kovtonyuk (2003), Efficient calculation of line spectral frequencies based on new method for solution of transcendental equations,
 * ICASSP 2003, 457--460
 * 		g [0 .. g_order]
 * 		work [0.. g_order + 1 + (numberOfDerivatives + 1) * 5]
 * 		root [1 .. (g_order+1)/2]
 */
static void Roots_fromPolynomial (Roots me, Polynomial g, integer numberOfDerivatives, double *work) {
	if (numberOfDerivatives < 3) {
		Melder_throw (U"Number of derivatives should be at least 3.");
	}
	double xmin = -1.0, xmax = 1.0;
	integer numberOfRootsFound = 0;
	integer g_order = g -> numberOfCoefficients - 1;
	double *gabs = work, *fact = gabs + g_order + 1, *p2 = fact + numberOfDerivatives + 1;
	double *derivatives = p2 + numberOfDerivatives + 1, *constraints = derivatives + numberOfDerivatives + 1;
	double *intervals = constraints + numberOfDerivatives + 1;
	
	/* Fill vectors with j! and 2^j only once */
	fact [0] = p2 [0] = 1.0;
	for (integer j = 1; j <= numberOfDerivatives; j ++) {
		fact [j] = fact [j - 1] * j; // j!
		p2 [j] = p2 [j - 1] * 2.0; // 2^j
	}
	
	/* The constraints M [j] (Semenov et al. eq. (8)) can be calculated by taking absolute values of 
	 * the polynomial coefficients and evaluating the polynomial and the derivatives at x = 1.0
	 */
	for (integer k = 0; k <= g_order; k ++) {
		gabs [k] = fabs (g -> coefficients [k + 1]);
	}
	evaluatePolynomialAndDerivatives (gabs, g_order, 1.0, constraints, numberOfDerivatives);
	intervals [0] = 1.0;
	while (numberOfRootsFound < g_order || xmin == xmax) {
		double dsum1 = 0.0, dsum2 = 0.0;
		double xmid = (xmin + xmax) / 2.0;
		evaluatePolynomialAndDerivatives (g, g_order, xmid, derivatives, numberOfDerivatives);
		double fxmid = derivatives [0], fdxmin = derivatives [1];
		integer j = 1;
		bool rootsOnIntervalPossible_f = true, rootsOnIntervalPossible_df = true;
		while (j <= numberOfDerivatives && (rootsOnIntervalPossible_f || rootsOnIntervalPossible_df)) {
			intervals [j] = intervals [j - 1] * (xmax - xmin);
			integer k = j - 1;
			if (j > 1) {   // start at first derivative
				dsum1 += fabs (derivatives [k]) * intervals [k] / (p2 [k] * fact [k]);
			}
			if (j > 2) {   // start at second derivative
				dsum2 += fabs (derivatives [k]) * intervals [k - 1] / (p2 [k - 1] * fact [k - 1]);
				if (rootsOnIntervalPossible_f) {
					double testValue1 = dsum1 + constraints [j] * intervals [j] / (p2 [j] * fact [j]);
					rootsOnIntervalPossible_f = ! (fxmid + testValue1 < 0.0 || fxmid - testValue1 > 0.0);
				}
				if (rootsOnIntervalPossible_df) {
					double testValue2 = dsum2 + constraints [j] * intervals [j - 1] / (p2 [j - 1] * fact [j - 1]);
					rootsOnIntervalPossible_df = ! (fdxmin + testValue2 < 0.0 || fdxmin - testValue2 > 0.0);
				}
			}
			j++;
		}
		if (rootsOnIntervalPossible_f) {
			if (rootsOnIntervalPossible_df) {   // f(x) uncertain && f'(x) uncertain: bisect
				xmax = xmid;
			} else {   // f(x) uncertain; f'(x) certain
				double fxmin = evaluatePolynomial (g, g_order, xmin);
				double fxmax = evaluatePolynomial (g, g_order, xmax);
				if (fxmin * fxmax <= 0.0) {
					double root;
					NUMnrbis (dpoly, xmin, xmax, &poly, &root);
					roots [++numberOfRootsFound] = root;
				} else {
					xmin = xmax; xmax = 1.0;
				}
			}
		} else {
			xmin = xmax; xmax = 1.0;
		}
	}	
}
#endif

void LPCFrameIntoLineSpectralFrequenciesFrame_init (mutableLPCFrameIntoLineSpectralFrequenciesFrame me, constLPC input,
	mutableLineSpectralFrequencies output)
{
	SampledFrameIntoSampledFrame_init (me, output);
	LPCFrameIntoSampledFrame_init (me, input);
	my outputLSF = output;
	const integer numberOfCoefficients = input -> maxnCoefficients + 1;
	my gsum = Polynomial_create (-2.0, 2.0, numberOfCoefficients); // large enough
	my gdif = Polynomial_create (-2.0, 2.0, numberOfCoefficients);
	my roots = Roots_create (numberOfCoefficients / 2);
}

autoLPCFrameIntoLineSpectralFrequenciesFrame LPCFrameIntoLineSpectralFrequenciesFrame_create (constLPC input, 
	mutableLineSpectralFrequencies output)
{
	try {
		autoLPCFrameIntoLineSpectralFrequenciesFrame me = Thing_new (LPCFrameIntoLineSpectralFrequenciesFrame);
		LPCFrameIntoLineSpectralFrequenciesFrame_init (me.get(), input, output);
		return me;
		// gridSize not initialized in here
	} catch (MelderError) {
		Melder_throw (U"LPCFrameIntoLineSpectralFrequenciesFrame not created.");
	}
}

/****************** LineSpectralFrequencies to LPC *********************************/

Thing_implement (LineSpectralFrequenciesFrameIntoLPCFrame, SampledFrameIntoSampledFrame, 0);

void structLineSpectralFrequenciesFrameIntoLPCFrame :: allocateOutputFrames () {
	for (integer iframe = 1; iframe <= outputLPC -> nx; iframe ++) {
		LPC_Frame lpcFrame = & outputLPC -> d_frames [iframe];
		LineSpectralFrequencies_Frame lsfFrame = & inputLSF -> d_frames [iframe];
		LPC_Frame_init (lpcFrame, lsfFrame -> numberOfFrequencies);
	}
}


bool structLineSpectralFrequenciesFrameIntoLPCFrame :: inputFrameToOutputFrame () {
	LPC_Frame lpcFrame = & outputLPC -> d_frames [currentFrame];
	VEC a = lpcFrame -> a.get();
	LineSpectralFrequencies_Frame lsfFrame = & inputLSF -> d_frames [currentFrame];
	const double maximumFrequency = inputLSF -> maximumFrequency;
	const integer numberOfFrequencies = lsfFrame -> numberOfFrequencies;
	integer numberOfOmegas = (numberOfFrequencies + 1) / 2;
	/*
		Reconstruct Fs (z)
		Use lpcFrame -> a as a buffer whose size changes!!!
	*/
	for (integer i = 1; i <= numberOfOmegas; i ++) {
		const double omega = lsfFrame -> frequencies [2 * i - 1] / maximumFrequency * NUMpi;
		a [i] = -2.0 * cos (omega);
	}
	Polynomial_initFromProductOfSecondOrderTerms (fs.get(), a.part (1, numberOfOmegas));
	/*
		Reconstruct Fa (z)
	*/
	numberOfOmegas = numberOfFrequencies / 2;
	for (integer i = 1; i <= numberOfOmegas; i ++) {
		const double omega = lsfFrame -> frequencies [2 * i] / maximumFrequency * NUMpi;
		a [i] = -2.0 * cos (omega);
	}
	Polynomial_initFromProductOfSecondOrderTerms (fa.get(), a.part (1, numberOfOmegas));
	
	if (numberOfFrequencies % 2 == 0) {
		Polynomial_multiply_firstOrderFactor (fs.get(), -1.0);   // * (z + 1)
		Polynomial_multiply_firstOrderFactor (fa.get(), 1.0);   // * (z - 1)
	} else {
		Polynomial_multiply_secondOrderFactor (fa.get(), 1.0);   // * (z^2 - 1)
	}
	Melder_assert (fs -> numberOfCoefficients == fa -> numberOfCoefficients);
	/*
		A(z) = (Fs(z) + Fa(z) / 2
	*/
	for (integer i = 1; i <= fs -> numberOfCoefficients - 2; i ++)
		a [lsfFrame -> numberOfFrequencies - i + 1] = 0.5 * (fs -> coefficients [i + 1] + fa -> coefficients [i + 1]);
	frameAnalysisInfo = 0;
	return true;
}

void LineSpectralFrequenciesFrameIntoLPCFrame_init (
	mutableLineSpectralFrequenciesFrameIntoLPCFrame me, constLineSpectralFrequencies inputLSF, mutableLPC outputLPC) 
{
	SampledFrameIntoSampledFrame_init (me, outputLPC);
	my inputLSF = inputLSF;
	my outputLPC = outputLPC;
	my fs = Polynomial_create (-1.0, 1.0, inputLSF -> maximumNumberOfFrequencies + 2);
	my fa = Polynomial_create (-1.0, 1.0, inputLSF -> maximumNumberOfFrequencies + 2);
}

autoLineSpectralFrequenciesFrameIntoLPCFrame LineSpectralFrequenciesFrameIntoLPCFrame_create 
	(constLineSpectralFrequencies inputLSF, mutableLPC outputLPC)
{
	try {
		autoLineSpectralFrequenciesFrameIntoLPCFrame me = Thing_new (LineSpectralFrequenciesFrameIntoLPCFrame);
		LineSpectralFrequenciesFrameIntoLPCFrame_init (me.get(), inputLSF, outputLPC);
		return me;
	} catch (MelderError) {
		Melder_throw (inputLSF, U" no LineSpectralFrequenciesFrameIntoLPCFrame created.");
	}	
}

/* End of file LPCFrameAndLineSpectralFrequenciesFrame.cpp */
