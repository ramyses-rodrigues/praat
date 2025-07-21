/* LPCFrameIntoLineSpectralFrequenciesFrame.cpp
 *
 * Copyright (C) 2024-2025 David Weenink
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

#include "melder.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame.h"
#include "Roots_and_Formant.h"

#include "oo_DESTROY.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_COPY.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_EQUAL.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_CAN_WRITE_AS_ENCODING.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_WRITE_TEXT.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_WRITE_BINARY.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_READ_TEXT.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_READ_BINARY.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"
#include "oo_DESCRIPTION.h"
#include "LPCFrameIntoLineSpectralFrequenciesFrame_def.h"

Thing_implement (LPCFrameIntoLineSpectralFrequenciesFrame, LPCFrameIntoSampledFrame, 0);

void structLPCFrameIntoLineSpectralFrequenciesFrame :: allocateOutputFrames () {
	for (integer iframe = 1; iframe <= output -> nx; iframe ++) {
		LineSpectralFrequencies_Frame lsff = & lineSpectralFrequencies -> d_frames [iframe];
		LineSpectralFrequencies_Frame_init (lsff, lineSpectralFrequencies -> maximumNumberOfFrequencies);
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

void structLPCFrameIntoLineSpectralFrequenciesFrame :: LPCFrameIntoLineSpectralFrequenciesFrame () {
	LPC_Frame inputFrame = & inputlpc -> d_frames [currentFrame];
	LineSpectralFrequencies_Frame outputFrame = & lineSpectralFrequencies -> d_frames [currentFrame];
	Melder_assert (inputFrame -> nCoefficients == inputFrame -> a.size); // check invariant
	const maximumFrequency = output -> 
	/*
		Construct Fs and Fa
		divide out the zeros
		transform to polynomial equations gsum and gdif of half the order
	*/
	LPC_Frame_into_Polynomial_sum (inputFrame, gsum.get());
	const integer half_order_gsum = gsum -> numberOfCoefficients - 1;
	LPC_Frame_into_Polynomial_dif (inputFrame, gdif.get());
	const integer half_order_gdif = gdif -> numberOfCoefficients - 1;
	
	integer numberOfBisections = 0, numberOfRootsFound = 0;
	while (numberOfRootsFound  < half_order_gsum && numberOfBisections < 10) {
		numberOfRootsFound = Polynomial_into_Roots_searchOnGrid (gsum.get(), roots.get(), gridSize);
		gridSize *= 0.5;
		numberOfBisections++;
	}
	
	Melder_require (numberOfBisections < 10,
		U"Too many bisections.");
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
}

void LPCFrameIntoLineSpectralFrequenciesFrame_init (mutableLPCFrameIntoLineSpectralFrequenciesFrame me, constLPC input,
	mutableLineSpectralFrequencies output)
{
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



/* End of file LPCFrameIntoLineSpectralFrequenciesFrame.cpp */

