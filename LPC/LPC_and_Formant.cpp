/* LPC_and_Formant.cpp
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

#include "LPC_and_Formant.h"
#include "LPC_and_Polynomial.h"
#include "NUM2.h"
#include "Roots_and_Formant.h"
#include "SampledIntoSampled2.h"

static void LPC_Frame_into_Polynomial (constLPC_Frame me, mutablePolynomial p) {
	/*
		The lpc coefficients are a[1..nCoefficients], a[0] == 1, is not stored
		For the polynomial we therefore need one extra coefficient since
		the a's are stored in reverse order in the polynomial and a[0]
		represents the highest power and is stored into the last position
		of the polynomial.
	*/
	Melder_assert (my nCoefficients  == my a.size); // check invariant
	const integer highestPolynomialCoefficientNumber = my nCoefficients + 1;
	Melder_assert (p -> _capacity >= highestPolynomialCoefficientNumber);
	
	p -> coefficients.resize (highestPolynomialCoefficientNumber);
	p -> numberOfCoefficients = p -> coefficients.size; // maintain invariant
	p -> coefficients [highestPolynomialCoefficientNumber] = 1.0;
	for (integer icof = 1; icof <= my nCoefficients; icof ++)
		p -> coefficients [icof] = my a [highestPolynomialCoefficientNumber - icof];
}

void Formant_Frame_init (Formant_Frame me, integer numberOfFormants) {
	if (numberOfFormants > 0)
		my formant = newvectorzero <structFormant_Formant> (numberOfFormants);
	my numberOfFormants = my formant.size; // maintain invariant
}

Thing_implement (LPCFrameIntoFormantFrame2, SampledFrameIntoSampledFrame2, 0);

void structLPCFrameIntoFormantFrame2 :: initBasicLPCFrameAndFormant (constLPC inputLPC, mutableFormant outputFormant, double margin) {
	LPCFrameIntoFormantFrame2_Parent :: initBasic (inputLPC, outputFormant);
	our inputLPC = inputLPC;
	our outputFormant = outputFormant;
	our margin = margin;
	our order = inputLPC -> maxnCoefficients;
}

void structLPCFrameIntoFormantFrame2 :: initHeap () {
	LPCFrameIntoFormantFrame2_Parent :: initHeap ();
	bufferSize = order * order + order + order + 11 * order;
	buffer = raw_VEC (bufferSize);		
	p = Polynomial_create (-1.0, 1.0, order);
	roots = Roots_create (order);
}

bool structLPCFrameIntoFormantFrame2 :: inputFrameIntoOutputFrame (integer iframe) {
	Formant_Frame formantFrame = & outputFormant -> frames [iframe];
	LPC_Frame inputLPCFrame = & inputLPC -> d_frames [iframe];
	formantFrame -> intensity = inputLPCFrame -> gain;
	integer frameAnalysisInfo = 0;
	if (inputLPCFrame -> nCoefficients == 0) {
		formantFrame -> numberOfFormants = 0;
		formantFrame -> formant.resize (formantFrame -> numberOfFormants); // maintain invariant
		frameAnalysisInfo = 1;	
		return true;
	}
	frameAnalysisInfo = 0;
	const double samplingFrequency = 1.0 / inputLPC -> samplingPeriod;
	LPC_Frame_into_Polynomial (inputLPCFrame, p.get());
	Polynomial_into_Roots (p.get(), roots.get(), buffer.get());
	Roots_fixIntoUnitCircle (roots.get());
	Roots_into_Formant_Frame (roots.get(), formantFrame, samplingFrequency, margin);
	return true;
}

autoLPCFrameIntoFormantFrame2 LPCFrameIntoFormantFrame2_create (constLPC inputLPC, mutableFormant outputFormant, double margin) {
	try {
		autoLPCFrameIntoFormantFrame2 me = Thing_new (LPCFrameIntoFormantFrame2);
		my initBasicLPCFrameAndFormant (inputLPC, outputFormant, margin);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Cannot create LPCFrameIntoFormantFrame.");
	}
}

void LPC_into_Formant (constLPC me, mutableFormant thee, double margin) {
	SampledIntoSampled2_requireEqualDomainsAndSampling (me, thee);
	autoLPCFrameIntoFormantFrame2 frameIntoFrame = LPCFrameIntoFormantFrame2_create (me, thee, margin);
	SampledIntoSampled_mt (frameIntoFrame.get(), 40);
	Formant_sort (thee);
}

autoFormant LPC_to_Formant (constLPC me, double margin) {
	try {
		/*
			In very extreme case all roots of the lpc polynomial might be real.
			A real root gives either a frequency at 0 Hz or at the Nyquist frequency.
			If margin > 0 these frequencies are filtered out and the number of formants can never exceed
			(my maxnCoefficients+1) / 2.
		*/
		const integer maximumNumberOfFormants = ( margin == 0.0 ? my maxnCoefficients : (my maxnCoefficients + 1) / 2 );
		Melder_require (my maxnCoefficients < 100,
			U"We cannot find the roots of a polynomial of order > 99.");
		autoFormant thee = Formant_create (my xmin, my xmax, my nx, my dx, my x1, maximumNumberOfFormants);
		for (integer iframe = 1; iframe <= thy nx; iframe ++) {
			Formant_Frame_init (& thy frames [iframe], maximumNumberOfFormants);
		}
		LPC_into_Formant (me, thee.get(), margin);
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no Formant created.");
	}
}

void Formant_Frame_scale (Formant_Frame me, double scale) {
	for (integer iformant = 1; iformant <= my numberOfFormants; iformant ++) {
		my formant [iformant]. frequency *= scale;
		my formant [iformant]. bandwidth *= scale;
	}
}

void LPC_Frame_into_Formant_Frame (constLPC_Frame me, Formant_Frame thee, double samplingPeriod, double margin) {
	Melder_assert (my nCoefficients == my a.size); // check invariant
	thy intensity = my gain;
	if (my nCoefficients == 0) {
		thy formant.resize (0);
		thy numberOfFormants = thy formant.size; // maintain invariant
		return;
	}
	autoPolynomial p = LPC_Frame_to_Polynomial (me);
	autoRoots r = Polynomial_to_Roots (p.get());
	Roots_fixIntoUnitCircle (r.get());
	Roots_into_Formant_Frame (r.get(), thee, 1.0 / samplingPeriod, margin);
}

void Formant_Frame_into_LPC_Frame (constFormant_Frame me, LPC_Frame thee, double samplingPeriod) {
	if (my numberOfFormants < 1)
		return;
	const double nyquistFrequency = 0.5 / samplingPeriod;
	integer numberOfPoles = 2 * my numberOfFormants;
	autoVEC lpc = zero_VEC (numberOfPoles + 2);   // all odd coefficients have to be initialized to zero
	lpc [2] = 1.0;
	integer m = 2;
	for (integer iformant = 1; iformant <= my numberOfFormants; iformant ++) {
		const double formantFrequency = my formant [iformant]. frequency;
		if (formantFrequency > nyquistFrequency)
			continue;
		/*
			D(z): 1 + p z^-1 + q z^-2
		*/
		const double r = exp (- NUMpi * my formant [iformant]. bandwidth * samplingPeriod);
		const double p = - 2.0 * r * cos (NUM2pi * formantFrequency * samplingPeriod);
		const double q = r * r;
		/*
			By setting the two extra elements (0, 1) in the lpc vector we can avoid boundary testing;
			lpc [3..n+2] come to contain the coefficients.
		*/
		for (integer j = m + 2; j > 2; j --)
			lpc [j] += p * lpc [j - 1] + q * lpc [j - 2];
		m += 2;
	}
	if (thy nCoefficients < numberOfPoles)
		numberOfPoles = thy nCoefficients;
	for (integer i = 1; i <= numberOfPoles; i ++)
		thy a [i] = lpc [i + 2];
	thy gain = my intensity;
}

autoLPC Formant_to_LPC (constFormant me, double samplingPeriod) {
	try {
		autoLPC thee = LPC_create (my xmin, my xmax, my nx, my dx, my x1, 2 * my maxnFormants, samplingPeriod);

		for (integer i = 1; i <= my nx; i ++) {
			const Formant_Frame f = & my frames [i];
			const LPC_Frame lpcFrame = & thy d_frames [i];
			const integer numberOfCoefficients = 2 * f -> numberOfFormants;
			LPC_Frame_init (lpcFrame, numberOfCoefficients);
			Formant_Frame_into_LPC_Frame (f, lpcFrame, samplingPeriod);
		}
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC created.");
	}
}

/* End of file LPC_and_Formant.cpp */
