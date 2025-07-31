/* LPC_and_LineSpectralFrequencies.cpp
 *
 * Copyright (C) 2016-2020, 2025 David Weenink
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 djmw 20160421  Initial version
*/

#include "LPC_and_LineSpectralFrequencies.h"
#include "LPCFrameAndLineSpectralFrequenciesFrame.h"
#include "NUM2.h"
#include "Polynomial.h"
#include "Roots.h"
#include "SampledIntoSampled.h"
#include "SampledIntoSampledStatus.h"

void LPC_into_LineSpectralFrequencies (constLPC me, mutableLineSpectralFrequencies outputLSF, double gridSize) {
	SampledIntoSampled_requireEqualDomainsAndSampling (me, outputLSF);
	autoLPCFrameIntoLineSpectralFrequenciesFrame ws = LPCFrameIntoLineSpectralFrequenciesFrame_create (me, outputLSF);
	if (gridSize <= 0.0)
		gridSize = 0.02;
	ws -> gridSize = gridSize;
	autoLPCAndLineSpectralFrequenciesStatus status = LPCAndLineSpectralFrequenciesStatus_create (outputLSF -> nx);
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, outputLSF, ws.move(), status.move());
	const integer numberOfErrorFrames = SampledIntoSampled_analyseThreaded (sis.get());
	if (numberOfErrorFrames > 0)
		Melder_warning (U"LineSpectralFrequencies_into_LPC: ", numberOfErrorFrames, U" frames have issues.");
}

autoLineSpectralFrequencies LPC_to_LineSpectralFrequencies (constLPC me, double gridSize) {
	try {
		const double nyquistFrequency = 0.5 / my samplingPeriod;
		autoLineSpectralFrequencies outputLSF = LineSpectralFrequencies_create (my xmin, my xmax, my nx, my dx, my x1, my maxnCoefficients, nyquistFrequency);
		LPC_into_LineSpectralFrequencies (me, outputLSF.get(), gridSize);
		return outputLSF;
	} catch (MelderError) {
		Melder_throw (me, U": no LineSpectralFrequencies created.");
	}
}

/**************************** LineSpectralFrequencies to LPC **************************************************/

void LineSpectralFrequencies_into_LPC (constLineSpectralFrequencies me, mutableLPC outputLPC) {
	SampledIntoSampled_requireEqualDomainsAndSampling (me, outputLPC);
	autoLineSpectralFrequenciesFrameIntoLPCFrame ws = LineSpectralFrequenciesFrameIntoLPCFrame_create (me, outputLPC);
	autoLPCAndLineSpectralFrequenciesStatus status = LPCAndLineSpectralFrequenciesStatus_create (outputLPC -> nx);
	autoSampledIntoSampled sis = SampledIntoSampled_create (me, outputLPC, ws.move(), status.move());
	const integer numberOfErrorFrames = SampledIntoSampled_analyseThreaded (sis.get());
	if (numberOfErrorFrames > 0)
		Melder_warning (U"LineSpectralFrequencies_into_LPC: ", numberOfErrorFrames, U" frames have issues.");
}

autoLPC LineSpectralFrequencies_to_LPC (constLineSpectralFrequencies me) {
	try {
		autoLPC thee = LPC_create (my xmin, my xmax, my nx, my dx, my x1, my maximumNumberOfFrequencies, 0.5 / my maximumFrequency);
		LineSpectralFrequencies_into_LPC (me, thee.get());	
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": no LPC created from LineSpectralFrequencies.");
	}
}

/* End of file LPC_and_LineSpectralFrequencies.cpp */
