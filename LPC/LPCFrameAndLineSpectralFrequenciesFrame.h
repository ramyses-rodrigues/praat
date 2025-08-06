#ifndef _LPCFrameAndLineSpectralFrequenciesFrame_h_
#define _LPCFrameAndLineSpectralFrequenciesFrame_h_
/* LPCFrameAndLineSpectralFrequenciesFrame.h
 *
 * Copyright (C) 2024-2025 David Weenink
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

#include "LPC.h"
#include "LineSpectralFrequencies.h"
#include "Polynomial.h"
#include "Roots.h"
#include "LPCFrameIntoSampledFrame.h"

#include "LPCFrameAndLineSpectralFrequenciesFrame_def.h"

void LPCFrameIntoLineSpectralFrequenciesFrame_init (mutableLPCFrameIntoLineSpectralFrequenciesFrame me, constLPC input,
	mutableLineSpectralFrequencies output);

autoLPCFrameIntoLineSpectralFrequenciesFrame LPCFrameIntoLineSpectralFrequenciesFrame_create (constLPC input, 
	mutableLineSpectralFrequencies output);

void LineSpectralFrequenciesFrameIntoLPCFrame_init
	(mutableLineSpectralFrequenciesFrameIntoLPCFrame me, constLineSpectralFrequencies inputLSF, mutableLPC outputLPC);

autoLineSpectralFrequenciesFrameIntoLPCFrame LineSpectralFrequenciesFrameIntoLPCFrame_create (
	constLineSpectralFrequencies inputLSF, mutableLPC outputLPC	
);
	
#endif /*_LPCFrameAndLineSpectralFrequenciesFrame_h_ */
