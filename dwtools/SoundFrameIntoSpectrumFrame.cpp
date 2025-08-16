/* SoundFrameIntoSpectrumFrame.cpp
 *
 * Copyright (C) 2025 David Weenink
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

#include "SoundFrameIntoSpectrumFrame.h"
#include "NUM2.h"

#include "oo_DESTROY.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_COPY.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_EQUAL.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_CAN_WRITE_AS_ENCODING.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_WRITE_TEXT.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_WRITE_BINARY.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_READ_TEXT.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_READ_BINARY.h"
#include "SoundFrameIntoSpectrumFrame_def.h"
#include "oo_DESCRIPTION.h"
#include "SoundFrameIntoSpectrumFrame_def.h"

Thing_implement (SoundFrameIntoSpectrumFrame, SoundFrameIntoSampledFrame, 0);


void structSoundFrameIntoSpectrumFrame :: inputFrameIntoOutputFrame () {}

void structSoundFrameIntoSpectrumFrame :: Sound_into_Spectrum () {

	if (numberOfChannels == 1)
		fftData.part (1, soundFrameSize)  <<=  my z.row (1);
	else {
		/*
			Multiple channels: take the average.
			*/
		for (integer ichan = 1; ichan <= numberOfChannels; ichan ++)
			fftData.part (1, soundFrameSize)  +=  my z.row (ichan);
		fftData.part (1, soundFrameSize)  *=  1.0 / numberOfChannels;
	}
	fftData.part (soundFrameSize + 1, numberOfFourierSamples)  <<=  0.0;
	NUMfft_forward (fourierTable.get(), fftData.get());

	const VEC re = spectrum -> z.row (1);
	const VEC im = spectrum -> z.row (2);
	const double scaling = output -> dx;
	re [1] = fftData [1] * scaling;
	im [1] = 0.0;
	for (integer i = 2; i < numberOfFrequencies; i ++) {
		re [i] = fftData [i + i - 2] * scaling;   // fftData [2], fftData [4], ...
		im [i] = fftData [i + i - 1] * scaling;   // fftData [3], fftData [5], ...
	}
	if ((numberOfFourierSamples & 1) != 0) {
		if (numberOfFourierSamples > 1) {
			re [numberOfFrequencies] = fftData [numberOfFourierSamples - 1] * scaling;
			im [numberOfFrequencies] = fftData [numberOfFourierSamples] * scaling;
		}
	} else {
		re [numberOfFrequencies] = fftData [numberOfFourierSamples] * scaling;
		im [numberOfFrequencies] = 0.0;
	}
}

void SoundFrameIntoSpectrumFrame_init (mutableSoundFrameIntoSpectrumFrame me, constSound input, mutableSampled output,
	double effectiveAnalysisWidth, kSound_windowShape windowShape, bool fast)
{
	SoundFrameIntoSampledFrame_init (me, input, output, effectiveAnalysisWidth, windowShape);
	my fast = fast;
	my numberOfFourierSamples = ( fast ? Melder_iroundUpToPowerOfTwo (my frameAsSound -> nx) : my frameAsSound -> nx );
	my numberOfFrequencies = numberOfFourierSamples / 2 + 1;
	my fourierTable = NUMFourierTable_create (numberOfFourierSamples);
	my fftData = raw_VEC (numberOfFourierSamples);
	my spectrum = Spectrum_create (0.5 / my frameAsSound -> dx, my numberOfFrequencies);
	my spectrum -> dx = 1.0 / (my frameAsSound -> dx * my numberOfFourierSamples);
}

/* End of file SoundFrameIntoSpectrumFrame.cpp */
