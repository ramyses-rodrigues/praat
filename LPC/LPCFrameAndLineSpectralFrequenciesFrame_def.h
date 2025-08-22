/* LPCFrameAndLineSpectralFrequenciesFrame_def.h
 *
 * Copyright (C) 2016-2020 David Weenink
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#define ooSTRUCT LPCFrameIntoLineSpectralFrequenciesFrame
oo_DEFINE_CLASS (LPCFrameIntoLineSpectralFrequenciesFrame, LPCFrameIntoSampledFrame)

	oo_UNSAFE_BORROWED_TRANSIENT_MUTABLE_OBJECT_REFERENCE (LineSpectralFrequencies, outputLSF)
	oo_INTEGER (maximumNumberOfFrequencies)
	oo_OBJECT (Polynomial, 0, gsum)
	oo_OBJECT (Polynomial, 0, gdif)
	oo_OBJECT (Roots, 0, roots)
	oo_DOUBLE (gridSize)

	#if oo_DECLARING
		bool inputFrameToOutputFrame ()
			override;   // the analysis, also sets 'frameAnalysisInfo'
		void allocateOutputFrames ()
			override;
	#endif

oo_END_CLASS (LPCFrameIntoLineSpectralFrequenciesFrame)
#undef ooSTRUCT

#define ooSTRUCT LineSpectralFrequenciesFrameIntoLPCFrame
oo_DEFINE_CLASS (LineSpectralFrequenciesFrameIntoLPCFrame, SampledFrameIntoSampledFrame)

	oo_UNSAFE_BORROWED_TRANSIENT_CONST_OBJECT_REFERENCE (LineSpectralFrequencies, inputLSF)
	oo_UNSAFE_BORROWED_TRANSIENT_MUTABLE_OBJECT_REFERENCE (LPC, outputLPC)
	oo_OBJECT (Polynomial, 0, fs)
	oo_OBJECT (Polynomial, 0, fa)

	#if oo_DECLARING
		bool inputFrameToOutputFrame ()
			override;   // the analysis, sets 'frameAnalysisInfo'
		void allocateOutputFrames ()
			override;
	#endif

oo_END_CLASS (LineSpectralFrequenciesFrameIntoLPCFrame)
#undef ooSTRUCT


/* Endof file LPCFrameAndLineSpectralFrequenciesFrame_def.h */
