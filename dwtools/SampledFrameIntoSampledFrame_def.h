/* SampledFrameIntoSampledFrame_def.h
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
	
/*
	A separate deep copy of the SampledIntoSampled is needed for each thread
*/

#define ooSTRUCT SampledFrameIntoSampledFrame
oo_DEFINE_CLASS (SampledFrameIntoSampledFrame, Daata)
	/*
		We need a reference to the output object because we need to know the number of frames to analyse
		which equals output -> nx;
	*/
	oo_UNSAFE_BORROWED_TRANSIENT_MUTABLE_OBJECT_REFERENCE (Sampled, output)
	/*
		if (updateStatus): we need the reference to the status object.
	*/
	oo_UNSAFE_BORROWED_TRANSIENT_MUTABLE_OBJECT_REFERENCE (SampledIntoSampledStatus, status)
	oo_INTEGER (startFrame)
	oo_INTEGER (currentFrame)				// the frame we are working on
	oo_INTEGER (frameAnalysisInfo)			// signals different "error" conditions etc in a frame analysis
	oo_BOOLEAN (updateStatus)				// gather extensive status and error info?
	oo_INTEGER (framesErrorCount)
	oo_DOUBLE (tol1)

	#if oo_DECLARING
		virtual void getInputFrame ();
		virtual bool inputFrameToOutputFrame ();   // the analysis, also sets 'frameAnalysisInfo'
		virtual void saveOutputFrame ();
		virtual void allocateOutputFrames ();
	#endif

oo_END_CLASS (SampledFrameIntoSampledFrame)
#undef ooSTRUCT

/* End of file SampledFrameIntoSampledFrame_def.h */
 
