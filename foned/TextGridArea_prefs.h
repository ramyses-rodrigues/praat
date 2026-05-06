/* TextGridArea_prefs.h
 *
 * Copyright (C) 2013,2015,2016,2019,2022 Paul Boersma
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

Prefs_begin (TextGridArea)

	ClassPrefs_overrideBool       (TextGridArea, picture_garnish,          1, true)

	InstancePrefs_addBool         (TextGridArea, useTextStyles,            1, false)
	InstancePrefs_addDouble       (TextGridArea, fontSize,                 1, U"18")
	InstancePrefs_addEnum         (TextGridArea, alignment,                1, kGraphics_horizontalAlignment, DEFAULT)
	InstancePrefs_addBool         (TextGridArea, shiftDragMultiple,        1, true)
	InstancePrefs_addEnum         (TextGridArea, showNumberOf,             1, kTextGridArea_showNumberOf, DEFAULT)
	InstancePrefs_addEnum         (TextGridArea, greenMethod,              1, kMelder_string, DEFAULT)
	InstancePrefs_addString       (TextGridArea, greenString,              1, U"some text here for green paint")
	ClassPrefs_addBool            (TextGridArea, picture_showBoundaries,   1, true)
	ClassPrefs_addBool            (TextGridArea, picture_pitch_speckle,    1, false)
	InstancePrefs_addString       (TextGridArea, align_language,           1, U"English")
	InstancePrefs_addBool         (TextGridArea, align_includeWords,       1, true)
	InstancePrefs_addBool         (TextGridArea, align_includePhonemes,    1, false)
	InstancePrefs_addBool         (TextGridArea, align_allowSilences,      1, false)
	InstancePrefs_addString       (TextGridArea, transcribe_model,         1, theSpeechRecognizerDefaultModelName)
	InstancePrefs_addString       (TextGridArea, transcribe_language,      1, theSpeechRecognizerDefaultLanguageName)
	InstancePrefs_addBool		  (TextGridArea, transcribe_includeWords,  1, theSpeechRecognizerDefaultIncludeWords)
	InstancePrefs_addBool		  (TextGridArea, transcribe_includeDiarization, 1, theSpeechRecognizerDefaultIncludeDiarization)
	InstancePrefs_addBool		  (TextGridArea, transcribe_useVad,		   1, theSpeechRecognizerDefaultUseVad)
	InstancePrefs_addDouble       (TextGridArea, transcribe_vadThreshold,  1, theVadDefaultThresholdStr)
	InstancePrefs_addDouble       (TextGridArea, transcribe_vadMinNonSpeech, 1, theVadDefaultMinNonSpeechDurationStr)
	InstancePrefs_addDouble       (TextGridArea, transcribe_vadMinSpeech,  1, theVadDefaultMinSpeechDurationStr)
	InstancePrefs_addDouble       (TextGridArea, transcribe_vadPadding,    1, theVadDefaultSpeechPadStr)
	InstancePrefs_addInteger      (TextGridArea, diarize_maxSimultaneiosSpeakers, 1, theDiarizationMaxSimultaneousSpeakersStr)
	InstancePrefs_addInteger      (TextGridArea, diarize_numSpeakers,	   1, theDiarizationNumSpeakersStr)
	InstancePrefs_addInteger      (TextGridArea, diarize_maxSpeakers,	   1, theDiarizationMaxSpeakersStr)
	InstancePrefs_addInteger      (TextGridArea, diarize_minSpeakers,	   1, theDiarizationMinSpeakersStr)
	InstancePrefs_addDouble       (TextGridArea, diarize_clusterThreshold, 1, theDiarizationClusterThresholdStr)
	InstancePrefs_addInteger      (TextGridArea, diarize_segmentationOverlap, 1, theDiarizationSegmentationOverlapStr)
	InstancePrefs_addString		  (TextGridArea, diarize_nonSpeechLabel,   1, theDiarizationDefaultNonSpeechLabel)
	InstancePrefs_addString		  (TextGridArea, diarize_speechLabel,      1, theDiarizationDefaultSpeechLabel)

Prefs_end (TextGridArea)

/* End of file TextGridArea_prefs.h */
