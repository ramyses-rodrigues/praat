#ifndef _SpeechRecognizer_h_
#define _SpeechRecognizer_h_
/* SpeechRecognizer.h
 *
 * Copyright (C) 2025,2026 Anastasia Shchupak
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#include "Sound.h"
#include <set>

struct structSpeechRecognizer;
struct whisper_context;
struct whisper_vad_context;
struct whisper_vad_segments;
struct diarize_context;

/*
	Default Whisper model parameters.
*/
inline constexpr conststring32 theSpeechRecognizerDefaultModelName = U"ggml-base.bin";
inline constexpr conststring32 theSpeechRecognizerDefaultLanguageName = U"Autodetect language";

/*
	Default Silero-VAD parameters.
*/
inline constexpr double theVadDefaultThreshold = 0.5;
inline constexpr conststring32 theVadDefaultThresholdStr = U"0.5";   // for UI
inline constexpr double theVadDefaultMinNonSpeechDuration = 0.1;
inline constexpr conststring32 theVadDefaultMinNonSpeechDurationStr = U"0.1";   // for UI
inline constexpr double theVadDefaultMinSpeechDuration = 0.25;
inline constexpr conststring32 theVadDefaultMinSpeechDurationStr = U"0.25";   // for UI
inline constexpr double theVadDefaultSpeechPad = 0.03;
inline constexpr conststring32 theVadDefaultSpeechPadStr = U"0.03";   // for UI
inline constexpr conststring32 theVadDefaultSpeechLabel = U"speech";   // for UI
inline constexpr conststring32 theVadDefaultNonSpeechLabel = U"non-speech";   // for UI

/*
	Default transcription parameters.
*/
inline constexpr bool theSpeechRecognizerDefaultIncludeWords = true;
inline constexpr bool theSpeechRecognizerDefaultIncludeDiarization = false;
inline constexpr bool theSpeechRecognizerDefaultUseVad = true;

/*
	Default diarization parameters.
*/
inline constexpr integer theDiarizationMaxSimultaneousSpeakers = 3;
inline constexpr conststring32 theDiarizationMaxSimultaneousSpeakersStr = U"3";   // for UI
inline constexpr integer theDiarizationNumSpeakers = 0;
inline constexpr conststring32 theDiarizationNumSpeakersStr = U"0";   // for UI
inline constexpr integer theDiarizationMaxSpeakers = 0;
inline constexpr conststring32 theDiarizationMaxSpeakersStr = U"0";   // for UI
inline constexpr integer theDiarizationMinSpeakers = 0;
inline constexpr conststring32 theDiarizationMinSpeakersStr = U"0";   // for UI
inline constexpr double theDiarizationClusterThreshold = 0.7045654963945799;
inline constexpr conststring32 theDiarizationClusterThresholdStr = U"0.7045654963945799";   // for UI
inline constexpr integer theDiarizationSegmentationOverlap = 90;
inline constexpr conststring32 theDiarizationSegmentationOverlapStr = U"90";   // for UI


struct autoWhisperContext {
	whisper_context *ptr;

	autoWhisperContext (whisper_context *p = nullptr) : ptr(p) {}
	~autoWhisperContext ();

	autoWhisperContext (const autoWhisperContext&) = delete;
	autoWhisperContext& operator= (const autoWhisperContext&) = delete;

	autoWhisperContext (autoWhisperContext&& other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperContext& operator= (autoWhisperContext&& other) noexcept;

	[[nodiscard]]
	whisper_context *get () const { return ptr; }
};

struct autoWhisperVadContext {
	whisper_vad_context *ptr;

	autoWhisperVadContext (whisper_vad_context *p = nullptr) : ptr (p) {}
	~autoWhisperVadContext ();

	autoWhisperVadContext (const autoWhisperVadContext&) = delete;
	autoWhisperVadContext& operator= (const autoWhisperVadContext&) = delete;

	autoWhisperVadContext (autoWhisperVadContext&& other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperVadContext& operator= (autoWhisperVadContext&& other) noexcept;

	[[nodiscard]]
	whisper_vad_context *get () const { return ptr; }
};

struct autoWhisperVadSegments {
	whisper_vad_segments *ptr;

	autoWhisperVadSegments (whisper_vad_segments *p = nullptr) : ptr (p) {}
	~autoWhisperVadSegments ();

	autoWhisperVadSegments (const autoWhisperVadSegments&) = delete;
	autoWhisperVadSegments& operator= (const autoWhisperVadSegments&) = delete;

	autoWhisperVadSegments (autoWhisperVadSegments && other) noexcept : ptr (other.ptr) {
		other.ptr = nullptr;
	}
	autoWhisperVadSegments& operator= (autoWhisperVadSegments && other) noexcept;

	[[nodiscard]]
	whisper_vad_segments *get () const { return ptr; }
};

struct autoDiarizeContext {
	diarize_context *ptr;

	autoDiarizeContext (diarize_context *p = nullptr) : ptr(p) {}
	~autoDiarizeContext ();

	autoDiarizeContext (const autoDiarizeContext&) = delete;
	autoDiarizeContext& operator= (const autoDiarizeContext&) = delete;

	autoDiarizeContext (autoDiarizeContext&& other) noexcept : ptr(other.ptr) {
		other.ptr = nullptr;
	}
	autoDiarizeContext& operator= (autoDiarizeContext&& other) noexcept;

	[[nodiscard]]
	diarize_context *get () const { return ptr; }
};


struct SileroVadParams {
	double speechProbabilityThreshold = theVadDefaultThreshold;   // probability threshold to decide that sound is speech
	double minSpeechDuration = theVadDefaultMinSpeechDuration;   // min duration of a speech segment
	double minNonSpeechDuration = theVadDefaultMinNonSpeechDuration;   // min duration of a non-speech segment
	double speechPad = theVadDefaultSpeechPad;   // padding added before and after each speech segment
};

struct DiarizationParams {
	integer maxSimultaneousSpeakers = theDiarizationMaxSimultaneousSpeakers;
	integer numSpeakers = theDiarizationNumSpeakers;
	integer maxSpeakers = theDiarizationMaxSpeakers;
	integer minSpeakers = theDiarizationMinSpeakers;
	double clusterThreshold = theDiarizationClusterThreshold;
	integer segmentationOverlap = theDiarizationSegmentationOverlap;
};

struct SpeechSegment {
	autostring32 text;
	double tmin;
	double tmax;
};

struct WhisperTranscription {
	SpeechSegment fullTranscription;
	autovector <SpeechSegment> words;
	autovector <SpeechSegment> sentences;
};

inline std::set <structSpeechRecognizer *> theLivingSpeechRecognizers;
#include "SpeechRecognizer_def.h"

/*
	Functions to access lists of models and languages from everywhere.
*/
constSTRVEC theCurrentSpeechRecognizerModelNames ();
constSTRVEC theSpeechRecognizerLanguageNames ();

/*
	Class SpeechRecognizer functions.
*/
autoSpeechRecognizer SpeechRecognizer_create (conststring32 modelName, conststring32 languageName);
WhisperTranscription SpeechRecognizer_recognize (constSpeechRecognizer me, constSound sound,
		bool useVad, SileroVadParams const& sileroVadParams);

/*
	Silero-VAD functions.
*/
autovector <SpeechSegment> doSileroVad (constSound sound, SileroVadParams const& sileroVadParams,
		conststring32 nonSpeechLabel, conststring32 speechLabel);

/*
	Diarization functions.
*/
/*
	Run speaker diarization on Sound. Return a vector of timelines, one timeline per speaker.
	Each timeline covers [sound->xmin, sound->xmax] as alternating silent/active intervals,
	labeled `nonSpeechLabel` and `speechLabel` respectively.
*/
autovector <autovector <SpeechSegment>> doDiarization (constSound sound, DiarizationParams const& diarizationParams,
		conststring32 nonSpeechLabel, conststring32 speechLabel);


/* End of file SpeechRecognizer.h */
#endif
