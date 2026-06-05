/* manual_SpeechRecognizer.cpp
 *
 * Copyright (C) 2026 Anastasia Shchupak
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

#include "ManPagesM.h"

void manual_SpeechRecognizer_init (ManPages me);
void manual_SpeechRecognizer_init (ManPages me) {

MAN_PAGES_BEGIN R"~~~(
################################################################################
"Speech recognition"
© Anastasia Shchupak 2026-06-01

This tutorial describes how you can use the speech recognition tools in Praat.
It assumes that you are familiar with the @Intro, especially with @@Intro 7. Annotation@.

There are two speech recognition tools available in Praat: automatic transcription (turning speech into text) and
automatic speaker diarization (finding out who speaks when). Transcription requires at least one model file to be
installed in the Praat preferences folder; diarization works without any external models. Diarization can be
performed together with transcription, so that transcribed text is divided between speakers. Sections 1 and 2
describe transcription and diarization respectively. Section 3 then explains what can be done in an attempt to improve
the performance of these tools in Praat.

1. Automatic transcription
==========================
Transcription in Praat turns the speech in a @Sound object into written text. Usually you might want this text to be
placed into @TextGrid's intervals, aligning it with the Sound (this is described in subsection 1b). But you can also
obtain it as plain text (subsection 1c is about how to do that). Either way, you first need a Whisper model in ggml
format installed, and how to do that is explained in subsection 1a.

1a. Installing the models
=========================


1b. Transcribing into a TextGrid
================================


1c. Transcribing into plain text
================================


2. Automatic speaker diarization
================================


3. Performance
==============


################################################################################
"speech activity detection with Silero VAD"
© Anastasia Shchupak 2026-06-01

Praat uses the @@whisper.cpp@ implementation of the @@Silero VAD@ speech activity detector.
The pre-trained Silero VAD model weights are converted to ggml format and compiled into Praat,
so no external model files are required. The sound is automatically resampled to 16 kHz
(the sampling frequency expected by the Silero VAD model) before being processed by the model.

Purpose
=======
to detect which portions of a sound contain speech. The output is a list of speech segments,
each defined by a start and an end time.

Algorithm
=========
The Silero VAD model processes the audio in fixed frames of 512 samples (32 ms, since the sound
is resampled to 16 kHz). For each frame, it outputs a probability that the frame contains speech.
Based on this output, the list of speech segments is constructed as described below. This description reflects
the whisper.cpp implementation, which differs slightly from Silero's original implementation.

A speech segment begins at the first frame whose speech probability exceeds the ##Speech probability
threshold#. The segment continues as long as the speech probability stays above a lower threshold
(called the "negative threshold" in the Silero source code, equal to the speech probability threshold minus 0.15).
The segment ends only when the probability stays below the negative threshold for at least ##Min. gap between speech segments#.
When a segment ends, it is kept only if it is longer than ##Min. speech segment#.

After all segments have been formed, segments separated by a gap shorter than 0.2 s are merged (this is hardcoded
in whisper.cpp and not configurable). After that, padding is applied: each segment is extended on both sides
by ##Padding around speech segments#. If padding would cause two segments to overlap, they instead meet
at the midpoint of the gap between them.

The result is a list of speech segments.

Settings
========
##Speech probability threshold (0-1)
:   determines the sensitivity of the speech detector. Higher values make the detector less sensitive,
	meaning that a frame requires a higher speech probability to be considered part of a speech segment.
	This reduces false positives (non-speech incorrectly classified as speech),	but may cause some speech
	to be missed. Lower values make the detector more sensitive. The default of 0.5 works well for most use cases.

##Min. gap between speech segments (s)
:   the minimum duration of a gap between two speech segments. You might want to increase this value if
	short silences within speech (e.g. plosive closures) are splitting speech into multiple segments. Note that
	gaps shorter than 0.2 s are removed from the output (with their adjacent speech segments merged),
    regardless of this setting.

##Min. speech segment (s)
:   the minimum duration of a speech segment. Shorter segments are discarded.

##Padding around speech segments (s)
:   extends each detected speech segment by this amount on both sides. You might want to increase this value if
	speech onsets and offsets are being clipped.

Availability in Praat
=====================
Silero VAD speech activity detection is available in two ways in Praat:
- as part of transcription, see @@transcription with whisper.cpp@;
- standalone, see @@Sound: To TextGrid (speech activity, Silero)...@.

################################################################################
"Sound: To TextGrid (speech activity, Silero)..."
© Anastasia Shchupak 2026-06-01

A command that creates a @TextGrid with one interval tier from every selected @Sound object.
The interval tier contains speech and non-speech intervals with boundaries determined by the
Silero VAD model (for the algorithm and settings, see @@speech activity detection with Silero VAD@).
The labels of the intervals are specified by the settings ##Non-speech interval label# and ##Speech interval label#.

Settings
========
##Speech probability threshold (0-1)#

##Min. gap between speech segments (s)#

##Min. speech segment (s)#

##Padding around speech segments (s)
:   see @@speech activity detection with Silero VAD@.

##Non-speech interval label
:	the label assigned to intervals classified as non-speech in the resulting TextGrid.

##Speech interval label
:	the label assigned to intervals classified as speech in the resulting TextGrid.

################################################################################
"transcription with whisper.cpp"
© Anastasia Shchupak 2026-06-01

Praat performs speech transcription using @@whisper.cpp@. This page documents the transcription settings.
If you are new to speech recognition in Praat, see the @@Speech recognition@ tutorial first.

Behaviour
=========
When transcription is run on a @TextGrid interval, it modifies the TextGrid: intervals are split,
tiers may be added or renamed. The exact resulting TextGrid structure depends on the combination
of the ##Include words# and ##Include diarization# settings. See the @@Speech recognition@ tutorial
for the details of the transcription output under different combinations of these settings.

Settings
========
##Whisper model
:	determines which Whisper model is used.
	The list is populated with the `.bin` files found in the `whispercpp` subfolder
	of the `models` folder in the Praat preferences folder. See the @@Speech recognition@
	tutorial for details on how to install models.

##Language
:	determines the language to be used for transcription.
	Choose `Autodetect language` to let the model detect the language automatically.
	If you know the language you want to use for transcription, selecting it explicitly may improve transcription accuracy.
	Note that English-only models (those with `.en` in the name) can only be used with
	`Autodetect language` or `English`.

##Include words
:	if on, each transcribed word is given a start and an end time, computed using whisper.cpp's internal
	dynamic time warping (DTW) algorithm.

##Include diarization
:   if on, speaker diarization is run alongside transcription (see @@speaker diarization with adapted pyannote.audio@).
	The results of both are later combined to attribute portions of transcribed speech to different speakers.

##Allow silences
:	if on, @@speech activity detection with Silero VAD@ runs before transcription to identify speech regions.
	Only those regions are then passed to the Whisper model. This generally improves both speed and accuracy of transcription.
	Speed is improved by reducing the length of the audio sent to the model, and accuracy by preventing the model from
	hallucinating text for silent regions.

##Speech probability threshold (0-1)#

##Min. gap between speech segments (s)#

##Min. speech segment (s)#

##Padding around speech segments (s)
:   see @@speech activity detection with Silero VAD@.

##Fixed number of speakers...#

##... or range of numbers of speakers#

##Allow speakers overlap#

##Clustering threshold (0-2)#

##Segmentation step (0-1)
:   see @@speaker diarization with adapted pyannote.audio@.

Availability in Praat
=====================
Transcription with whisper.cpp is available in two ways in Praat:
- if you select a @Sound together with its @TextGrid and choose
  @@TextGrid & Sound: Transcribe interval...|Transcribe interval...@;
- via ##Transcribe interval# from the #Interval menu in the @TextGridEditor.
  The settings for this command are set via ##Transcription settings...# in the same menu
  and are remembered across Praat sessions.

################################################################################
"TextGrid & Sound: Transcribe interval..."
© Anastasia Shchupak 2026-06-01

This command transcribes the portion of the @Sound that corresponds to a specified interval of
the @TextGrid and writes the result back into the TextGrid.

Settings
========
##Tier number
:	the number of the tier containing the interval to be transcribed.

##Interval number
:	the number of the interval to be transcribed.

The remaining settings in this dialog control the transcription itself; see @@transcription with whisper.cpp@
for their meaning.

################################################################################
"speaker diarization with adapted pyannote.audio"
© Anastasia Shchupak 2026-06-01

Praat performs speaker diarization using a C++/ggml adaptation of @@pyannote.audio@'s
`pyannote/speaker-diarization-3.1` pipeline. The two neural models used by the pipeline,
`pyannote/segmentation-3.0` (for segmentation) and `wespeaker-voxceleb-resnet34-LM` (for
speaker embedding, see @@WeSpeaker@), are converted to ggml format and compiled into Praat,
so no external model files are required. The sound is automatically resampled to 16 kHz
(the sampling frequency expected by both models) before being processed.

Purpose
=======
to detect which portions of a sound contain speech and to attribute each portion to one or more speakers.
The output is a list of segments, each segment defined by a start time, an end time and a speaker identifier.
Speaker identifiers are natural numbers (1, 2, 3, ...); they are arbitrary labels that roughly follow the order
in which the speakers first appear.

Settings
========
##Fixed number of speakers...
:   if set to a positive integer, the algorithm aims to produce exactly this many speakers. Takes precedence
	over the ##... or range of numbers of speakers# setting below. This is a target, not a guarantee: if the Sound
	does not have that many distinct voices, fewer speakers may be produced (and a warning is shown).

##... or range of numbers of speakers
:   a minimum and a maximum for the number of speakers the algorithm may produce. You can set only the minimum,
	only the maximum, or both; a value of 0 for any of the bounds leaves that bound unconstrained.
	Note that if ##Fixed number of speakers...# is set, this setting will be ignored.

##Allow speakers overlap
:   if on, two speakers may be active at the same moment (not more than two though, as this is the limit of the
	segmentation model). If off, every moment is attributed to a single speaker, the one most active at that moment,
	according to the segmentation model.

##Clustering threshold (0-2)
:	if two groups at stage 3 of the algorithm are closer than this threshold, then they are merged. Lower values
	produce more speakers, higher values produce fewer speakers. When the number of speakers is constrained by ##Fixed
	number of speakers...# or by ##... or range of numbers of speakers#, this threshold no longer determines how
	many speakers are found (however, it still influences the process of group forming).

##Segmentation step (0-1)
:   the spacing between consecutive overlapping 10-second chunks, as a fraction of the chunk length (see stage 1 of
	the algorithm). A value of 1.0 means no overlap; lower values mean more overlap. Less overlap is faster (but
	generally less accurate), so this is the setting to lower when you want to speed diarization up.

Algorithm
=========
The algorithm has four stages. It is a port of @@pyannote.audio@'s `pyannote/speaker-diarization-3.1`
pipeline (see @@Bredin (2023)@) with some adaptations.

##1. Segmentation#. The sound is divided into overlapping 10-second chunks. The distance between the starts of two
consecutive chunks is a fraction of the chunk length defined by ##Segmentation step (0-1)#. For example, a
segmentation step of 0.1 makes this distance 1 second, so that consecutive chunks have a 90\% overlap.

Each chunk is then sent to the segmentation model, which divides it into 589 frames (each frame spanning about 0.017
s) and assigns a label to each frame in the following way. The model assumes that there are at most three speakers in a
chunk, of whom at most two can be active in the same frame. For every frame, it gives a probability to each of seven
possible combinations of active speakers: {∅}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}. The most probable combination
is taken as that frame's label. The speaker numbers 1, 2 and 3 are local to each chunk: speaker 1 in one chunk is not
necessarily the same person as speaker 1 in another. Linking these local speakers across chunks is the task of the
remaining stages.

##2. Speaker embeddings#. For each chunk, every speaker who actually speaks in it is turned into an %embedding (a
256-dimensional vector) by the embedding model. The model takes as input an audio piece assembled from the parts of
the chunk where the speaker is active, and outputs that vector. The model works in such a way that two embeddings
from the same speaker (from different chunks) tend to be closer to each other than two embeddings from two different
speakers. This property of the model makes the next stage possible.

##3. Clustering#. The previous stage produces one embedding for each active speaker in each chunk. Of these, only the
ones backed by enough non-overlapping speech are treated as "reliable"; the rest are set aside for now.

The reliable embeddings from all chunks are L2-normalized and grouped by similarity using agglomerative hierarchical
clustering: the two most similar groups are merged repeatedly, stopping when the closest remaining groups are further
apart than ##Clustering threshold (0-2)#. Similarity between two groups is measured by the distance between their
centres (more similar = less distance). Because the embeddings are L2-normalized, the distance between any two of
them lies between 0 (identical) and 2 (opposite).

Each final group is one speaker; the number of groups can be influenced by ##Fixed number of speakers...# and
##... or range of numbers of speakers# settings. Finally, the unreliable embeddings (those not used to form the
groups) are attached to their nearest groups. In this way, each chunk-local speaker is assigned to a global speaker.

##4. Reconstruction#. Each chunk consists of frames, and because the chunks overlap, each frame on the global
timeline is covered by several chunks. Such a frame receives, for every speaker, the average of the activations the
covering chunks gave that speaker (using the global speaker numbers from stage 3). It is then assigned to the speaker
with the highest average activation (or, where overlap is allowed and two people are talking, to the two highest).
Finally, for each speaker, every run of consecutive frames in which that speaker is active becomes one %segment. The
result is a list of segments, each attributed to one speaker.

Availability in Praat
=====================
Speaker diarization is available in two ways in Praat:
- as part of transcription, see @@transcription with whisper.cpp@;
- standalone, if you select a @Sound together with its @TextGrid and choose
  @@TextGrid & Sound: Diarize interval...|Diarize interval...@;
- standalone, via ##Diarize interval# from the #Interval menu in the @TextGridEditor.
  The settings for this command are set via ##Diarization settings...# in the same menu
  and are remembered across Praat sessions.

################################################################################
"TextGrid & Sound: Diarize interval..."
© Anastasia Shchupak 2026-06-01

This command performs speaker diarization on the portion of the @Sound that corresponds to a specified interval of
the @TextGrid and writes the result back into the TextGrid.

Settings
========
##Tier number
:	the number of the tier containing the interval to be diarized.

##Interval number
:	the number of the interval to be diarized.

The remaining settings in this dialog control the speaker diarization itself;
see @@speaker diarization with adapted pyannote.audio@ for their meaning.

################################################################################
"SpeechRecognizer"
© Anastasia Shchupak 2026-03-15

The SpeechRecognizer is one of the @@types of objects@ in Praat.
It performs automatic transcription (speech-to-text) on a @Sound object, using the @@whisper.cpp@ engine,
and therefore our SpeechRecognizer is an interface to whisper.cpp.

If you are new to speech recognition in Praat, see the @@Speech recognition@ tutorial first.

Commands
========

Creation:
,	@@Create SpeechRecognizer...@

Recognition:
,	@@SpeechRecognizer & Sound: Transcribe|Transcribe@

################################################################################
"Create SpeechRecognizer..."
© Anastasia Shchupak 2026-03-15

Creates the @@whisper.cpp@ speech recognizer.
If you are new to speech recognition in Praat, see the @@Speech recognition@ tutorial first.

Settings
========

##Whisper model
:	determines which Whisper model to use for recognition.
The list is populated with the ##.bin# files found in the ##whispercpp# subfolder
of the ##models# folder in the Praat preferences folder.
Models that contain ##.en# in their name are English-only; all other models are multilingual.

##Language
:	determines the language to be recognized.
Choose ##Autodetect language# to let the model detect the language automatically.
If you know the language of the audio, selecting it explicitly may improve recognition accuracy.
Note that English-only models (those with ##.en# in the name) can only be used with
##Autodetect language# or ##English#.

################################################################################
"SpeechRecognizer & Sound: Transcribe"
© Anastasia Shchupak 2026-03-15

Performs speech recognition on the selected @Sound using the selected @SpeechRecognizer,
and writes the recognized text to the Info window.

The sound is automatically resampled to 16 kHz (the sampling frequency expected by @@whisper.cpp@)
if its original sampling frequency differs.

The recognition uses the @@speech activity detection with Silero VAD@ built into @@whisper.cpp@
to skip non-speech parts of the audio, which improves both speed and accuracy.

The result is a flat text string containing the full transcription.

################################################################################
"Silero VAD"
© Anastasia Shchupak 2026-06-01

Silero VAD is a pre-trained speech activity detector created by Silero Team; see @@Silero Team (2024)@
for the model and documentation. VAD stands for "voice activity detection", but we use
"speech activity detection" instead, because Silero VAD detectsunvoiced parts of speech as well as voiced ones.
The Silero VAD model is compiled into Praat (see @Acknowledgments).

For how speech activity detection with Silero VAD is used in Praat, see the @@Speech recognition@ tutorial.
For the algorithm and settings, see @@speech activity detection with Silero VAD@.

################################################################################
"whisper.cpp"
© Anastasia Shchupak 2026-06-01

Whisper is an automatic speech recognition (ASR) system by OpenAI for transcribing speech audio into text;
see @@Radford et al. (2022)@ for the model architecture, training and evaluation. OpenAI's Whisper implementation
is in Python (PyTorch).

Praat uses whisper.cpp, a lightweight C/C++ port of Whisper built on top of the ggml tensor library for machine learning,
developed by Georgi Gerganov and many other contributors (see @Acknowledgments).
The original OpenAI Whisper models are converted to ggml format for use with whisper.cpp.

For how transcription is used in Praat, see the @@Speech recognition@ tutorial.
For the transcription settings, see @@transcription with whisper.cpp@.

################################################################################
"pyannote.audio"
© Anastasia Shchupak 2026-06-01

pyannote.audio is an automatic speaker diarization toolkit developed by Hervé Bredin and collaborators
(see @@Plaquet & Bredin (2023)@ and @@Bredin (2023)@).

Praat contains a C++/ggml adaptation of its `pyannote/speaker-diarization-3.1` pipeline.
The pipeline uses two neural models: `pyannote/segmentation-3.0` by pyannote.audio for segmentation
and `wespeaker-voxceleb-resnet34-LM` by @@WeSpeaker@ for speaker embedding.
The `pyannote/segmentation-3.0` weights have been converted to ggml format
and embedded into Praat (see @Acknowledgments).

For how speaker diarization is used in Praat, see the @@Speech recognition@ tutorial.
For the diarization settings, see @@speaker diarization with adapted pyannote.audio@.

################################################################################
"WeSpeaker"
© Anastasia Shchupak 2026-06-01

WeSpeaker is a speaker embedding learning toolkit developed by WeNet Community
(see @@Wang et al. (2023)@ and @@Wang et al. (2024)@). Praat uses one of its pretrained models,
the `wespeaker-voxceleb-resnet34-LM` embedding model, as part of @@pyannote.audio@'s
`pyannote/speaker-diarization-3.1` pipeline.

`wespeaker-voxceleb-resnet34-LM` was trained on the VoxCeleb2 dataset
(see @@Chung, Nagrani & Zisserman (2018)@ and @@Nagrani, Chung & Zisserman (2017)@).
The model weights have been converted to ggml format and embedded into Praat (see @Acknowledgments).
Praat contains a C++/ggml port of WeSpeaker's ResNet34 architecture with TSTP pooling,
used for inference on this model.

################################################################################
"Silero Team (2024)"
© Anastasia Shchupak 2026-06-01

Silero Team (2024). Silero VAD: pre-trained enterprise-grade voice activity detector (VAD),
number detector and language classifier [Computer software]. Version 6.2.0.

Available on `https://github.com/snakers4/silero-vad`.

################################################################################
"Radford et al. (2022)"
© Anastasia Shchupak 2026-06-01

Alec Radford, Jong Wook Kim, Tao Xu, Greg Brockman, Christine McLeavey & Ilya Sutskever (2022):
"Robust speech recognition via large-scale weak supervision."

Available on `https://cdn.openai.com/papers/whisper.pdf`.

################################################################################
"Plaquet & Bredin (2023)"
© Anastasia Shchupak 2026-06-01

Alexis Plaquet & Hervé Bredin (2023): "Powerset multi-class cross entropy loss for neural
speaker diarization." %%Proc. Interspeech 2023%.

Available on `https://www.isca-archive.org/interspeech_2023/plaquet23_interspeech.html`.

################################################################################
"Bredin (2023)"
© Anastasia Shchupak 2026-06-01

Hervé Bredin (2023): "pyannote.audio 2.1 speaker diarization pipeline: principle, benchmark,
and recipe." %%Proc. Interspeech 2023%.

Available on `https://www.isca-archive.org/interspeech_2023/bredin23_interspeech.html`.

################################################################################
"Wang et al. (2023)"
© Anastasia Shchupak 2026-06-01

Hongji Wang, Chengdong Liang, Shuai Wang, Zhengyang Chen, Binbin Zhang, Xu Xiang, Yanlei Deng
& Yanmin Qian (2023): "Wespeaker: a research and production oriented speaker embedding learning
toolkit." %%ICASSP 2023 - 2023 IEEE International Conference on Acoustics, Speech and Signal
Processing (ICASSP)%, 1–5.

Available on `https://doi.org/10.48550/arXiv.2210.17016`.

################################################################################
"Wang et al. (2024)"
© Anastasia Shchupak 2026-06-01

Shuai Wang, Zhengyang Chen, Bing Han, Hongji Wang, Chengdong Liang, Binbin Zhang, Xu Xiang,
Wen Ding, Johan Rohdin, Anna Silnova, Yanmin Qian & Haizhou Li (2024): "Advancing speaker
embedding learning: Wespeaker toolkit for research and production." %%Speech Communication%
##162#: 103104.

Available on `https://doi.org/10.1016/j.specom.2024.103104`.

################################################################################
"Chung, Nagrani & Zisserman (2018)"
© Anastasia Shchupak 2026-06-01

Joon Son Chung, Arsha Nagrani & Andrew Zisserman (2018): "VoxCeleb2: deep speaker recognition."
%%Proc. Interspeech 2018%.

Available on `https://www.isca-archive.org/interspeech_2018/chung18b_interspeech.html`.

################################################################################
"Nagrani, Chung & Zisserman (2017)"
© Anastasia Shchupak 2026-06-01

Arsha Nagrani, Joon Son Chung & Andrew Zisserman (2017): "VoxCeleb: a large-scale speaker
identification dataset." %%Proc. Interspeech 2017%.

Available on `https://www.isca-archive.org/interspeech_2017/nagrani17_interspeech.html`.

################################################################################
)~~~"
MAN_PAGES_END

}

/* End of file manual_SpeechRecognizer.cpp */
