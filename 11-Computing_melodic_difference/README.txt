Requirements (Ubuntu 18.04.6 LTS):
    python (3.7.10)
	pretty_midi (0.2.10)
	scipy (1.15.3)
	numpy (2.0.1)
    Procodon (Module in this work)

Directories:
	synIXR.SCRaMbLE.note_info: Note infos of SCRaMbLEd strains.
	synIXR.SCRaMbLE.melody: Audios of SCRaMbLEd strains in MIDI format.

Files:
	Canon.synIXR.150_bit_per_note.txt: Note infos of Canon stored in synIXR.
	Canon.synIXR.150_bit_per_note.mid: Audio of Canon stored in synIXR in MIDI format.

Script:
	01-Computing_audio_distance.py: Compute distance scores and count new melodic fragments of audios from SCRaMbLEd strains.
		Usage: python 01-Computing_audio_distance.py
