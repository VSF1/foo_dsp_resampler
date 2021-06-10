# foo_dsp_resampler
Foobar dsp resampler

Good quality and fast resampler based on code from SoX `rate` effect.

If you want to know what settings are best:
1. Read SoX FAQ, "What are the best 'rate' settings to resample a file and retain the highest quality?"
2. This post: http://www.hydrogenaudio.org/forums/index....st&p=626176 (an excerpt from SoX help)
3. Feel free to experiment and decide what's best for you.

Changelog:

0.8.7:
- fixed a bug introduced in 0.8.4 that may result in crashes.
- recompiled with MSVS 2019 (16.2.0) and SDK-2019-06-30.

0.8.6:
- fixed a bug introduced in 0.8.5 when "Target samplerate" was set to "Upsample" or "Downsample".

0.8.5: maintenance release
- minor code cleanup;
- minor update of project files;
- update FFmpeg FFT code.
No bugfixes / speed increase / new features; 'mod' versions were not updated.

0.8.4:
- minor fixes;
- changed default quality from 'Normal' to 'Best';
- added signal extrapolation to reduce a possibility of clicks during track change;
- recompiled with latest SDK (2018-10-11). Works with foobar2000 1.4.x; users of 1.3.x need to download and install MSVC redist (here).
