# multiplepassage

After a large (M7ish) quake:
* obtain waveform data from a broadband station
* plot a spectrogram (OK, actually a continuous wavelet transform)
* and optionally overlay theoretical Rayleigh wave arrival times
to demonstrate multiple passages of surface waves
(and make a pretty plot)

![M8.1 Kermadec quake at IU.ANMO](https://github.com/ewolin/multiplepassage/blob/main/multiplepassage.png?raw=true)



Example script, waveforms, and plot provided for March 2021 M8.1 earthquake near the Kermadec Islands, New Zealand recorded on GSN station IU.ANMO.

requires ObsPy:

https://docs.obspy.org/

https://github.com/obspy/obspy/wiki#installation


Predicted Rayleigh wave arrival times in t-vs-u.csv were digitized from Oliver (1962) Fig 2: https://pubs.geoscienceworld.org/ssa/bssa/article/52/1/81/101314
