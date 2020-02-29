#! /bin/sh

for f in output/animWave/*.svg; do
	echo "$f"
	convert "$f" "${f%.svg}".png
done

ffmpeg -f image2 -framerate 30 -i output/animWave/frame-%04d.png output/animWave/vid.mp4
