# Compiler variables
CC=nvcc
CLINKS= -larmadillo -lcudart
CFLAGS= -Xptxas -O3

all :
	make clean
	make tests
	make experiments
	make run

	make clean

	make process

	make analyze

	make finalize

Simulation.o : Simulation.cu Simulation.hpp
	$(CC) $(CFLAGS) -c Simulation.cu $(CLINKS)

experiments : WellMixed Chain SphericalColony

Chain : Chain.cu Simulation.o
	$(CC) $(CFLAGS) -o Chain Chain.cu Simulation.o $(CLINKS)

WellMixed : WellMixed.cu Simulation.o
	$(CC) $(CFLAGS) -o WellMixed WellMixed.cu Simulation.o $(CLINKS)

SphericalColony : SphericalColony.cu Simulation.o
	$(CC) $(CFLAGS) -o SphericalColony SphericalColony.cu Simulation.o $(CLINKS)

tests : TimeStepEstimation BendingAngleValidation EquilibrateTest

TimeStepEstimation : TimeStepEstimation.cu Simulation.o
	$(CC) $(CFLAGS) -o TimeStepEstimation TimeStepEstimation.cu Simulation.o $(CLINKS)

BendingAngleValidation : BendingAngleValidation.cu Simulation.o
	$(CC) $(CFLAGS) -o BendingAngleValidation BendingAngleValidation.cu Simulation.o $(CLINKS)

EquilibrateTest : EquilibrateTest.cu Simulation.o
	$(CC) $(CFLAGS) -o EquilibrateTest EquilibrateTest.cu Simulation.o $(CLINKS)

run : tests experiments

	### Run the tests ###
	# Test 1: Angle between cells, volume scaling and visual examples
	#
	./BendingAngleValidation
	#
	# Test 2: Determine the maximum size of the timestep
	#
	./TimeStepEstimation
	#
	# Test 3: Determine the effect of equilbrating
	#
	./EquilibrateTest
	#
	### Run the adsorption experiments ###
	#
	./WellMixed
	./Chain
	./SphericalColony

process :

	### Visualize the chains ###
	#
	python3 ../analysis/CompareChains.py
	#
	### Characterize the chains ###
	#
	python3 ../analysis/OverlapDetection.py
	#
	python3 ../analysis/RandomWalkScaling.py
	#
	python3 ../analysis/PersistenceLengthDetection.py
	#
	### Analyze the shielding effect ###
	#
	python3 ../analysis/HitDetection.py
	#

analyze :

	### Analyze overlaps ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/Overlaps.m'); exit" > /dev/null 2>&1
	#
	### Analyze convergence ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/TimeStep.m'); exit" > /dev/null 2>&1
	#
	### Analyze bending angle and volume scaling ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/VolumeScaling.m'); exit" > /dev/null 2>&1
	#
	### Analyze the scaling of the chains ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/RandomWalkScaling.m'); exit" > /dev/null 2>&1
	#
	### Analyze the adsorption rate ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/AdsorptionRate.m'); exit" > /dev/null 2>&1
	#
	### Analyze the persistence length of the chains ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/PersistenceLength.m'); exit"
	#
	### Analyze the steady-stateness of the adsorption measurements ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/EquilibrateTest.m'); exit" > /dev/null 2>&1
	#
	### Analyze the shielding effect ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/QuantifyShielding.m'); exit" > /dev/null 2>&1
	#
	### Analyze the lysis events ###
	#
	matlab -nodesktop -nosplash -r "run('../analysis/SecondaryInfections.m'); exit" > /dev/null 2>&1
	#

finalize :

	### Prepare illustration
	convert ../figures/Figure_1/Fig1_Export.png -flatten ../figures/Figure_1/Fig1.tif
	#
	convert ../figures/Figure_1/Fig1.tif -gravity northwest -stroke none -fill black -pointsize 300 -annotate +80+80   "A" ../figures/Figure_1/Fig1.tif
	convert ../figures/Figure_1/Fig1.tif -gravity northwest -stroke none -fill black -pointsize 300 -annotate +80+1700 "B" ../figures/Figure_1/Fig1.tif
	#
	convert ../figures/Figure_1/Fig1.tif -density 900  ../figures/Figure_1/Fig1.tif
	convert ../figures/Figure_1/Fig1.tif -compress lzw ../figures/Figure_1/Fig1.tif
	#
	convert ../figures/Figure_2/Fig2_Export.png -flatten ../figures/Figure_2/Fig2.tif
	#
	convert ../figures/Figure_2/Fig2.tif -density 900  ../figures/Figure_2/Fig2.tif
	convert ../figures/Figure_2/Fig2.tif -compress lzw ../figures/Figure_2/Fig2.tif
	#
	### Visualize three example chains ###
	mkdir -p ../figures/Figure_4 > /dev/null 2>&1
	#
	cp data/BendingAngleValidation/theta_0.300_pi/repeat_6/Render.png ../figures/Figure_4/theta_0_300_pi.tif
	cp data/BendingAngleValidation/theta_0.600_pi/repeat_6/Render.png ../figures/Figure_4/theta_0_600_pi.tif
	cp data/BendingAngleValidation/theta_0.900_pi/repeat_2/Render.png ../figures/Figure_4/theta_0_900_pi.tif
	cp data/SphericalColony/N_100/repeat_0/Render.png ../figures/Figure_4/SphericalColony.tif
	#
	convert ../figures/Figure_4/theta_0_300_pi.tif  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_300_pi.tif
	convert ../figures/Figure_4/theta_0_600_pi.tif  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_600_pi.tif
	convert ../figures/Figure_4/theta_0_900_pi.tif  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_900_pi.tif
	convert ../figures/Figure_4/SphericalColony.tif -crop 1800x1400+0+200 +repage ../figures/Figure_4/SphericalColony.tif
	#
	convert ../figures/Figure_4/theta_0_300_pi.tif  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.15 π"   ../figures/Figure_4/theta_0_300_pi.tif
	convert ../figures/Figure_4/theta_0_600_pi.tif  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.30 π"   ../figures/Figure_4/theta_0_600_pi.tif
	convert ../figures/Figure_4/theta_0_900_pi.tif  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.45 π"   ../figures/Figure_4/theta_0_900_pi.tif
	convert ../figures/Figure_4/SphericalColony.tif -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Spherical colony" ../figures/Figure_4/SphericalColony.tif
	#
	convert ../figures/Figure_4/theta_0_300_pi.tif ../figures/Figure_4/theta_0_600_pi.tif ../figures/Figure_4/theta_0_900_pi.tif ../figures/Figure_4/SphericalColony.tif +append ../figures/Figure_4/Fig4.tif
	convert ../figures/Figure_4/Fig4.tif -density 900  ../figures/Figure_4/Fig4.tif
	convert ../figures/Figure_4/Fig4.tif -compress lzw ../figures/Figure_4/Fig4.tif
	#
	rm ../figures/Figure_4/theta_0_300_pi.tif ../figures/Figure_4/theta_0_600_pi.tif ../figures/Figure_4/theta_0_900_pi.tif ../figures/Figure_4/SphericalColony.tif
	#
	### Analyze the adsorption rate images ###
	#
	convert ../figures/Figure_5/Fig5a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_5/Fig5a.tif
	convert ../figures/Figure_5/Fig5b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_5/Fig5b.tif
	convert ../figures/Figure_5/Fig5c.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "C" ../figures/Figure_5/Fig5c.tif
	#
	convert ../figures/Figure_5/Fig5a.tif ../figures/Figure_5/Fig5b.tif ../figures/Figure_5/Fig5c.tif +append ../figures/Figure_5/Fig5.tif
	convert ../figures/Figure_5/Fig5.tif -density 900  ../figures/Figure_5/Fig5.tif
	convert ../figures/Figure_5/Fig5.tif -compress lzw ../figures/Figure_5/Fig5.tif
	#
	rm ../figures/Figure_5/Fig5a.tif ../figures/Figure_5/Fig5b.tif ../figures/Figure_5/Fig5c.tif
	#
	convert ../figures/Figure_S8/FigS8a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_S8/FigS8a.tif
	convert ../figures/Figure_S8/FigS8b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_S8/FigS8b.tif
	#
	convert -size 5250x3938 xc:white ../figures/Figure_S8/FigS8a.tif -composite ../figures/Figure_S8/FigS8a.tif
	#
	convert ../figures/Figure_S8/FigS8a.tif ../figures/Figure_S8/FigS8b.tif +append ../figures/Figure_S8/FigS8.tif
	convert ../figures/Figure_S8/FigS8.tif -density 900  ../figures/Figure_S8/FigS8.tif
	convert ../figures/Figure_S8/FigS8.tif -compress lzw ../figures/Figure_S8/FigS8.tif
	#
	rm ../figures/Figure_S8/FigS8a.tif ../figures/Figure_S8/FigS8b.tif
	#
	### Analyze the shielding effect images ###
	#
	convert ../figures/Figure_6/Fig6a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_6/Fig6a.tif
	convert ../figures/Figure_6/Fig6b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_6/Fig6b.tif
	convert ../figures/Figure_6/Fig6c.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "C" ../figures/Figure_6/Fig6c.tif
	#
	convert ../figures/Figure_6/Fig6a.tif ../figures/Figure_6/Fig6b.tif ../figures/Figure_6/Fig6c.tif +append ../figures/Figure_6/Fig6.tif
	convert ../figures/Figure_6/Fig6.tif -density 900  ../figures/Figure_6/Fig6.tif
	convert ../figures/Figure_6/Fig6.tif -compress lzw ../figures/Figure_6/Fig6.tif
	#
	rm ../figures/Figure_6/Fig6a.tif ../figures/Figure_6/Fig6b.tif ../figures/Figure_6/Fig6c.tif
	#
	### Analyze the secondary infection images ###
	#
	convert ../figures/Figure_7/Fig7a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_7/Fig7a.tif
	convert ../figures/Figure_7/Fig7b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_7/Fig7b.tif
	convert ../figures/Figure_7/Fig7c.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "C" ../figures/Figure_7/Fig7c.tif
	convert ../figures/Figure_7/Fig7d.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "D" ../figures/Figure_7/Fig7d.tif
	#
	convert ../figures/Figure_7/Fig7c.tif ../figures/Figure_7/Fig7d.tif -append ../figures/Figure_7/Fig7cd.tif
	#
	convert ../figures/Figure_7/Fig7a.tif ../figures/Figure_7/Fig7b.tif ../figures/Figure_7/Fig7cd.tif +append ../figures/Figure_7/Fig7.tif
	convert ../figures/Figure_7/Fig7.tif -density 900  ../figures/Figure_7/Fig7.tif
	convert ../figures/Figure_7/Fig7.tif -compress lzw ../figures/Figure_7/Fig7.tif
	#
	rm ../figures/Figure_7/Fig7a.tif ../figures/Figure_7/Fig7b.tif ../figures/Figure_7/Fig7c.tif ../figures/Figure_7/Fig7d.tif ../figures/Figure_7/Fig7cd.tif
	#
	convert ../figures/Figure_8/Fig8a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_8/Fig8a.tif
	convert ../figures/Figure_8/Fig8b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_8/Fig8b.tif
	#
	convert ../figures/Figure_8/Fig8a.tif ../figures/Figure_8/Fig8b.tif +append ../figures/Figure_8/Fig8.tif
	convert ../figures/Figure_8/Fig8.tif -density 900  ../figures/Figure_8/Fig8.tif
	convert ../figures/Figure_8/Fig8.tif -compress lzw ../figures/Figure_8/Fig8.tif
	#
	rm ../figures/Figure_8/Fig8a.tif ../figures/Figure_8/Fig8b.tif
	#
	### Chain characterization ###
	#
	convert ../figures/Figure_S3/FigS3a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_S3/FigS3a.tif
	convert ../figures/Figure_S3/FigS3b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_S3/FigS3b.tif
	#
	convert ../figures/Figure_S3/FigS3a.tif ../figures/Figure_S3/FigS3b.tif +append ../figures/Figure_S3/FigS3.tif
	convert ../figures/Figure_S3/FigS3.tif -density 900  ../figures/Figure_S3/FigS3.tif
	convert ../figures/Figure_S3/FigS3.tif -compress lzw ../figures/Figure_S3/FigS3.tif
	#
	rm ../figures/Figure_S3/FigS3a.tif ../figures/Figure_S3/FigS3b.tif
	#
	convert ../figures/Figure_S4/FigS4a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_S4/FigS4a.tif
	convert ../figures/Figure_S4/FigS4b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_S4/FigS4b.tif
	#
	convert ../figures/Figure_S4/FigS4a.tif ../figures/Figure_S4/FigS4b.tif +append ../figures/Figure_S4/FigS4.tif
	convert ../figures/Figure_S4/FigS4.tif -density 900  ../figures/Figure_S4/FigS4.tif
	convert ../figures/Figure_S4/FigS4.tif -compress lzw ../figures/Figure_S4/FigS4.tif
	#
	rm ../figures/Figure_S4/FigS4a.tif ../figures/Figure_S4/FigS4b.tif
	#
	convert ../figures/Figure_S5/FigS5a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_S5/FigS5a.tif
	convert ../figures/Figure_S5/FigS5b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_S5/FigS5b.tif
	convert ../figures/Figure_S5/FigS5c.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "C" ../figures/Figure_S5/FigS5c.tif
	#
	convert ../figures/Figure_S5/FigS5a.tif ../figures/Figure_S5/FigS5b.tif ../figures/Figure_S5/FigS5c.tif +append ../figures/Figure_S5/FigS5.tif
	convert ../figures/Figure_S5/FigS5.tif -density 900  ../figures/Figure_S5/FigS5.tif
	convert ../figures/Figure_S5/FigS5.tif -compress lzw ../figures/Figure_S5/FigS5.tif
	#
	rm ../figures/Figure_S5/FigS5a.tif ../figures/Figure_S5/FigS5b.tif ../figures/Figure_S5/FigS5c.tif
	#
	### Equilibriate testing ###
	#
	convert ../figures/Figure_S9/FigS9a.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "A" ../figures/Figure_S9/FigS9a.tif
	convert ../figures/Figure_S9/FigS9b.tif -gravity northwest -stroke none -fill black -pointsize 384 -annotate +0+0 "B" ../figures/Figure_S9/FigS9b.tif
	#
	convert ../figures/Figure_S9/FigS9a.tif ../figures/Figure_S9/FigS9b.tif +append ../figures/Figure_S9/FigS9.tif
	convert ../figures/Figure_S9/FigS9.tif -density 900  ../figures/Figure_S9/FigS9.tif
	convert ../figures/Figure_S9/FigS9.tif -compress lzw ../figures/Figure_S9/FigS9.tif
	#
	rm ../figures/Figure_S9/FigS9a.tif ../figures/Figure_S9/FigS9b.tif
	#
	### Make a PNG copy of .tif images ###
	find .. -name '*.tif' -print0 | xargs -0 -r mogrify -resize 25% -format png

clean :
	rm -f Simulation.o
	rm -f Chain
	rm -f WellMixed
	rm -f SphericalColony
	rm -f TimeStepEstimation
	rm -f BendingAngleValidation
	rm -f EquilibrateTest
