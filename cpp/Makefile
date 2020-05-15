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
	./BendingAngleValidation
	# Test 2: Determine the maximum size of the timestep
	./TimeStepEstimation
	# Test 3: Determine the effect of equilbrating
	./EquilibrateTest
	### Run the adsorption experiments ###
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
	matlab -nodesktop -nosplash -r "run('../analysis/PersistenceLength.m'); exit" > /dev/null 2>&1
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

	### Visualize three example chains ###
	mkdir -p ../figures/Figure_4 > /dev/null 2>&1
	#
	cp data/BendingAngleValidation/theta_0.300_pi/repeat_6/Render.png ../figures/Figure_4/theta_0_300_pi.png
	cp data/BendingAngleValidation/theta_0.600_pi/repeat_6/Render.png ../figures/Figure_4/theta_0_600_pi.png
	cp data/BendingAngleValidation/theta_0.900_pi/repeat_2/Render.png ../figures/Figure_4/theta_0_900_pi.png
	cp data/SphericalColony/N_100/repeat_0/Render.png ../figures/Figure_4/SphericalColony.png
	#
	convert ../figures/Figure_4/theta_0_300_pi.png  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_300_pi.png
	convert ../figures/Figure_4/theta_0_600_pi.png  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_600_pi.png
	convert ../figures/Figure_4/theta_0_900_pi.png  -crop 1800x1400+0+200 +repage ../figures/Figure_4/theta_0_900_pi.png
	convert ../figures/Figure_4/SphericalColony.png -crop 1800x1400+0+200 +repage ../figures/Figure_4/SphericalColony.png
	#
	convert ../figures/Figure_4/theta_0_300_pi.png  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.15 π"   ../figures/Figure_4/theta_0_300_pi.png
	convert ../figures/Figure_4/theta_0_600_pi.png  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.30 π"   ../figures/Figure_4/theta_0_600_pi.png
	convert ../figures/Figure_4/theta_0_900_pi.png  -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Θ / 2 = 0.45 π"   ../figures/Figure_4/theta_0_900_pi.png
	convert ../figures/Figure_4/SphericalColony.png -gravity north -stroke none -fill black -pointsize 200 -font FreeSans -annotate 0 "Spherical colony" ../figures/Figure_4/SphericalColony.png
	#
	convert ../figures/Figure_4/theta_0_300_pi.png ../figures/Figure_4/theta_0_600_pi.png ../figures/Figure_4/theta_0_900_pi.png ../figures/Figure_4/SphericalColony.png +append ../figures/Figure_4/Fig4.png
	#
	rm ../figures/Figure_4/theta_0_300_pi.png ../figures/Figure_4/theta_0_600_pi.png ../figures/Figure_4/theta_0_900_pi.png ../figures/Figure_4/SphericalColony.png
	#
	### Analyze the adsorption rate images ###
	#
	convert ../figures/Figure_5/Fig5a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_5/Fig5a.png
	convert ../figures/Figure_5/Fig5b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_5/Fig5b.png
	convert ../figures/Figure_5/Fig5c.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "C" ../figures/Figure_5/Fig5c.png
	#
	convert ../figures/Figure_5/Fig5a.png ../figures/Figure_5/Fig5b.png ../figures/Figure_5/Fig5c.png +append ../figures/Figure_5/Fig5.png
	#
	rm ../figures/Figure_5/Fig5a.png ../figures/Figure_5/Fig5b.png ../figures/Figure_5/Fig5c.png
	#
	convert ../figures/Figure_S8/FigS8a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_S8/FigS8a.png
	convert ../figures/Figure_S8/FigS8b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_S8/FigS8b.png
	#
	convert -size 950x656 xc:white ../figures/Figure_S8/FigS8a.png -composite ../figures/Figure_S8/FigS8a.png
	#
	convert ../figures/Figure_S8/FigS8a.png ../figures/Figure_S8/FigS8b.png +append ../figures/Figure_S8/FigS8.png
	#
	rm ../figures/Figure_S8/FigS8a.png ../figures/Figure_S8/FigS8b.png
	#
	### Analyze the shielding effect images ###
	#
	convert ../figures/Figure_6/Fig6a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0     "A" ../figures/Figure_6/Fig6a.png
	convert ../figures/Figure_6/Fig6b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate +50+0 "B" ../figures/Figure_6/Fig6b.png
	convert ../figures/Figure_6/Fig6c.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate +50+0 "C" ../figures/Figure_6/Fig6c.png
	#
	convert ../figures/Figure_6/Fig6a.png ../figures/Figure_6/Fig6b.png ../figures/Figure_6/Fig6c.png +append ../figures/Figure_6/Fig6.png
	#
	rm ../figures/Figure_6/Fig6a.png ../figures/Figure_6/Fig6b.png ../figures/Figure_6/Fig6c.png
	#
	### Analyze the secondary infection images ###
	#
	convert ../figures/Figure_7/Fig7a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_7/Fig7a.png
	convert ../figures/Figure_7/Fig7b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_7/Fig7b.png
	convert ../figures/Figure_7/Fig7c.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "C" ../figures/Figure_7/Fig7c.png
	convert ../figures/Figure_7/Fig7d.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "D" ../figures/Figure_7/Fig7d.png
	#
	convert ../figures/Figure_7/Fig7c.png ../figures/Figure_7/Fig7d.png -append ../figures/Figure_7/Fig7cd.png
	#
	convert ../figures/Figure_7/Fig7a.png ../figures/Figure_7/Fig7b.png ../figures/Figure_7/Fig7cd.png +append ../figures/Figure_7/Fig7.png
	#
	rm ../figures/Figure_7/Fig7a.png ../figures/Figure_7/Fig7b.png ../figures/Figure_7/Fig7c.png ../figures/Figure_7/Fig7d.png ../figures/Figure_7/Fig7cd.png
	#
	convert ../figures/Figure_8/Fig8a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_8/Fig8a.png
	convert ../figures/Figure_8/Fig8b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_8/Fig8b.png
	#
	convert ../figures/Figure_8/Fig8a.png ../figures/Figure_8/Fig8b.png +append ../figures/Figure_8/Fig8.png
	#
	rm ../figures/Figure_8/Fig8a.png ../figures/Figure_8/Fig8b.png
	#
	### Chain characterization ###
	#
	convert ../figures/Figure_S3/FigS3a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_S3/FigS3a.png
	convert ../figures/Figure_S3/FigS3b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_S3/FigS3b.png
	#
	convert ../figures/Figure_S3/FigS3a.png ../figures/Figure_S3/FigS3b.png +append ../figures/Figure_S3/FigS3.png
	#
	rm ../figures/Figure_S3/FigS3a.png ../figures/Figure_S3/FigS3b.png
	#
	convert ../figures/Figure_S4/FigS4a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_S4/FigS4a.png
	convert ../figures/Figure_S4/FigS4b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_S4/FigS4b.png
	#
	convert ../figures/Figure_S4/FigS4a.png ../figures/Figure_S4/FigS4b.png +append ../figures/Figure_S4/FigS4.png
	#
	rm ../figures/Figure_S4/FigS4a.png ../figures/Figure_S4/FigS4b.png
	#
	convert ../figures/Figure_S5/FigS5a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_S5/FigS5a.png
	convert ../figures/Figure_S5/FigS5b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_S5/FigS5b.png
	convert ../figures/Figure_S5/FigS5c.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "C" ../figures/Figure_S5/FigS5c.png
	#
	convert ../figures/Figure_S5/FigS5a.png ../figures/Figure_S5/FigS5b.png ../figures/Figure_S5/FigS5c.png +append ../figures/Figure_S5/FigS5.png
	#
	rm ../figures/Figure_S5/FigS5a.png ../figures/Figure_S5/FigS5b.png ../figures/Figure_S5/FigS5c.png
	#
	### Equilibriate testing ###
	#
	convert ../figures/Figure_S9/FigS9a.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" ../figures/Figure_S9/FigS9a.png
	convert ../figures/Figure_S9/FigS9b.png -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" ../figures/Figure_S9/FigS9b.png
	#
	convert ../figures/Figure_S9/FigS9a.png ../figures/Figure_S9/FigS9b.png +append ../figures/Figure_S9/FigS9.png
	#
	rm ../figures/Figure_S9/FigS9a.png ../figures/Figure_S9/FigS9b.png
	#

clean :
	rm -f Simulation.o
	rm -f Chain
	rm -f WellMixed
	rm -f SphericalColony
	rm -f TimeStepEstimation
	rm -f BendingAngleValidation
	rm -f EquilibrateTest
