#!/bin/bash/

srcDir=/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/CosmicShearRB/Cosmo/montepython_light
dataDir=/disks/shear15/ssli/KV450/Cosmo


cd $srcDir/montepython

rm -r $dataDir/kv450_cf/

mkdir $dataDir/kv450_cf/

cp $srcDir/montepython/likelihoods/kv450_cf_likelihood_public/kv450_cf_likelihood_public.data  \
   $dataDir/kv450_cf/


python3 MontePython.py run \
	`# supply relative path from working directory (wd) to param-file` \
	-p $srcDir/input/kv450_cf_best.param \
	`# supply relative path from wd to output folder` \
	-o $dataDir/kv450_cf/ \
	`# supply relative path from wd to correctly set config-file (otherwise default.conf from MontePython will be used)` \
	--conf $srcDir/kv450_cf.conf \
	`# choose the MultiNest sampler (nested sampling)` \
	-m NS \
	`# set an arbitrary but large number of steps (run should converge well before!)` \
	--NS_max_iter 10000000 \
	`# flag for using importance nested sampling (we did not use it, but it might be faster when using it!)` \
	--NS_importance_nested_sampling False \
	`# for parameter estimation use 0.8 (0.3 recommended for more accurate evidences)` \
	--NS_sampling_efficiency 0.8 \
	`# the more live points the smoother the contours, empirical number, experiment (depends also on hardware available)` \
	--NS_n_live_points 1000 \
	`# run will finish/is converged if ln(Z_i) - ln(Z_j) <= NS_evidence_tolerance for i>j (0.1 is for evidences)` \
	--NS_evidence_tolerance 0.5