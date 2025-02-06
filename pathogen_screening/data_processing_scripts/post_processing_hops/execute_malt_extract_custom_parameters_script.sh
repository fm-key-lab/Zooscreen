#! /bin/bash

cd /ptmp/iclight/eager/screening/zooscreen
# general parameters
~/post_screening_malt_updates/post_eager_hops_postprocessing -O /ptmp/iclight/malt_reanalysis/zooscreen_final/zooscreen -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s SE -t 2 -F -S
# run allowance of no read dist cutoff for samples with many reads on a node
mkdir conditional_parameters
cp -r ancient conditional_parameters
cp -r default conditional_parameters
~/post_screening_malt_updates/custom_postprocessing.AMPS.r -r /ptmp/iclight/eager/screening/zooscreen/conditional_parameters -n ~/zooscreen/zooscreen_final_production_run.txt -t 2 -m def_anc -b 1000 -d 0.1 -c 0 -e 0.8 -a 0.5 -s SE -f

cd /ptmp/iclight/eager/screening/200423_tartu_data
# general parameters
~/post_screening_malt_updates/post_eager_hops_postprocessing -O /ptmp/iclight/malt_reanalysis/zooscreen_final/200423_tartu_data -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s PE -t 2 -F -S
# run allowance of no read dist cutoff for samples with many reads on a node
mkdir conditional_parameters
cp -r ancient conditional_parameters
cp -r default conditional_parameters
~/post_screening_malt_updates/custom_postprocessing.AMPS.r -r /ptmp/iclight/eager/screening/200423_tartu_data/conditional_parameters -n ~/zooscreen/zooscreen_final_production_run.txt -t 2 -m def_anc -b 1000 -d 0.1 -c 0 -e 0.8 -a 0.5 -s PE -f

cd /ptmp/iclight/eager/screening/tartu_screening_names_corrected
# general parameters
~/post_screening_malt_updates/post_eager_hops_postprocessing -O /ptmp/iclight/malt_reanalysis/zooscreen_final/tartu_screening_1 -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s PE -t 2 -F -S
# run allowance of no read dist cutoff for samples with many reads on a node
mkdir conditional_parameters
cp -r ancient conditional_parameters
cp -r default conditional_parameters
~/post_screening_malt_updates/custom_postprocessing.AMPS.r -r /ptmp/iclight/eager/screening/tartu_screening_names_corrected/conditional_parameters -n ~/zooscreen/zooscreen_final_production_run.txt -t 2 -m def_anc -b 1000 -d 0.1 -c 0 -e 0.8 -a 0.5 -s PE -f

