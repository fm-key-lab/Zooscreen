#! /bin/bash

# run AFTER eager, to apply post processing maltextract to samples

cd /ptmp/iclight/eager/screening/zooscreen
~/post_screening_malt_updates/post_screening_test -O /ptmp/iclight/malt_reanalysis/zooscreen_final/zooscreen -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s SE -t 2 -F -S

cd /ptmp/iclight/eager/screening/200423_tartu_data
~/post_screening_malt_updates/post_screening_test -O /ptmp/iclight/malt_reanalysis/zooscreen_final/200423_tartu_data -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s PE -t 2 -F -S

cd /ptmp/iclight/eager/screening/tartu_screening_names_corrected
~/post_screening_malt_updates/post_screening_test -O /ptmp/iclight/malt_reanalysis/zooscreen_final/tartu_screening_1 -I . -n ~/zooscreen/zooscreen_final_production_run.txt -b 50 -i 0.6 -d 0.1 -e 0.8 -a 0.5 -s PE -t 2 -F -S
