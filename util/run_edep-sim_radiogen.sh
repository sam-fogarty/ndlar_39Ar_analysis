physics_list=QGSP_BERT_LIV
#  0: 39Ar betas
#  1: 85Kr betas
#  2: 40K betas
#  3: 40K gammas
#  4: 60Co betas
#  5: 60Co gammas
#  6: 232Th betas
#  7: 232Th gammas
#  8: 232Th alphas
#  9: 238U betas
# 10: 238U gammas
# 11: 238U alphas
chosen_radiogen=0
DESCR=A
MACRO=A
N_LABEL=50k
N=50000
OUTDIR=/sdf/data/neutrino/sfogarty/edep_files/radiogen_uniform
if [ ${chosen_radiogen} -eq 0 ]; then
	DESCR=Ar39_betas_2x2_spatiallyuniform_${N_LABEL}
	MACRO=../../ndlar_39Ar_analysis/macros/2x2/Argon_39_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 1 ]; then
        DESCR=Kr85_betas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Krypton_85_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 2 ]; then
        DESCR=K40_betas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Potassium_40_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 3 ]; then
        DESCR=K40_gammas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Potassium_40_gammas_2x2.mac
elif [ ${chosen_radiogen} -eq 4 ]; then
        DESCR=Co60_betas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Cobalt_60_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 5 ]; then
        DESCR=Co60_gammas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Cobalt_60_gammas_2x2.mac
elif [ ${chosen_radiogen} -eq 6 ]; then
        DESCR=Th232_betas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Thorium_232_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 7 ]; then
        DESCR=Th232_gammas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Thorium_232_gammas_2x2.mac
elif [ ${chosen_radiogen} -eq 8 ]; then
        DESCR=Th232_alphas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Thorium_232_alphas_2x2.mac
elif [ ${chosen_radiogen} -eq 9 ]; then
        DESCR=U238_betas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Uranium_238_betas_2x2.mac
elif [ ${chosen_radiogen} -eq 10 ]; then
        DESCR=U238_gammas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Uranium_238_gammas_2x2.mac
elif [ ${chosen_radiogen} -eq 11 ]; then
        DESCR=U238_alphas_2x2_spatiallyuniform_${N_LABEL}
        MACRO=../../ndlar_39Ar_analysis/macros/2x2/Uranium_238_alphas_2x2.mac
fi
for i in {1..16}; do
	FILENAME=${DESCR}_${i}
	OUTFILE=${OUTDIR}/${FILENAME}.root
	edep-sim -e ${N} -p $physics_list -o ${OUTFILE} $MACRO
	rm -f /sdf/data/neutrino/sfogarty/edep_files/${FILENAME}.h5
	python3 dumpTree.py ${OUTFILE} ${OUTDIR}/${FILENAME}.h5
done
