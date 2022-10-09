qsubst () { 
	qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:$1 -l vf=$2,num_proc=$1 $3
} 
qsubsuper () {
	qsub -cwd -l vf=$2,num_proc=$1 -P P20Z10200N0206_super -binding linear:$1 -q st_supermem.q $3
}
qsubst 10 20G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_ALL.sh
qsubsuper 5 100G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_ALL.sh
qsubst 10 20G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_DNA.sh
qsubsuper 5 100G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_DNA.sh
qsubst 10 20G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_RNA.sh
qsubsuper 5 100G /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/ea_gam/output/2022-07-28/fit_RNA.sh
