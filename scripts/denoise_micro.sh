
ID=$1
input=$2
output=$3

# ID="sub-1000083_ses-2"
# input="../WMH_micro_spatial/maps_UKB_space"
# output="./maps_UKB_space_anlm_all"
# tmpdir="tmp"

tmpdir=$(mktemp -d)

# All but MD
micro_long=('dti_FA' 'NODDI_ICVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')

for i in ${!micro_long[@]}
do
    minc_anlm --short --mt 1 --beta 0.7 --clobber ${input}/${ID}_${micro_long[i]}_UKB.mnc ${output}/${ID}_${micro_long[i]}_UKB_anlm.mnc
done

# MD
mincmath -clobber -mult -const 1000 ${input}/${ID}_dti_MD_UKB.mnc ${tmpdir}/${ID}_MD_times1000.mnc
minc_anlm --short --mt 1 --beta 0.7 --clobber ${tmpdir}/${ID}_MD_times1000.mnc ${tmpdir}/${ID}_MD_times1000_anlm.mnc
mincmath -clobber -div ${tmpdir}/${ID}_MD_times1000_anlm.mnc -const 1000 ${tmpdir}/${ID}_MD_anlm.mnc

scp ${tmpdir}/${ID}_MD_anlm.mnc ${output}/${ID}_dti_MD_UKB_anlm.mnc

rm -rf ${tmpdir}
