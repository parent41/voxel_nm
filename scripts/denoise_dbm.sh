
ID=$1
input=$2
output=$3

# ID="sub-1000083_ses-2"
# input="../UKB/DBM_2mm"
# output="./maps_UKB_space_anlm_all"
# tmpdir="tmp"

tmpdir=$(mktemp -d)

metrics=('jacobians_rel' 'jacobians_abs')

for i in ${!metrics[@]}
do
    minc_anlm --short --mt 1 --beta 0.7 --clobber ${input}/${ID}_${metrics[i]}_2mm.mnc ${output}/${ID}_${metrics[i]}_2mm_anlm.mnc
done

rm -rf ${tmpdir}
