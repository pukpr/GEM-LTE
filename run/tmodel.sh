
NAME=`python3 lookup_index.py $1`

echo $NAME

if [ "$1" -gt 0 ]; then
  #./slh_avg_diff.awk ~/Downloads/rlr_monthly/data/$1.rlrdata >"$NAME" 
  # awk 'BEGIN{FS=";"}{print $1 " " $2}' <  ~/Downloads/rlr_monthly/data/$1.rlrdata >"$NAME" 
   python3 massage.py ~/Downloads/rlr_monthly/data/$1.rlrdata 2 >"$NAME"
  #awk '{getline x < "1.txt"; split(x, a); print $1, ($2 + a[2]) / 2}' 2.txt > "$NAME"
else
  echo "No filter"
fi


env QUADEXCLUDE=$2 ../obj/index_regress "$NAME" "$3 $4" | tail -n 4 | tee metrics.txt


python3 plot.py $1 $2 $3 $4 $5

