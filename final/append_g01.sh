card=$1
part=$2
while [ "1" = "1" ]
do
        num=$(nvidia-smi|grep "^|    $card"|wc -l)
        if [ $num -eq 0 ];then
			/home/liuqiao/anaconda2/bin/python2.7	regression_gexp_trans.py $1 4 0.1 1 $part 
			exit
        fi
        echo "Waiting..."
		sleep 5s

    done
