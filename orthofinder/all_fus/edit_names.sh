FILE=$1

sed -i "s/>/>${FILE}_/g" $FILE
