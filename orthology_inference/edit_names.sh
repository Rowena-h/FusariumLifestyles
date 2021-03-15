FILES=$(ls *.faa)

for FILE in $FILES
do
	sed -i "s/>${FILE}_/>/g" ${FILE}
done
