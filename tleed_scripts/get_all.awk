awk '
	/optimization parameters / {
	lin=NR}
	{if (NR==lin+1 && lin>0 ||NR >lin+2 && NR <lin+17 && lin >0 || NR==lin+21 && lin>0)
	print}
' $*
	

