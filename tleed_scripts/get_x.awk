awk '
	/     1     1     1     1     1     2     2     2     2     2     2     2     2     2/ {
	lin=NR}
	{if (NR >lin+1 && NR <lin+16 && lin >0 || NR==lin+20 && lin>0)
#print 1     1     1     1     1     2     2     2     2     2     2     2     2     2
print}
' $*
	

