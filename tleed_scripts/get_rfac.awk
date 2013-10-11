awk ' /rfactor =/{lin=NR}
	{if (NR ==lin+2 && lin > 0)
	print $1}
' $* #|sort -nr 

