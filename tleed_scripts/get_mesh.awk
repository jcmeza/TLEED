awk ' BEGIN {inc=1;
	mesh[inc]=1.0000}
	/rfactor =/{lin=NR
		rec=rec+1}
	{if (NR ==lin+2 && lin > 0)
		rfac[rec]=$1
		m[rec]=mesh[inc]}
	/mesh size/{inc=inc+1
		mesh[inc]=$6} 
	END {for (i=2;i<rec;i++)
		print rfac[i], m[i-1]
}
' $* #|sort -nr 


